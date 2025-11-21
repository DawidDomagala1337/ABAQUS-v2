"""
Compute radon exhalation from a cuboidal building element using the steady-state
radon diffusion model from Nazaroff and Nero (1988) as summarised in the
provided notes (Eqs. 1–10).

The model solves the steady-state diffusion equation with radioactive decay and
an internal radium source (Eq. 1) inside a cuboid of half-dimensions a, b, h
(thicknesses in the x, y, z directions). The pore-space concentration at all
surfaces is set to zero (Eqs. 2a–2c). A numerical Gauss–Seidel solver is used to
approximate the concentration field directly from Eq. 1; outward fluxes on the
six faces are then obtained with Fick's law (Eqs. 8a–8c). The mean exhalation
flux J_b is the area-weighted average flux over all faces (Eq. 9).

Inputs are limited to the physical parameters shown in the source equations:
- a, b, h: half-dimensions of the building element [m]
- D: radon diffusion coefficient in the material [m^2/s]
- lambda_: decay constant of radon [1/s]
- Q: radium content (emission factor) [Bq/kg]
- rho: bulk density of the material [kg/m^3]
- epsilon: porosity of the material [dimensionless]

Optional controls allow you to adjust the finite-difference grid resolution and
solver tolerances. The defaults are chosen to give a quick yet stable estimate.
"""

from __future__ import annotations

import argparse
import math
import sys
from typing import List, Tuple


Grid = List[List[List[float]]]

# Eq. (1): D * (d2C/dx2 + d2C/dy2 + d2C/dz2) - lambda * C = -lambda * Q * rho / epsilon
# Eq. (3): C(x, y, z) = Σ_{m=0..∞} Σ_{n=0..∞} A_mn cos(k_n x) cos(k_m y), which is the
#          separable cosine expansion shown in the provided figure. The code does not
#          evaluate the series directly; the Gauss–Seidel scheme below numerically solves
#          Eq. (1) for the same boundary-value problem instead of summing the coefficients
#          A_mn.
# Eq. (5): C(x, y, z) = Σ_{m=0..∞} Σ_{n=0..∞} 4C_∞ (-1)^{n+mn} d_mn^2 / l^2 · cos(k_n x) cos(k_m y)
#          · [ cosh( ωs h/d_mn ) - cosh( ωs z/d_mn ) ] / cosh( ωs h/d_mn ) is the fully expanded
#          cosine–hyperbolic solution for the zero-boundary cuboid. The implemented solver
#          targets the same solution space numerically rather than computing the closed-form
#          coefficients A_mn and d_mn explicitly.
# Eq. (6): l = sqrt(D / lambda) is the diffusion length used below for reporting; k_n = (n + 1/2)π / a
#          and k_m = (m + 1/2)π / b appear in the series definitions of k-space separable modes, and
#          1 / d_mn^2 = k_n^2 + k_m^2 + 1 / l^2 defines the combined decay term. The Gauss–Seidel
#          solver does not explicitly assemble these series coefficients, but the same diffusion
#          length l is calculated when summarising the run.
# Eq. (7): C(z) = C_inf * [1 - cosh(z / l) / cosh(h / l)] is the 1-D slab concentration profile; we
#          do not evaluate this closed form directly, but the z-directed dependence emerges from the
#          full 3-D Gauss–Seidel solution on the grid spanning [-h, h].
# Eq. (8a–c): J = -epsilon * D * (dC/dn) on each outward normal n; the finite-difference evaluation is
#          performed in face_fluxes() where each face uses the interior neighbour to approximate dC/dn.
# Eq. (10): J_b = lambda * Q * epsilon * f * sum_{m,n>=0} [4 a^2 b^2 / (ab + bh + ha)] * { (h/d_mn)
#           tanh(h/d_mn) [ (b/a) / ((n+1/2)^2 pi^2) + (a/b) / ((m+1/2)^2 pi^2) ] + (h/d_mn) tanh(h/d_mn)
#           (ab / d_mn^2) / ((n+1/2)^2 pi^2 (m+1/2)^2 pi^2) }. The implementation below does not sum this
#           series explicitly; instead, average_exhalation() evaluates the area-weighted mean flux (Eq. 9)
#           from the numerical face fluxes produced by the Gauss–Seidel solution of Eq. (1).


def linspace(start: float, stop: float, num: int) -> List[float]:
    if num < 2:
        return [start]
    step = (stop - start) / (num - 1)
    return [start + i * step for i in range(num)]


def build_grid(a: float, b: float, h: float, nodes: Tuple[int, int, int]) -> Tuple[List[float], List[float], List[float]]:
    """Return 1-D coordinate arrays covering [-a, a], [-b, b], [-h, h]."""
    nx, ny, nz = nodes
    x = linspace(-a, a, nx)
    y = linspace(-b, b, ny)
    z = linspace(-h, h, nz)
    return x, y, z


def zeros_grid(nx: int, ny: int, nz: int) -> Grid:
    return [[[0.0 for _ in range(nz)] for _ in range(ny)] for _ in range(nx)]


def solve_concentration(
    a: float,
    b: float,
    h: float,
    D: float,
    lambda_: float,
    Q: float,
    rho: float,
    epsilon: float,
    nodes: Tuple[int, int, int],
    max_iter: int,
    tol: float,
    progress_cb=None,
) -> Tuple[Grid, float, float, float]:
    """Solve Eq. (1) on a regular grid with zero concentration at the boundaries."""

    nx, ny, nz = nodes
    x, y, z = build_grid(a, b, h, nodes)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]

    source = lambda_ * Q * rho / epsilon
    denom = lambda_ + 2 * D * (1 / dx**2 + 1 / dy**2 + 1 / dz**2)

    C = zeros_grid(nx, ny, nz)

    for it in range(max_iter):
        max_delta = 0.0
        for i in range(1, nx - 1):
            for j in range(1, ny - 1):
                for k in range(1, nz - 1):
                    new_val = (
                        D
                        * (
                            (C[i + 1][j][k] + C[i - 1][j][k]) / dx**2
                            + (C[i][j + 1][k] + C[i][j - 1][k]) / dy**2
                            + (C[i][j][k + 1] + C[i][j][k - 1]) / dz**2
                        )
                        + source
                    ) / denom
                    delta = abs(new_val - C[i][j][k])
                    if delta > max_delta:
                        max_delta = delta
                    C[i][j][k] = new_val
        if progress_cb is not None:
            progress_cb(it + 1, max_iter)
        if max_delta < tol:
            break

    return C, dx, dy, dz


def face_flux_mean(face_values: List[List[float]]) -> float:
    total = 0.0
    count = 0
    for row in face_values:
        for v in row:
            total += v
            count += 1
    return total / count


def face_fluxes(C: Grid, dx: float, dy: float, dz: float, D: float, epsilon: float) -> Tuple[float, float, float, float, float, float]:
    """Compute outward flux on each face using Eqs. (8a)–(8c)."""

    flux_x_neg = [[epsilon * D * C[1][j][k] / dx for k in range(len(C[0][0]))] for j in range(len(C[0]))]
    flux_x_pos = [[epsilon * D * C[-2][j][k] / dx for k in range(len(C[0][0]))] for j in range(len(C[0]))]

    flux_y_neg = [[epsilon * D * C[i][1][k] / dy for k in range(len(C[0][0]))] for i in range(len(C))]
    flux_y_pos = [[epsilon * D * C[i][-2][k] / dy for k in range(len(C[0][0]))] for i in range(len(C))]

    flux_z_neg = [[epsilon * D * C[i][j][1] / dz for j in range(len(C[0]))] for i in range(len(C))]
    flux_z_pos = [[epsilon * D * C[i][j][-2] / dz for j in range(len(C[0]))] for i in range(len(C))]

    return (
        face_flux_mean(flux_x_neg),
        face_flux_mean(flux_x_pos),
        face_flux_mean(flux_y_neg),
        face_flux_mean(flux_y_pos),
        face_flux_mean(flux_z_neg),
        face_flux_mean(flux_z_pos),
    )


def average_exhalation(
    fluxes: Tuple[float, float, float, float, float, float],
    a: float,
    b: float,
    h: float,
) -> float:
    """Area-weighted mean flux across the six faces (Eq. 9)."""

    # Eq. (9): J_b = [∬(J_x^- + J_x^+) dy dz + ∬(J_y^- + J_y^+) dz dx
    #              + ∬(J_z^- + J_z^+) dx dy] / [8(ab + bh + ha)]
    # The face areas below (area_x_faces etc.) correspond to the denominator, while the
    # weighted_flux sums the opposing-face integrals from the discrete mean fluxes.
    # Eq. (10) gives an equivalent series expression for J_b using modal coefficients
    # d_mn; the numerical average here is the practical evaluation of that series via
    # the finite-difference fluxes.

    fx_neg, fx_pos, fy_neg, fy_pos, fz_neg, fz_pos = fluxes

    area_x_faces = 2 * (2 * b) * (2 * h)  # two faces normal to x
    area_y_faces = 2 * (2 * a) * (2 * h)
    area_z_faces = 2 * (2 * a) * (2 * b)

    total_area = area_x_faces + area_y_faces + area_z_faces

    weighted_flux = (
        fx_neg * area_x_faces / 2
        + fx_pos * area_x_faces / 2
        + fy_neg * area_y_faces / 2
        + fy_pos * area_y_faces / 2
        + fz_neg * area_z_faces / 2
        + fz_pos * area_z_faces / 2
    )

    return weighted_flux / total_area


def run_model(
    a: float,
    b: float,
    h: float,
    D: float,
    lambda_: float,
    Q: float,
    rho: float,
    epsilon: float,
    nodes: Tuple[int, int, int],
    max_iter: int,
    tol: float,
    progress_cb=None,
):
    """Run the diffusion solver and return summary metrics for UI/CLI callers."""

    for name, value in [
        ("a", a),
        ("b", b),
        ("h", h),
        ("D", D),
        ("lambda_", lambda_),
        ("Q", Q),
        ("rho", rho),
        ("max_iter", max_iter),
        ("tol", tol),
    ]:
        if value <= 0:
            raise ValueError(f"{name} must be greater than zero")

    if epsilon <= 0:
        raise ValueError("epsilon must be greater than zero")
    if epsilon < 0.005 or epsilon > 1:
        raise ValueError("epsilon must be between 0.005 and 1")

    if any(n < 3 for n in nodes):
        raise ValueError("Each node count must be at least 3 to include boundaries.")
    if any(n > 50 for n in nodes):
        raise ValueError("Each node count must be 50 or fewer to keep the solver stable in the UI.")

    C, dx, dy, dz = solve_concentration(
        a,
        b,
        h,
        D,
        lambda_,
        Q,
        rho,
        epsilon,
        nodes,
        max_iter,
        tol,
        progress_cb=progress_cb,
    )

    fluxes = face_fluxes(C, dx, dy, dz, D, epsilon)
    Jb = average_exhalation(fluxes, a, b, h)

    # Eq. (6): diffusion length l = sqrt(D / lambda)
    diffusion_length = math.sqrt(D / lambda_)

    return {
        "fluxes": fluxes,
        "Jb": Jb,
        "diffusion_length": diffusion_length,
        "grid_spacing": (dx, dy, dz),
        "nodes": nodes,
    }


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Predict radon exhalation from a cuboidal building element using Eqs. 1–10.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--example",
        action="store_true",
        help=(
        "Use the published example parameters (a=0.3 m, b=0.15 m, h=0.3 m, D=8e-7 m^2/s, "
        "lambda=2.1e-6 1/s, Q=20 Bq/kg, rho=2400 kg/m^3, epsilon=0.25, nodes=5 5 5)."
        ),
    )
    parser.add_argument("a", type=float, nargs="?", help="Half-thickness along x [m]")
    parser.add_argument("b", type=float, nargs="?", help="Half-breadth along y [m]")
    parser.add_argument("h", type=float, nargs="?", help="Half-height along z [m]")
    parser.add_argument("D", type=float, nargs="?", help="Radon diffusion coefficient [m^2/s]")
    parser.add_argument("lambda_", type=float, nargs="?", help="Radon decay constant [1/s]")
    parser.add_argument("Q", type=float, nargs="?", help="Radium content (emanation factor) [Bq/kg]")
    parser.add_argument("rho", type=float, nargs="?", help="Material bulk density [kg/m^3]")
    parser.add_argument("epsilon", type=float, nargs="?", help="Porosity (pore fraction)")
    parser.add_argument(
        "--nodes",
        type=int,
        nargs=3,
        metavar=("NX", "NY", "NZ"),
        default=None,
        help="Grid nodes in x, y, z (>=3).",
    )
    parser.add_argument("--tol", type=float, default=1e-6, help="Convergence tolerance for Gauss–Seidel.")
    parser.add_argument("--max-iter", type=int, default=5000, help="Maximum Gauss–Seidel iterations.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])

    example_params = dict(
        a=0.3,
        b=0.15,
        h=0.3,
        D=8e-7,
        lambda_=2.1e-6,
        Q=20.0,
        rho=2400.0,
        epsilon=0.25,
        nodes=(5, 5, 5),
    )

    # Either consume user-provided positional arguments or fall back to the
    # documented example when --example is given. A clear error is raised if
    # neither path supplies the full parameter set.
    if args.example:
        a = args.a if args.a is not None else example_params["a"]
        b = args.b if args.b is not None else example_params["b"]
        h = args.h if args.h is not None else example_params["h"]
        D = args.D if args.D is not None else example_params["D"]
        lambda_ = args.lambda_ if args.lambda_ is not None else example_params["lambda_"]
        Q = args.Q if args.Q is not None else example_params["Q"]
        rho = args.rho if args.rho is not None else example_params["rho"]
        epsilon = args.epsilon if args.epsilon is not None else example_params["epsilon"]
        nodes = tuple(args.nodes) if args.nodes is not None else example_params["nodes"]
    else:
        missing = [
            name
            for name, val in [
                ("a", args.a),
                ("b", args.b),
                ("h", args.h),
                ("D", args.D),
                ("lambda_", args.lambda_),
                ("Q", args.Q),
                ("rho", args.rho),
                ("epsilon", args.epsilon),
            ]
            if val is None
        ]
        if missing:
            raise SystemExit(
                "Missing parameters: "
                + ", ".join(missing)
                + ". Either provide all positional arguments or use --example."
            )
        a, b, h, D, lambda_, Q, rho, epsilon = (
            args.a,
            args.b,
            args.h,
            args.D,
            args.lambda_,
            args.Q,
            args.rho,
            args.epsilon,
        )
        nodes = tuple(args.nodes) if args.nodes is not None else (5, 5, 5)

    try:
        result = run_model(
            a,
            b,
            h,
            D,
            lambda_,
            Q,
            rho,
            epsilon,
            nodes,
            args.max_iter,
            args.tol,
        )
    except ValueError as exc:
        raise SystemExit(str(exc))

    fluxes = result["fluxes"]
    Jb = result["Jb"]
    diffusion_length = result["diffusion_length"]

    print("Radon exhalation model (Eqs. 1–10)")
    print(f"Diffusion length l = {diffusion_length:.4e} m")
    print("Mean outward fluxes on faces [Bq m^-2 s^-1]:")
    labels = ["x-", "x+", "y-", "y+", "z-", "z+"]
    for label, flux in zip(labels, fluxes):
        print(f"  {label}: {flux:.4e}")
    print(f"Area-weighted mean exhalation J_b = {Jb:.4e} Bq m^-2 s^-1")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

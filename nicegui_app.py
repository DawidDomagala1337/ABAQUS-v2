"""
NiceGUI front-end for the radon exhalation solver in ``radon_exhalation.py``.

Launch with ``python nicegui_app.py`` after installing NiceGUI, enter the geometry
and material parameters, and view the diffusion length plus face and mean fluxes
computed by the shared solver.
"""

from __future__ import annotations

import importlib.util
from typing import Tuple

from radon_exhalation import run_model


def _load_nicegui():
    spec = importlib.util.find_spec("nicegui")
    if spec is None:
        raise SystemExit(
            "NiceGUI is not installed. Install it with 'pip install nicegui' to use the web UI."
        )
    from nicegui import ui

    return ui


def launch_ui(example_params: dict[str, float | Tuple[int, int, int]]):
    ui = _load_nicegui()

    def number_input(label: str, default: float, digits: int = 6, *, min_value=None, max_value=None):
        return ui.number(label, value=default, format=f"%.{digits}g", min=min_value, max=max_value)

    def integer_input(label: str, default: int):
        return ui.number(label, value=default, format="%d", step=1, min=3)

    with ui.header().classes("items-center justify-between"):
        ui.label("Radon exhalation (Eqs. 1–10) – NiceGUI front-end")
        ui.label("Fill parameters and run the solver")

    with ui.row().classes("w-full"):
        with ui.card().classes("w-1/2"):
            ui.label("Geometry [m]").classes("text-bold")
            a_input = number_input("a (half-thickness, x)", example_params["a"])
            b_input = number_input("b (half-breadth, y)", example_params["b"])
            h_input = number_input("h (half-height, z)", example_params["h"])

            ui.separator()
            ui.label("Material properties").classes("text-bold")
            D_input = number_input(
                "D (diffusion coefficient) [m²/s]", example_params["D"], digits=4
            )
            lambda_input = number_input(
                "λ (decay constant) [1/s]", example_params["lambda_"], digits=4
            )
            Q_input = number_input(
                "Q (radium content) [Bq/kg]", example_params["Q"], digits=4
            )
            rho_input = number_input(
                "ρ (bulk density) [kg/m³]", example_params["rho"], digits=4
            )
            epsilon_input = number_input(
                "ε (porosity)", example_params["epsilon"], digits=3, min_value=0.005, max_value=1
            )

        with ui.card().classes("w-1/2"):
            ui.label("Solver controls").classes("text-bold")
            nx_input = integer_input("Nodes in x", example_params["nodes"][0])
            ny_input = integer_input("Nodes in y", example_params["nodes"][1])
            nz_input = integer_input("Nodes in z", example_params["nodes"][2])

            max_iter_input = integer_input("Max iterations", example_params["max_iter"])
            tol_input = number_input("Tolerance", example_params["tol"], digits=2)

            result_area = ui.column().classes("gap-2 mt-4")

    async def fill_defaults():
        a_input.value = example_params["a"]
        b_input.value = example_params["b"]
        h_input.value = example_params["h"]
        D_input.value = example_params["D"]
        lambda_input.value = example_params["lambda_"]
        Q_input.value = example_params["Q"]
        rho_input.value = example_params["rho"]
        epsilon_input.value = example_params["epsilon"]
        nx_input.value, ny_input.value, nz_input.value = example_params["nodes"]
        max_iter_input.value = example_params["max_iter"]
        tol_input.value = example_params["tol"]

    async def run_solver():
        result_area.clear()
        try:
            nodes: Tuple[int, int, int] = (
                int(nx_input.value),
                int(ny_input.value),
                int(nz_input.value),
            )
            output = run_model(
                float(a_input.value),
                float(b_input.value),
                float(h_input.value),
                float(D_input.value),
                float(lambda_input.value),
                float(Q_input.value),
                float(rho_input.value),
                float(epsilon_input.value),
                nodes,
                int(max_iter_input.value),
                float(tol_input.value),
            )
        except Exception as exc:
            ui.notify(str(exc), color="negative")
            return

        flux_labels = ["x-", "x+", "y-", "y+", "z-", "z+"]
        fluxes = output["fluxes"]

        with result_area:
            ui.label("Results").classes("text-lg text-bold")
            ui.label(f"Diffusion length l = {output['diffusion_length']:.4e} m")
            ui.label(
                f"Area-weighted mean exhalation J_b = {output['Jb']:.4e} Bq m^-2 s^-1"
            )
            flux_rows = [
                {"face": name, "flux": f"{value:.4e}"}
                for name, value in zip(flux_labels, fluxes)
            ]
            ui.table(
                columns=[
                    {"name": "face", "label": "Face", "field": "face"},
                    {
                        "name": "flux",
                        "label": "Mean outward flux [Bq m^-2 s^-1]",
                        "field": "flux",
                    },
                ],
                rows=flux_rows,
            )

    with ui.footer().classes("items-center justify-between"):
        ui.button("Use example parameters", on_click=fill_defaults)
        ui.button("Run solver", color="primary", on_click=run_solver)

    ui.run(title="Radon exhalation GUI", reload=False)


def main():
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
        max_iter=5000,
        tol=1e-6,
    )
    launch_ui(example_params)


if __name__ == "__main__":
    main()

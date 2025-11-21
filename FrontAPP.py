"""NiceGUI front-end for the radon exhalation solver.

Launch with ``python FrontAPP.py`` and open http://127.0.0.1:8001 to enter
model parameters, run the calculation, and view diffusion length, face fluxes,
and area-weighted mean exhalation with live progress updates.
"""

from __future__ import annotations

import asyncio
import math
from pathlib import Path
from urllib.parse import quote

from nicegui import ui

import radon_exhalation


LAMBDA_CONST = 2.1e-6
DIAGRAM_DATA_URI = "data:image/svg+xml;utf8," + quote(Path("static_cuboid.svg").read_text())


def parse_float(value: str, name: str) -> float:
    """Parse a float value from a string and raise a helpful error."""
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(f"{name} must be a number") from exc
    if parsed <= 0:
        raise ValueError(f"{name} must be greater than zero")
    return parsed


def parse_int(
    value: str, name: str, *, min_value: int = 1, max_value: int | None = None
) -> int:
    """Parse an int value ensuring grid sizes stay valid for the solver."""
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if parsed < min_value:
        raise ValueError(f"{name} must be at least {min_value}")
    if max_value is not None and parsed > max_value:
        raise ValueError(f"{name} must be at most {max_value}")
    return parsed


@ui.page("/")
def main_page() -> None:
    """Render the input form and results display for the solver."""

    ui.markdown("## Radon exhalation model inputs")

    with ui.row().classes("flex-wrap gap-6 items-start"):
        with ui.column().classes("gap-2 max-w-4xl"):
            # Input fields grouped for clarity
            with ui.row().classes("flex-wrap gap-4"):
                a = ui.input("a (half-thickness x) [m]", value="0.5").props(
                    "type=number step=any min=1e-12"
                )
                b = ui.input("b (half-breadth y) [m]", value="0.5").props(
                    "type=number step=any min=1e-12"
                )
                h = ui.input("h (half-height z) [m]", value="0.5").props(
                    "type=number step=any min=1e-12"
                )
            with ui.row().classes("flex-wrap gap-4"):
                D = ui.input("D (diffusion coeff.) [m^2/s]", value="8e-7").props(
                    "type=number step=any min=1e-12"
                )
                lambda_display = (
                    ui.input("lambda (decay const.) [1/s]", value=f"{LAMBDA_CONST}")
                    .props("type=number step=any readonly")
                    .tooltip("Physical constant; fixed in the model")
                )
                Q = ui.input("Q (radium content) [Bq/kg]", value="20").props(
                    "type=number step=any min=1e-12"
                )
            with ui.row().classes("flex-wrap gap-4"):
                rho = ui.input("rho (bulk density) [kg/m^3]", value="2400").props(
                    "type=number step=any min=1e-12"
                )
                epsilon = ui.input("epsilon (porosity)", value="0.25").props(
                    "type=number step=any min=0.005 max=1"
                )
            with ui.row().classes("flex-wrap gap-4 items-end mt-2"):
                with ui.column().classes("min-w-56"):
                    ui.label("Grid nodes").classes("text-sm text-gray-700")
                    ui.label("3–50 per dimension").classes("text-xs text-gray-500")
                node_props = "type=number step=1 min=3 max=50"
                node_width = "min-w-44"
                nx = ui.input("nodes in x", value="5").props(node_props).classes(
                    node_width
                )
                ny = ui.input("nodes in y", value="5").props(node_props).classes(
                    node_width
                )
                nz = ui.input("nodes in z", value="5").props(node_props).classes(
                    node_width
                )
            with ui.row().classes("flex-wrap gap-4"):
                max_iter = ui.input("max iterations", value="5000").props(
                    "type=number step=1 min=1"
                )
                tol = ui.input("tolerance", value="1e-6").props(
                    "type=number step=any min=1e-12"
                )

        with ui.card().classes("items-center gap-3 p-3 min-w-[320px]"):
            ui.image(DIAGRAM_DATA_URI).props("fit=contain").style(
                "max-width: 360px;"
            )
            ui.label(
                "Schematic diagram of a cuboidal shaped building material and the coordinate system used for the model development."
            ).classes("text-sm text-gray-700 text-center")

    message = ui.label().classes("text-red-700")

    with ui.dialog() as progress_dialog, ui.card():
        ui.markdown("### Solving...")
        progress_bar = ui.linear_progress(value=0).props("striped indeterminate=false").classes("w-64")
        progress_label = ui.label("0% complete")

    async def on_calculate() -> None:
        message.text = ""
        result_area.content = ""
        try:
            parsed_nodes = (
                parse_int(nx.value, "nodes x", min_value=3, max_value=50),
                parse_int(ny.value, "nodes y", min_value=3, max_value=50),
                parse_int(nz.value, "nodes z", min_value=3, max_value=50),
            )
            epsilon_val = parse_float(epsilon.value, "epsilon")
            if epsilon_val < 0.005 or epsilon_val > 1:
                raise ValueError("epsilon must be between 0.005 and 1")

            params = dict(
                a=parse_float(a.value, "a"),
                b=parse_float(b.value, "b"),
                h=parse_float(h.value, "h"),
                D=parse_float(D.value, "D"),
                lambda_=LAMBDA_CONST,
                Q=parse_float(Q.value, "Q"),
                rho=parse_float(rho.value, "rho"),
                epsilon=epsilon_val,
                nodes=parsed_nodes,
                max_iter=parse_int(max_iter.value, "max iterations"),
                tol=parse_float(tol.value, "tolerance"),
            )
        except Exception as exc:  # noqa: BLE001 - surface user input errors in UI
            message.text = str(exc)
            return

        progress_bar.value = 0
        progress_label.text = "0% complete"
        progress_dialog.open()
        queue: asyncio.Queue[float] = asyncio.Queue()
        loop = asyncio.get_running_loop()

        def quantize_progress(fraction: float, current: float) -> float:
            """Clamp progress to non-decreasing 5% steps."""

            stepped = math.floor(fraction / 0.05 + 0.5) * 0.05
            return min(1.0, max(current, stepped))

        def progress_cb(done: int, total: int) -> None:
            raw_fraction = max(0.0, min(1.0, done / total))
            stepped = quantize_progress(raw_fraction, 0.0)
            asyncio.run_coroutine_threadsafe(queue.put(stepped), loop)

        async def run_solver():
            """Call run_model with progress support when available."""

            try:
                return await asyncio.to_thread(
                    radon_exhalation.run_model, progress_cb=progress_cb, **params
                )
            except TypeError as exc:
                # Backward-compatibility: fall back if run_model lacks progress_cb
                if "progress_cb" not in str(exc):
                    raise
                return await asyncio.to_thread(radon_exhalation.run_model, **params)

        solver_task = asyncio.create_task(run_solver())

        display_frac = 0.0
        last_event = loop.time()
        try:
            while not solver_task.done():
                try:
                    frac = await asyncio.wait_for(queue.get(), timeout=0.2)
                    display_frac = quantize_progress(frac, display_frac)
                    last_event = loop.time()
                except asyncio.TimeoutError:
                    # If the solver does not emit progress (older versions or long steps),
                    # gently advance the bar so users see activity.
                    if loop.time() - last_event > 0.5 and display_frac < 0.95:
                        display_frac = quantize_progress(display_frac + 0.05, display_frac)
                        last_event = loop.time()
                progress_bar.value = display_frac
                progress_label.text = f"{display_frac * 100:.0f}% complete"
            result = await solver_task
        except Exception as exc:  # noqa: BLE001 - surface user input errors in UI
            message.text = str(exc)
            progress_dialog.close()
            return

        progress_bar.value = 1
        progress_label.text = "100% complete"
        progress_dialog.close()

        fluxes = result["fluxes"]
        diffusion_length = result["diffusion_length"]
        Jb = result["Jb"]

        lines = [
            f"Diffusion length l = {diffusion_length:.4e} m",
            "Mean outward fluxes on faces [Bq m^-2 s^-1]:",
            f"  x-: {fluxes[0]:.4e}",
            f"  x+: {fluxes[1]:.4e}",
            f"  y-: {fluxes[2]:.4e}",
            f"  y+: {fluxes[3]:.4e}",
            f"  z-: {fluxes[4]:.4e}",
            f"  z+: {fluxes[5]:.4e}",
            f"Area-weighted mean exhalation J_b = {Jb:.4e} Bq m^-2 s^-1",
        ]
        result_area.content = f"```\n{'\n'.join(lines)}\n```"

    ui.button("Calculate", on_click=on_calculate).props("color=primary no-caps")
    with ui.card().classes("w-full mt-4"):
        ui.markdown("### Results")
        result_area = ui.markdown("Enter parameters and press **Calculate**.").classes("mt-2")
    ui.label(
        "Server runs on 127.0.0.1:8001. Keep the process active while using the page."
    ).classes("text-sm text-gray-500 mt-1")

ui.run(host="127.0.0.1", port=8001)

from pathlib import Path
from time import perf_counter

from config import read_parameters
from mesh import generate_mesh
from solver import Solver
from visualization import (
    plot_stress,
    plot_displacement,
    plot_stress_vonMises,
    save_results,
)


def main():
    config = read_parameters(Path(""))
    output = Path("")
    output.mkdir(parents=True, exist_ok=True)
    times = {}

    # Generate mesh
    start = perf_counter()
    mesh = generate_mesh(config)
    times["mesh"] = perf_counter() - start

    # Solve system
    solver = Solver(mesh, config)
    start = perf_counter()
    solver.solve()
    times["solve"] = perf_counter() - start

    # Compute stresses
    start = perf_counter()
    solver.compute_stress()
    solver.compute_von_mises()
    times["stress"] = perf_counter() - start

    # Disply results
    save_results(solver, output, times)
    plot_displacement( solver, component="magnitude", scale=1000, show=False, save_path=output / "displacement.png")
    plot_stress_vonMises(solver, show=False, save_path=output / "von_mises.png")

    return solver


if __name__ == "__main__":
    main()
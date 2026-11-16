from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

from elements import T6Element


def plot_stress(solver, component="xx", show=True, save_path=None):
    """Plot element-centroid stresses as a discontinuous field."""
    if solver.stresses is None:
        raise ValueError("Compute stresses before plotting.")

    components = {"xx": 0, "yy": 1, "xy": 2}

    if component not in components:
        raise ValueError("Component must be 'xx', 'yy', or 'xy'.")

    mesh = solver.mesh
    values = solver.stresses[:, components[component]]

    # T6 elements use their three corner nodes for visualization.
    triangles = mesh.elements[:, :3]

    fig, ax = plt.subplots(figsize=(10, 5))

    plot = ax.tripcolor(
        mesh.nodes[:, 0],
        mesh.nodes[:, 1],
        triangles,
        facecolors=values,
        shading="flat",
        cmap="viridis"
    )

    fig.colorbar(plot, ax=ax, label=f"Stress σ_{component} (Pa)")

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"{mesh.element_type} - Stress σ_{component}")

    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300)

    if show:
        plt.show()

    plt.close(fig)


def plot_displacement(solver, component="magnitude", scale=1,
                      show=True, save_path=None):
    """Plot the deformed mesh colored by nodal displacement."""
    if solver.u is None:
        raise ValueError("Solve the system before plotting displacements.")

    mesh = solver.mesh
    u = solver.u.reshape(-1, 2)

    components = {
        "x": u[:, 0],
        "y": u[:, 1],
        "magnitude": np.linalg.norm(u, axis=1)
    }

    if component not in components:
        raise ValueError("Component must be 'x', 'y', or 'magnitude'.")

    # Deformed node coordinates
    deformed_nodes = mesh.nodes + scale * u

    # T6 elements use corner nodes for visualization
    triangles = mesh.elements[:, :3]

    fig, ax = plt.subplots(figsize=(10, 5))

    plot = ax.tripcolor(
        deformed_nodes[:, 0],
        deformed_nodes[:, 1],
        triangles,
        components[component],
        shading="gouraud",
        cmap="jet"
    )

    fig.colorbar(plot, ax=ax, label=f"Displacement {component} (m)")

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(
        f"{mesh.element_type} - Displacement {component} "
        f"(scale = {scale:g})"
    )

    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300)

    if show:
        plt.show()

    plt.close(fig)


def plot_stress_vonMises(solver, show=True, save_path=None):
    """Plot the maximum sampled von Mises stress per element."""
    if solver.von_mises is None:
        raise ValueError("Compute von Mises stresses before plotting.")

    mesh = solver.mesh
    triangles = mesh.elements[:, :3]

    fig, ax = plt.subplots(figsize=(10, 5))

    plot = ax.tripcolor(
        mesh.nodes[:, 0],
        mesh.nodes[:, 1],
        triangles,
        facecolors=solver.von_mises,
        shading="flat",
        cmap="jet"
    )

    fig.colorbar(plot, ax=ax, label="Von Mises stress (Pa)")

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"{mesh.element_type} - Maximum von Mises stress")

    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300)

    if show:
        plt.show()

    plt.close(fig)


def save_results(solver, experiment_path, times):
    """Save mesh statistics, maximum responses and computation times."""
    if solver.u is None or solver.von_mises is None:
        raise ValueError("Solve and compute von Mises stresses before saving results.")

    mesh = solver.mesh

    displacement = solver.u.reshape(-1, 2)
    max_displacement = np.max(np.linalg.norm(displacement, axis=1))

    max_von_mises = np.max(solver.von_mises)

    total_time = sum(times.values())

    results = (
        "FINITE ELEMENT ANALYSIS RESULTS\n"
        "========================================\n\n"
        f"Element type:       {mesh.element_type}\n"
        f"Number of nodes:    {mesh.number_of_nodes}\n"
        f"Number of elements: {mesh.number_of_elements}\n"
        f"Number of DOFs:     {2 * mesh.number_of_nodes}\n\n"
        "MAXIMUM RESPONSE\n"
        "----------------------------------------\n"
        f"Maximum displacement: {max_displacement:.6e} m\n"
        f"Maximum von Mises:     {max_von_mises:.6e} Pa\n\n"
        "COMPUTATIONAL TIME\n"
        "----------------------------------------\n"
        f"Mesh generation:    {times['mesh']:.6f} s\n"
        f"Assembly and solve: {times['solve']:.6f} s\n"
        f"Stress computation: {times['stress']:.6f} s\n"
        f"Total:              {total_time:.6f} s\n"
    )

    experiment_path = Path(experiment_path)
    experiment_path.mkdir(parents=True, exist_ok=True)

    output_path = experiment_path / "results.txt"
    output_path.write_text(results, encoding="utf-8")

    print(results)
    return output_path


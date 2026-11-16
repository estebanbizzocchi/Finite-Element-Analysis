# Finite Element Analysis of a Plate with Circular Holes

A Python implementation of two-dimensional, linear-elastic finite element analysis under plane stress. The project studies a rectangular plate subjected to equal and opposite horizontal tensile tractions, with optional circular holes, and compares linear (T3) and quadratic (T6) triangular elements.

The repository includes the solver, saved numerical experiments, CSV summaries, plotting scripts, and a project report in [`Report.pdf`](Report.pdf).

## Features

- Gmsh geometry construction and triangular mesh generation.
- Three-node T3 and six-node T6 elements; six-point quadrature for T6 stiffness integration.
- Sparse global stiffness assembly and solution using SciPy.
- Displacement, centroid stress, and sampled von Mises stress calculations.
- Saved displacement and stress plots, mesh statistics, and computation times.
- Archived mesh-refinement studies for several hole radii and a plate without holes.

## Physical model

The plate occupies `0 <= x <= length`, `0 <= y <= width`, with constant thickness. The material is homogeneous, isotropic, and linearly elastic, with small displacements and a plane-stress constitutive law.

Horizontal tractions of `-traction` and `+traction` act on the left and right edges. The other external edges and hole boundaries are traction-free. Three displacement degrees of freedom remove rigid-body motion: both components at the bottom-left corner and the vertical component at the bottom-right corner.

Use consistent SI units: metres for geometry, thickness, and mesh size; pascals for Young's modulus, traction, and stresses. Displacements are reported in metres.

## Installation

The supplied environment uses Python **3.14.3**. `requirements.txt` pins the five direct dependencies to the versions installed in that environment. Other Python/package combinations have not been verified. Gmsh includes a native library, so its package must support your operating system and Python environment.

From the project directory, create a new virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, use `py -m venv .venv`, then activate it in PowerShell with `.venv\Scripts\Activate.ps1`.

Do not upload or copy the existing `.venv` folder between computers. Recreate it from the dependency file.

## Run an analysis

Run from the directory containing `config.yaml`:

```bash
python main.py
```

The program reads `config.yaml` from the **current working directory** and writes these files into that same directory:

| File | Content |
| --- | --- |
| `results.txt` | Mesh size statistics, maximum displacement, maximum sampled von Mises stress, and timings |
| `displacement.png` | Displacement magnitude on a deformed mesh, with deformation magnified by 1000 |
| `von_mises.png` | Maximum sampled von Mises stress per element |

Repeated runs overwrite these three files. The displacement scale changes the drawing only, not the reported displacement. Timing totals include mesh generation, assembly/solution, and stress evaluation, but exclude plotting and file writing.

For a noninteractive run on macOS/Linux:

```bash
MPLBACKEND=Agg python main.py
```

## Configuration

The supplied example is:

```yaml
geometry:
  length: 10.0
  width: 4.0
  thickness: 0.01
  holes:
    - center: [5.0, 2.0]
      radius: 1.0

material:
  young_modulus: 210.0e9
  poisson_ratio: 0.3

traction: 1.0e6

mesh:
  element_type: T6
  element_size: 0.2
```

Set `element_type` to `T3` or `T6`. Smaller `element_size` values generally produce finer, more expensive meshes. Use `holes: []` for a solid plate; additional entries define additional circular cutouts.

Geometry dimensions, thickness, Young's modulus, and mesh size must be positive and finite. Poisson's ratio must satisfy `-1 < poisson_ratio < 0.5`. Hole radii must be positive and centers must contain two finite coordinates. The configuration checks do not enforce hole containment or separation: use non-overlapping holes fully inside the plate for the intended model.

## Project structure

```text
main.py                 Analysis entry point
config.py               YAML loading and parameter validation
config.yaml             Default analysis parameters
elements.py             T3 and T6 element formulations
mesh.py                 Gmsh geometry, mesh extraction, and constraints
solver.py               Assembly, solution, and stress recovery
visualization.py        Plots and text result export
plots/                  Scripts for study figures
experiments/            Archived configurations, results, and stress images
results/                CSV summaries by hole radius
Report.pdf              Project report
requirements.txt        Direct dependency versions
```

## Saved studies and plots

`experiments/` contains a no-hole study and studies with hole radii of 0.1, 0.2, 0.5, 1.0, 1.2, 1.5, and 1.8. Subdirectories such as `T3_0.2` identify element type and mesh size. The archived stress images are named `vM.png`; current runs produce `von_mises.png`.

The `results/results<RADIUS>.csv` files contain element type, element size, maximum von Mises stress, and total time. The current entry point does **not** automatically regenerate these CSV files or execute the complete parameter sweep.

To generate the study figures from the stored arrays in the plotting scripts:

```bash
python plots/plot_SCF.py
python plots/plot_convergence.py
```

These scripts use **hard-coded arrays**, not the CSV files or the output of the latest simulation. `plot_SCF.py` creates `stress_concentration.pdf`. `plot_convergence.py` currently selects radius 0.5 and creates `stress_vs_mesh_0.5.pdf` and `accuracy_vs_time_0.5.pdf`. They save into the current working directory and display figures; close each figure to continue, or use `MPLBACKEND=Agg` for noninteractive execution.

For a new convergence study, vary `element_type` and `element_size`, preserve each run's configuration and outputs, and update the plotting data together. Changing only the `radius` variable in `plot_convergence.py` changes the labels, not the underlying data.

## Interpretation and limitations

- Stress components are evaluated at element centroids. T3 stresses are constant within each element. For T6, the reported von Mises value is the maximum sampled at the centroid and six element nodes, not an exact continuous maximum.
- Plots use the three corner nodes of each triangle, including T6 elements; they do not display the full curved quadratic geometry.
- The stress-concentration plot uses maximum von Mises stress divided by the applied nominal traction. Its reference line at 3 represents the ideal infinite-plate, circular-hole limit; it is not an exact reference for every finite-width case.
- Convergence errors in the plotting script are measured against the finest stored T6 result, not an analytical solution.
- This implementation does not model plasticity, large deformation, dynamics, or three-dimensional effects.
- There is no automated test suite or command-line configuration interface in the supplied project. Archived results are study records rather than automated regression tests.


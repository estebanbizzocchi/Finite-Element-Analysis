import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from elements import T3Element, T6Element


class Solver:
    """Assemble and solve the 2D plane-stress elasticity problem."""

    def __init__(self, mesh, config):
        self.mesh = mesh
        self.config = config

        E = config.young_modulus
        nu = config.poisson_ratio

        # Plane-stress constitutive matrix
        self.D = E / (1 - nu**2) * np.array([
            [1, nu, 0],
            [nu, 1, 0],
            [0, 0, (1 - nu) / 2]
        ])

        self.K = None
        self.F = None
        self.u = None
        self.stresses = None
        self.von_mises = None

    @staticmethod
    def get_global_dofs(connectivity):
        """Map node indices to displacement DOFs."""
        connectivity = np.asarray(connectivity, dtype=int)

        dofs = np.empty(2 * len(connectivity), dtype=int)
        dofs[0::2] = 2 * connectivity
        dofs[1::2] = 2 * connectivity + 1

        return dofs

    def assemble_stiffness(self):
        """Assemble the global sparse stiffness matrix."""
        n_dofs = 2 * self.mesh.number_of_nodes
        K = sp.lil_matrix((n_dofs, n_dofs))

        Element = (
            T3Element if self.mesh.element_type == "T3"
            else T6Element
        )

        for connectivity in self.mesh.elements:
            nodes = self.mesh.nodes[connectivity]
            element = Element(nodes)

            Ke = element.stiffness_matrix(
                self.D,
                self.config.thickness
            )

            dofs = self.get_global_dofs(connectivity)

            for i, I in enumerate(dofs):
                for j, J in enumerate(dofs):
                    K[I, J] += Ke[i, j]

        self.K = K.tocsr()
        return self.K

    def assemble_forces(self):
        """Assemble equal and opposite tensile boundary tractions."""
        T = self.config.traction
        h = self.config.thickness

        F = np.zeros(2 * self.mesh.number_of_nodes)

        if self.mesh.element_type == "T3":
            weights = np.array([1/2, 1/2])
        else:
            weights = np.array([1/6, 4/6, 1/6])

        for side, sign in (("left", -1), ("right", 1)):
            for edge in self.mesh.get_boundary_edges(side):
                p1 = self.mesh.nodes[edge[0]]
                p2 = self.mesh.nodes[edge[-1]]

                L = np.linalg.norm(p2 - p1)

                F[2 * edge] += sign * T * h * L * weights

        self.F = F
        return self.F

    def apply_boundary_conditions(self):
        """Eliminate constrained displacement DOFs."""
        fixed_dofs = self.mesh.get_fixed_dofs()

        all_dofs = np.arange(2 * self.mesh.number_of_nodes)
        free_dofs = np.setdiff1d(all_dofs, fixed_dofs)

        Kff = self.K[free_dofs, :][:, free_dofs]
        Ff = self.F[free_dofs]

        return Kff, Ff, free_dofs

    def solve(self):
        """Solve the reduced system and reconstruct displacements."""
        self.assemble_stiffness()
        self.assemble_forces()

        Kff, Ff, free_dofs = self.apply_boundary_conditions()

        u = np.zeros(2 * self.mesh.number_of_nodes)
        u[free_dofs] = spsolve(Kff, Ff)

        self.u = u
        return self.u

    def compute_stress(self):
        """Compute and store stresses at element centroids."""
        if self.u is None:
            raise ValueError("Solve the system before computing stresses.")

        stresses = np.zeros((self.mesh.number_of_elements, 3))

        Element = (
            T3Element if self.mesh.element_type == "T3"
            else T6Element
        )

        for i, connectivity in enumerate(self.mesh.elements):
            element = Element(self.mesh.nodes[connectivity])
            dofs = self.get_global_dofs(connectivity)
            ue = self.u[dofs]

            if self.mesh.element_type == "T3":
                B = element.strain_matrix()
            else:
                B = element.strain_matrix(1/3, 1/3)

            stresses[i] = self.D @ B @ ue

        self.stresses = stresses
        return self.stresses

    def compute_von_mises(self):
        """Compute the maximum sampled von Mises stress in each element."""
        if self.u is None:
            raise ValueError("Solve the system before computing stresses.")

        def von_mises(stress):
            sx, sy, txy = stress
            return np.sqrt(sx**2 + sy**2 - sx * sy + 3 * txy**2)

        values = np.zeros(self.mesh.number_of_elements)

        if self.mesh.element_type == "T3":
            if self.stresses is None:
                self.compute_stress()

            values = np.array([von_mises(s) for s in self.stresses])

        else:
            points = [
                (1/3, 1/3),  # Centroid
                (1, 0),      # Node 1
                (0, 1),      # Node 2
                (0, 0),      # Node 3
                (0.5, 0.5),  # Node 4
                (0, 0.5),    # Node 5
                (0.5, 0)     # Node 6
            ]

            for i, connectivity in enumerate(self.mesh.elements):
                element = T6Element(self.mesh.nodes[connectivity])
                dofs = self.get_global_dofs(connectivity)
                ue = self.u[dofs]

                values[i] = max(
                    von_mises(self.D @ element.strain_matrix(L1, L2) @ ue)
                    for L1, L2 in points
                )

        self.von_mises = values
        return self.von_mises
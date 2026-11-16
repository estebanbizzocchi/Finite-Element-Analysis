import numpy as np


class T3Element:
    """Three-node linear triangular finite element."""

    def __init__(self, nodes):
        self.nodes = np.asarray(nodes, dtype=float)

        if self.nodes.shape != (3, 2):
            raise ValueError("T3 element requires three 2D nodes.")

    def area(self, signed=False):
        """Return the signed or absolute triangle area."""
        (xA, yA), (xB, yB), (xC, yC) = self.nodes

        A = 0.5 * (
            xA * (yB - yC)
            + xB * (yC - yA)
            + xC * (yA - yB)
        )

        return A if signed else abs(A)

    def strain_matrix(self):
        """Return the constant strain-displacement matrix (3 x 6)."""
        (xA, yA), (xB, yB), (xC, yC) = self.nodes

        area = self.area(signed=True)

        if abs(area) < 1e-14:
            raise ValueError("Degenerate T3 element.")

        b1 = (yB - yC) / (2 * area)
        b2 = (yC - yA) / (2 * area)
        b3 = (yA - yB) / (2 * area)

        c1 = (xC - xB) / (2 * area)
        c2 = (xA - xC) / (2 * area)
        c3 = (xB - xA) / (2 * area)

        return np.array([
            [b1,  0, b2,  0, b3,  0],
            [ 0, c1,  0, c2,  0, c3],
            [c1, b1, c2, b2, c3, b3]
        ])

    def stiffness_matrix(self, D, thickness):
        """Compute the element stiffness matrix (6 x 6)."""
        A = self.area()
        B = self.strain_matrix()

        return thickness * A * (B.T @ D @ B)


class T6Element:
    """Six-node quadratic isoparametric triangular finite element."""

    # Six-point quadrature; weights are normalized to sum to 1.
    GAUSS_POINTS = (
        (0.445948490915965, 0.445948490915965, 0.223381589678011),
        (0.445948490915965, 0.108103018168070, 0.223381589678011),
        (0.108103018168070, 0.445948490915965, 0.223381589678011),
        (0.091576213509771, 0.091576213509771, 0.109951743655322),
        (0.091576213509771, 0.816847572980458, 0.109951743655322),
        (0.816847572980458, 0.091576213509771, 0.109951743655322)
    )

    def __init__(self, nodes):
        self.nodes = np.asarray(nodes, dtype=float)

        if self.nodes.shape != (6, 2):
            raise ValueError("T6 element requires six 2D nodes.")

    def shape_derivatives(self, L1, L2):
        """Shape-function derivatives with respect to L1 and L2."""
        L3 = 1 - L1 - L2

        return np.array([
            [4*L1-1, 0, 1-4*L3, 4*L2, -4*L2, 4*(L3-L1)],
            [0, 4*L2-1, 1-4*L3, 4*L1, 4*(L3-L2), -4*L1]
        ])

    def jacobian(self, L1, L2):
        """Return the Jacobian of the isoparametric mapping."""
        dN = self.shape_derivatives(L1, L2)
        return dN @ self.nodes

    def strain_matrix(self, L1, L2):
        """Return the strain-displacement matrix (3 x 12)."""
        dN = self.shape_derivatives(L1, L2)
        J = self.jacobian(L1, L2)

        dNxy = np.linalg.solve(J, dN)

        B = np.zeros((3, 12))

        B[0, 0::2] = dNxy[0]
        B[1, 1::2] = dNxy[1]
        B[2, 0::2] = dNxy[1]
        B[2, 1::2] = dNxy[0]

        return B

    def stiffness_matrix(self, D, thickness):
        """Integrate the element stiffness using six-point quadrature."""
        Ke = np.zeros((12, 12))

        for L1, L2, weight in self.GAUSS_POINTS:
            J = self.jacobian(L1, L2)
            detJ = np.linalg.det(J)

            if abs(detJ) < 1e-14:
                raise ValueError(
                    "Invalid T6 element: singular Jacobian."
                )

            B = self.strain_matrix(L1, L2)

            Ke += weight * abs(detJ) * (B.T @ D @ B)

        return 0.5 * thickness * Ke

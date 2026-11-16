import numpy as np


class Mesh:
    def __init__(self, element_type, nodes, elements, boundary_edges):
        """
        Store a T3 or T6 triangular mesh.

        Parameters
        ----------
        element_type : str
            "T3" or "T6".
        nodes : array_like, shape (N, 2)
            Node coordinates.
        elements : array_like, shape (M, 3) or (M, 6)
            Element connectivity (zero-based).
        boundary_edges : dict
            Left and right boundary edges, indexed by side.
            T3: [node1, node2]
            T6: [node1, midpoint, node2]
        """
        self.element_type = element_type
        self.nodes = np.asarray(nodes, dtype=float)
        self.elements = np.asarray(elements, dtype=int)

        n = 2 if element_type == "T3" else 3
        self.boundary_edges = {
            side: np.asarray(edges, dtype=int).reshape(-1, n)
            for side, edges in boundary_edges.items()
        }

        self.validate()

    @property
    def number_of_nodes(self):
        return len(self.nodes)

    @property
    def number_of_elements(self):
        return len(self.elements)

    def validate(self):
        """Check mesh data consistency."""
        if self.element_type not in ("T3", "T6"):
            raise ValueError("Element type must be T3 or T6.")

        n = 3 if self.element_type == "T3" else 6

        if self.nodes.ndim != 2 or self.nodes.shape[1] != 2:
            raise ValueError("Nodes must have shape (N, 2).")

        if not np.all(np.isfinite(self.nodes)):
            raise ValueError("Node coordinates must be finite.")

        if self.elements.ndim != 2 or self.elements.shape[1] != n:
            raise ValueError(f"Elements must have shape (M, {n}).")

        if self.number_of_nodes == 0 or self.number_of_elements == 0:
            raise ValueError("Mesh cannot be empty.")

        if np.any(self.elements < 0) or np.any(
            self.elements >= self.number_of_nodes
        ):
            raise ValueError("Invalid node index in elements.")

        for side in ("left", "right"):
            if side not in self.boundary_edges:
                raise ValueError(f"Missing {side} boundary.")

            edges = self.boundary_edges[side]

            if len(edges) == 0:
                raise ValueError(f"No edges on {side} boundary.")

            if np.any(edges < 0) or np.any(
                edges >= self.number_of_nodes
            ):
                raise ValueError(
                    f"Invalid node index on {side} boundary."
                )

    def get_boundary_edges(self, side):
        """Return boundary edges on the specified side."""
        if side not in self.boundary_edges:
            raise ValueError(f"Unknown boundary: {side}")

        return self.boundary_edges[side]

    def get_fixed_dofs(self):
        """
        Constrain three DOFs to eliminate rigid-body motion:
        ux = uy = 0 at bottom-left,
        uy = 0 at bottom-right.
        """
        x = self.nodes[:, 0]
        y = self.nodes[:, 1]

        bottom_left = np.argmin(x + y)
        bottom_right = np.argmax(x - y)

        return np.array([
            2 * bottom_left,
            2 * bottom_left + 1,
            2 * bottom_right + 1
        ], dtype=int)

    def __str__(self):
        return (
            f"Mesh type: {self.element_type}\n"
            f"Nodes: {self.number_of_nodes}\n"
            f"Elements: {self.number_of_elements}"
        )


def generate_mesh(config):
    """
    Generate a T3 or T6 triangular mesh using Gmsh.
    """
    import gmsh

    config.validate()
    gmsh.initialize()

    try:
        gmsh.model.add("plate")

        # ----------------------------------
        # 1. Geometry
        # ----------------------------------

        rectangle = gmsh.model.occ.addRectangle(
            0, 0, 0,
            config.length,
            config.width
        )

        holes = []

        for hole in config.holes:
            x, y = hole["center"]
            radius = hole["radius"]

            disk = gmsh.model.occ.addDisk(
                x, y, 0, radius, radius
            )

            holes.append((2, disk))

        if holes:
            gmsh.model.occ.cut(
                [(2, rectangle)],
                holes,
                removeObject=True,
                removeTool=True
            )

        gmsh.model.occ.synchronize()

        surfaces = gmsh.model.getEntities(2)

        if not surfaces:
            raise ValueError(
                "The holes have removed the entire plate."
            )

        
        # ----------------------------------
        # 2. Identify loaded boundaries
        # ----------------------------------

        boundary_curves = {"left": [], "right": []}

        tol = 1e-8 * max(config.length, config.width)

        for _, tag in gmsh.model.getEntities(1):
            x, y, z = gmsh.model.occ.getCenterOfMass(1, tag)

            if abs(x) < tol:
                boundary_curves["left"].append(tag)

            elif abs(x - config.length) < tol:
                boundary_curves["right"].append(tag)

        for side, curves in boundary_curves.items():
            if not curves:
                raise ValueError(
                    f"No geometry curves found on {side} boundary."
                )


        # ----------------------------------
        # 3. Mesh generation
        # ----------------------------------

        gmsh.option.setNumber(
            "Mesh.MeshSizeMin",
            config.element_size
        )
        gmsh.option.setNumber(
            "Mesh.MeshSizeMax",
            config.element_size
        )

        gmsh.option.setNumber("Mesh.ElementOrder", 1)

        gmsh.model.mesh.generate(2)

        if config.element_type == "T6":
            gmsh.model.mesh.setOrder(2)

        # ----------------------------------
        # 4. Extract nodes
        # ----------------------------------

        node_tags, coordinates, _ = (
            gmsh.model.mesh.getNodes()
        )

        nodes = coordinates.reshape(-1, 3)[:, :2]

        tag_to_index = {
            int(tag): i
            for i, tag in enumerate(node_tags)
        }

        # ----------------------------------
        # 5. Extract triangular elements
        # ----------------------------------

        element_types, _, element_nodes = (
            gmsh.model.mesh.getElements(dim=2)
        )

        gmsh_type = (
            2 if config.element_type == "T3" else 9
        )

        n = 3 if config.element_type == "T3" else 6
        elements = []

        for element_type, connectivity in zip(
            element_types, element_nodes
        ):
            if element_type != gmsh_type:
                continue

            connectivity = connectivity.reshape(-1, n)

            for element in connectivity:
                elements.append([
                    tag_to_index[int(tag)]
                    for tag in element
                ])

        if not elements:
            raise ValueError(
                "No triangular elements were generated."
            )

        elements = np.asarray(elements, dtype=int)

        # ----------------------------------
        # 6. Extract boundary edges
        # ----------------------------------

        boundary_edges = {
            "left": [],
            "right": []
        }

        line_type = (
            1 if config.element_type == "T3" else 8
        )

        n_edge_nodes = (
            2 if config.element_type == "T3" else 3
        )

        for side, curves in boundary_curves.items():
            for curve in curves:
                types, _, connectivities = (
                    gmsh.model.mesh.getElements(
                        dim=1,
                        tag=curve
                    )
                )

                for element_type, connectivity in zip(
                    types, connectivities
                ):
                    if element_type != line_type:
                        continue

                    connectivity = connectivity.reshape(
                        -1, n_edge_nodes
                    )

                    for edge in connectivity:
                        edge = [
                            tag_to_index[int(tag)]
                            for tag in edge
                        ]

                        # Gmsh T6 line ordering:
                        # [endpoint1, endpoint2, midpoint]
                        if config.element_type == "T6":
                            edge = [edge[0], edge[2], edge[1]]

                        boundary_edges[side].append(edge)

        # ----------------------------------
        # 7. Return mesh
        # ----------------------------------

        return Mesh(
            config.element_type,
            nodes,
            elements,
            boundary_edges
        )

    finally:
        gmsh.finalize()

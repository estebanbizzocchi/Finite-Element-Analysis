from pathlib import Path
import math
import yaml


class Config:
    def __init__(
        self, length, width, thickness, holes,  # Geometry
        young_modulus, poisson_ratio,           # Material
        traction,                               # Loading
        element_type, element_size              # Mesh
    ):
        # Geometry
        self.length = length
        self.width = width
        self.thickness = thickness
        self.holes = holes

        # Material
        self.young_modulus = young_modulus
        self.poisson_ratio = poisson_ratio

        # Loading
        self.traction = traction

        # Mesh
        self.element_type = element_type
        self.element_size = element_size

        self.validate()

    def validate(self):
        """Validate the configuration parameters."""

        # Geometry
        if not all(
            math.isfinite(v) and v > 0
            for v in (self.length, self.width, self.thickness)
        ):
            raise ValueError(
                "Length, width and thickness must be positive and finite."
            )

        for i, hole in enumerate(self.holes):
            center = hole["center"]
            radius = hole["radius"]

            if len(center) != 2:
                raise ValueError(
                    f"Hole {i}: center must contain [x, y]."
                )

            if not all(math.isfinite(x) for x in center):
                raise ValueError(
                    f"Hole {i}: center coordinates must be finite."
                )

            if not math.isfinite(radius) or radius <= 0:
                raise ValueError(
                    f"Hole {i}: radius must be positive and finite."
                )

        # Material
        if not math.isfinite(self.young_modulus) or self.young_modulus <= 0:
            raise ValueError(
                "Young's modulus must be positive and finite."
            )

        if not math.isfinite(self.poisson_ratio) or not (
            -1.0 < self.poisson_ratio < 0.5
        ):
            raise ValueError(
                "Poisson's ratio must satisfy -1 < nu < 0.5."
            )

        # Loading
        if not math.isfinite(self.traction):
            raise ValueError("Traction must be finite.")

        # Mesh
        if self.element_type not in ("T3", "T6"):
            raise ValueError("Element type must be T3 or T6.")

        if not math.isfinite(self.element_size) or self.element_size <= 0:
            raise ValueError(
                "Element size must be positive and finite."
            )

    def __str__(self):
        """Return a readable representation of the configuration."""

        holes_str = "\n".join(
            f"    Hole {i + 1}: center={hole['center']}, "
            f"radius={hole['radius']}"
            for i, hole in enumerate(self.holes)
        )

        if not self.holes:
            holes_str = "    None"

        return (
            "CONFIGURATION\n"
            "------------------------------\n"
            "Geometry:\n"
            f"  Length:       {self.length}\n"
            f"  Width:        {self.width}\n"
            f"  Thickness:    {self.thickness}\n"
            f"  Holes:\n{holes_str}\n\n"
            "Material:\n"
            f"  Young modulus: {self.young_modulus:.3e} Pa\n"
            f"  Poisson ratio: {self.poisson_ratio}\n\n"
            "Loading:\n"
            f"  Traction:      {self.traction:.3e} Pa\n\n"
            "Mesh:\n"
            f"  Element type:  {self.element_type}\n"
            f"  Element size:  {self.element_size}\n"
            "------------------------------"
        )

def read_parameters(experiment_path):
    """
    Read config.yaml from an experiment directory
    and return a validated Config object.
    """

    config_path = Path(experiment_path) / "config.yaml"

    with open(config_path, "r", encoding="utf-8") as file:
        data = yaml.safe_load(file)

    geometry = data["geometry"]
    material = data["material"]
    mesh = data["mesh"]

    return Config(
        length=float(geometry["length"]),
        width=float(geometry["width"]),
        thickness=float(geometry["thickness"]),
        holes=geometry.get("holes", []),
        young_modulus=float(material["young_modulus"]),
        poisson_ratio=float(material["poisson_ratio"]),
        traction=float(data["traction"]),
        element_type=mesh["element_type"],
        element_size=float(mesh["element_size"])
    )

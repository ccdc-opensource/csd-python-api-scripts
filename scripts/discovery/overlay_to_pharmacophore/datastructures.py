"""
This script can be used for any purpose without limitation subject to the
conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
This permission notice and the following statement of attribution must be
included in all copies or substantial portions of this script.

"07/08/2026": created by the Cambridge Crystallographic Data Centre
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np


@dataclass
class FeatureTolerances:
    """
    Allowed feature types and their tolerances.

    Each field is a valid feature label; its value is the tolerance:
    a single weight, or a (parent, virtual_point) pair for projected features.
    """
    acceptor: float = 1.0
    acceptor_projected: tuple[float, float] = (0.8, 0.8)
    donor_projected: tuple[float, float] = (1.0, 1.0)
    hydrophobe: float = 1.0
    ring_planar_projected: tuple[float, float] = (1.0, 1.0)
    ring_non_planar: tuple[float, float] = (1.0, 1.0)
    halogen: float = 1.0

    def __getitem__(self, key: str) -> float | tuple[float, float]:
        return getattr(self, key)


class PharmFeaturePoint(np.ndarray):
    def __new__(
        cls,
        *coordinates: float | Iterable[float],
        label: Optional[str] = None,
        virtual_point: Optional[np.ndarray] = None,
    ):
        if len(coordinates) == 1:
            arr = np.asarray(coordinates[0], dtype=float)
        elif len(coordinates) == 3:
            arr = np.asarray(coordinates, dtype=float)
        else:
            raise TypeError("Coordinates must be either an iterable of length three, or three floats")
        if arr.shape != (3,):
            raise ValueError("Coordinates must be a 3-element array")
        obj = arr.view(cls)
        obj.label = label
        obj.virtual_point = virtual_point
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.label = getattr(obj, 'label', None)
        self.virtual_point = getattr(obj, 'virtual_point', None)

    def __repr__(self) -> str:
        return (
            f"PharmFeaturePoint({self.x}, {self.y}, {self.z}, "
            f"label={self.label}, virtual_point={self.virtual_point}), "
        )

    def __str__(self) -> str:
        return self.__repr__()

    @property
    def coordinates(self) -> np.ndarray:
        return np.asarray(self)

    @property
    def x(self):
        return self[0]

    @property
    def y(self):
        return self[1]

    @property
    def z(self):
        return self[2]

    @property
    def tolerance(self) -> float | tuple[float, float]:
        """
        Get the tolerance for the feature based on its label.
        Returns:
            A single tolerance for features with a single tolerance or a tuple for those with two tolerances.
        """
        if self.label is None:
            raise LookupError("Feature label must be set to determine tolerances.")
        return FeatureTolerances()[self.label]

    @property
    def weight_parent(self) -> float:
        if isinstance(self.tolerance, float):
            return self.tolerance
        elif isinstance(self.tolerance, tuple):
            return self.tolerance[0]
        else:
            raise ValueError("Incorrect tolerances loaded for feature point.")

    @property
    def weight_vp(self) -> float:
        return self.tolerance[1]



@dataclass
class OverlayData:
    input_folder: Path
    output_folder: Path
    # Pharmacophore file from the pharmacophores folder
    pharm_file: Path
    # Overlay solution file (the chosen one from the many produced)
    overlay_file: Path

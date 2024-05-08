#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2023-06-21: created by Pietro Sacchi, The Cambridge Crystallographic Data Centre
# 2024-04-15: modified by Alex Moldovan, The Cambridge Crystallographic Data Centre
#
from typing import TypeVar, TYPE_CHECKING

MorphologyBase = TypeVar('MorphologyBase')

if TYPE_CHECKING:
    from ccdc.morphology import MorphologyBase


class ShapeClassification:
    """
    Calculate parameters for shape classification as defined by Angelidakis et al. in
    Powder Technology, 396 (2022), 689-695
    """

    def __init__(self, major: float, medium: float, minor: float):
        self.minor_length = minor
        self.medium_length = medium
        self.major_length = major

    @classmethod
    def from_morphology(cls, morphology: MorphologyBase):
        """Return object from morphology"""
        return cls(major=morphology.oriented_bounding_box.major_length,
                   medium=morphology.oriented_bounding_box.median_length,
                   minor=morphology.oriented_bounding_box.minor_length)

    def shape_classification_data(self):
        """ returns data as dictionary for database """
        return {"shape_classification": self.shape_description,
                'S/M': self.minor_length / self.medium_length,
                'M/L': self.medium_length / self.major_length}

    @property
    def elongation(self) -> float:
        return (self.major_length * self.minor_length) / (
                self.major_length * self.minor_length + self.medium_length ** 2) - self.minor_length / (
                self.major_length + self.minor_length)

    @property
    def flatness(self) -> float:
        return self.medium_length ** 2 / (
                self.major_length * self.minor_length + self.medium_length ** 2) - self.minor_length / (
                self.major_length + self.minor_length)

    @property
    def compactness(self) -> float:
        return 2 * self.minor_length / (self.major_length + self.minor_length)

    @property
    def shape_description(self) -> str:
        """Return the classification of the shape"""
        if self.elongation >= 0.2:
            if self.flatness >= 0.2:
                return "Lath"
            else:
                return "Needle"
        if self.elongation <= 0.2:
            if self.flatness <= 0.2:
                return "Block"
            else:
                return "Plate"

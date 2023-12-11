#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#

"""
Script to classify the shape of a crystal morphology as belonging to one of four classes:

    - block
    - plate
    - lath
    - needle

The shape classification follows the definition of Angelidakis, Nadimi and Utili found in
Powder Technology, 396 (2022), 689-695 DOI: 10.1016/j.powtec.2021.11.027

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from cycler import cycler


def _elongation_line():
    """
    This function's purpose is exclusively that of calculating the lines that separate different regions of the plot.
    Shape classification is not based on this function.

    Finds the line corresponding to aspect ratios with elongation = 0.2
    """

    # set b = 1, then solve elongation = 0.2 for a using a solver
    # ac/(ac+1)-c/(a+c) = 0.2

    values = np.linspace(0.0001, 1, 100)

    def func(c):
        # these are the solutions of the equation to get a as a function of c
        sol1 = (c ** 2 + 1 - np.sqrt(c ** 4 + 98 * c ** 2 + 1)) / (8 * c)
        sol2 = (c ** 2 + 1 + np.sqrt(c ** 4 + 98 * c ** 2 + 1)) / (8 * c)
        return sol1, sol2

    # get the "elongation" line
    x_values = []
    y_values = []
    for c in values:
        x_values.append(c)
        # we keep positive solutions only
        _, a = func(c)
        y_values.append(1 / a)
    # x_values are the c/b ratios
    # y_values are the b/a ratios
    return x_values, y_values


def _flatness_line():
    """
    This function's purpose is exclusively that of calculating the lines that separate different regions of the plot.
    Shape classification is not based on this function.

    Finds the line corresponding to aspect ratios with flatness = 0.2
    """
    # set a = 1, then solve flatness = 0.2 for c using a solver
    # b^2/(c+b^2) - c/(c+1) = 0.2
    values = np.linspace(0.0001, 1, 100)

    def func(b):
        # these are the solutions of the equation to get c as a function of b
        sol1 = (- b ** 2 - 1 + np.sqrt(b ** 4 + 98 * b ** 2 + 1)) / 12
        sol2 = (- b ** 2 + 1 + np.sqrt(b ** 4 + 98 * b ** 2 + 1)) / 12
        return sol1, sol2

    # get the "flatness" line
    x_values = []
    y_values = []
    for b in values:
        y_values.append(b)
        # we keep positive solutions only
        c, _ = func(b)
        x_values.append(c / b)
    # x_values are the c/b ratios
    # y_values are the b/a ratios
    return x_values, y_values


def plot(shape):
    """
    Prepare a Zingg plot for a single morphology.

    :param shape: an instance of the ShapeClassifier class
    """
    fig, ax = plt.subplots()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xlabel("S / M")
    ax.set_ylabel("M / L")
    ax.text(0.01, 0.95, "PLATE")
    ax.text(0.85, 0.01, "NEEDLE")
    ax.text(0.85, 0.95, "BLOCK")
    ax.text(0.2, 0.2, "LATH")

    # calculate the dividing lines for the plot
    # line where elongation = 0.2
    el_line = elongation_line()
    # line where flatness = 0.2
    fl_line = flatness_line()
    ax.plot(el_line[0], el_line[1], lw=1, c="k")
    ax.plot(fl_line[0], fl_line[1], lw=1, c="k")
    # plot lines of classic Zingg plot
    ax.axvline(x=2 / 3, c="grey", lw=0.5, ls="--")
    ax.axhline(y=2 / 3, c="grey", lw=0.5, ls="--")
    # finally plot the data
    # can add additional arguments here to make prettier
    ax.scatter(self.minor_length / self.median_length, self.median_length / self.major_length)
    plt.show()


def plot_multiple_shapes(shapes):
    """
    Make a Zingg plo for multiple morphologies.

    :param shapes: a list of ShapeClassifier objects
    """

    plt.rcParams.update({'font.family': 'Arial'})
    fig, ax = plt.subplots()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xlabel("S / M")
    ax.set_ylabel("M / L")
    ax.text(0.01, 0.95, "PLATE")
    ax.text(0.85, 0.01, "NEEDLE")
    ax.text(0.85, 0.95, "BLOCK")
    ax.text(0.2, 0.2, "LATH")

    # calculate the dividing lines
    # line where elongation = 0.2
    el_line = _elongation_line()
    # line where flatness = 0.2
    fl_line = _flatness_line()
    ax.plot(el_line[0], el_line[1], lw=1, c="k")
    ax.plot(fl_line[0], fl_line[1], lw=1, c="k")
    # plot lines of classic Zingg plot
    ax.axvline(x=2 / 3, c="grey", lw=0.5, ls="--")
    ax.axhline(y=2 / 3, c="grey", lw=0.5, ls="--")
    # define a colormap
    colours = cm.bwr(np.linspace(0, 1, len([item for item in shapes if item != 0])))
    # finally plot the data
    ax.set_prop_cycle(cycler('color', colours))
    for shape in shapes:
        ax.scatter(shape.minor_length / shape.median_length, shape.median_length / shape.major_length)
    plt.show()


class ShapeClassifier:

    def __init__(self, morphology: object):
        """
        param morphology: an instance of ccdc.morphology classes
        """
        self.morphology = morphology
        self.bounding_box = self.morphology.oriented_bounding_box
        self.major_length, self.median_length, self.minor_length = self.bounding_box.lengths

    @property
    def elongation(self) -> float:
        return (self.major_length * self.minor_length) / (
            self.major_length * self.minor_length + self.median_length ** 2) - self.minor_length / (
                   self.major_length + self.minor_length)

    @property
    def flatness(self) -> float:
        return self.median_length ** 2 / (
            self.major_length * self.minor_length + self.median_length ** 2) - self.minor_length / (
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
        if self.elongation < 0.2:
            if self.flatness < 0.2:
                return "Block"
            else:
                return "Plate"


if __name__ == "__main__":
    # an example of how to use this script

    from ccdc.io import EntryReader
    from ccdc.morphology import BFDHMorphology

    # let's calculate a morphology!
    refcode = "VUKRAW"
    reader = EntryReader("CSD")
    entry = reader.entry(refcode)
    crystal = entry.crystal

    bfdh_morphology = BFDHMorphology(crystal)

    shape_1 = ShapeClassifier(bfdh_morphology)
    print(f"{refcode} BFDH morphology is classified as a: {shape_1.shape_description}")

    # let's try with another morphology!
    refcode = "AABHTZ"
    reader = EntryReader("CSD")
    entry = reader.entry(refcode)
    crystal = entry.crystal

    bfdh_morphology = BFDHMorphology(crystal)

    shape_2 = ShapeClassifier(bfdh_morphology)
    print(f"{refcode} BFDH morphology is classified as a: {shape_2.shape_description}")

    # now let's make a Zingg plot to visualise our particles' shapes
    plot_multiple_shapes([shape_1, shape_2])

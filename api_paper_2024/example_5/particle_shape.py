#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2021-11-29: created by Alex Moldovan, The Cambridge Crystallographic Data Centre
# 2024-04-15: modified by Alex Moldovan, The Cambridge Crystallographic Data Centre
#
import argparse
import warnings
from typing import TypeVar

from ccdc.io import CrystalReader
from ccdc.morphology import BFDHMorphology

MorphologyBase = TypeVar('MorphologyBase')
Crystal = TypeVar('Crystal')

try:
    from plotly import graph_objects as go
    from plotly.subplots import make_subplots
    from plotly import offline

except ModuleNotFoundError:
    warnings.warn("Plotly could not be found, please install plotly using `conda install plotly`")

from visualiser import PlotlyParticle, PlotlyZingg


def calculate_morphology(structure: Crystal) -> MorphologyBase:
    """Calculate the morphology of a crystal as an example we use the BFDH"""
    morphology = BFDHMorphology(structure)
    return morphology


def load_crystal(name: str) -> Crystal:
    """Load crystal. In this case we will use a refcode. But you can substitute any cif/mol2 method in here."""
    return CrystalReader('CSD').crystal(name)


def combine_graphs(zingg: PlotlyZingg, particle: PlotlyParticle) -> go.Figure:
    """ Puts the two graphs together"""
    combined_fig = make_subplots(rows=1, cols=2,
                                 specs=[[{'type': 'scatter'}, {'type': 'scatter3d'}]],
                                 column_widths=[0.6, 0.4])
    combined_fig.add_traces(zingg.fig.data, rows=1, cols=1)
    for shape in zingg.fig.layout.shapes:
        combined_fig.add_shape(shape, row=1, col=1)
    for annotation in zingg.fig.layout.annotations:
        combined_fig.add_annotation(annotation, row=1, col=1)
    combined_fig.add_traces(particle.fig.data, rows=1, cols=2)
    combined_fig.update_scenes(particle.fig.layout.scene)
    combined_fig.update_layout(yaxis=dict(title='M / L', range=[0, 1]),
                               xaxis=dict(title='S / M', range=[0, 1]),
                               template='simple_white',
                               font_family="Courier New",
                               font_size=20,
                               title="Shape Classification",
                               scene_aspectmode="data",
                               legend=dict(
                                   orientation="h",
                                   yanchor="bottom",
                                   xanchor="center",
                                   x=0.5,
                                   y=1
                               ))
    return combined_fig


def main(args: argparse.Namespace) -> None:
    """ Runs the classification and pieces the plot together """
    ref = args.refcode

    crystal = load_crystal(name=ref)
    morphology = calculate_morphology(structure=crystal)

    plotly_particle = PlotlyParticle(morphology=morphology)
    zingg_plot = PlotlyZingg(zone_opacity=0.1)
    zingg_plot.plot_shape(major_length=morphology.oriented_bounding_box.major_length,
                          medium_length=morphology.oriented_bounding_box.median_length,
                          minor_length=morphology.oriented_bounding_box.minor_length, name=ref)

    figure = combine_graphs(zingg=zingg_plot, particle=plotly_particle)
    offline.plot(figure_or_data=figure, filename=f"{ref}_shape_classification.html")


def arg_parser() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Particle Shape Classifier")
    parser.add_argument("refcode", help="Refcode referencing the CSD or in-house database")
    return parser.parse_args()


if __name__ == "__main__":
    args = arg_parser()
    main(args)

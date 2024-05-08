#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2021-11-29: created by Alex Moldovan, The Cambridge Crystallographic Data Centre
# 2024-04-01: modified by Alex Moldovan & Pietro Sacchi, The Cambridge Crystallographic Data Centre
#
import warnings
from typing import List, Tuple, TypeVar

import numpy as np

MorphologyBase = TypeVar('MorphologyBase')

try:
    from plotly import graph_objects as go
    from plotly.subplots import make_subplots

except ModuleNotFoundError:
    warnings.warn("Plotly could not be found, please install plotly using `conda install plotly`")


class PlotlyParticle:
    def __init__(self, morphology: MorphologyBase, colour: str = "rgba(175, 193, 242, 0.8)"):
        self.fig = self.generate_morphology_plot(morphology=morphology, colour=colour)
        self.fig.update_scenes(camera_projection_type="orthographic")
        self.fig.update_layout(scene_aspectmode='data')

    @staticmethod
    def generate_morphology_plot(morphology: MorphologyBase, colour='teal'):  # noqa:E194 #NOSONAR
        """Plots a basic morphology for visualisation"""
        fig = go.Figure()

        for facet in morphology.facets:
            verts = np.array([[*f] for f in facet.coordinates])
            verts = np.concatenate([verts, [verts[0]]])
            surface_axis = 2
            hkl = facet.miller_indices.hkl
            vector = np.array([*facet.plane.normal])
            if (hkl[0] != 0) and (hkl[1] != 0) and (hkl[2] == 0):
                surface_axis = 1
                if np.allclose(np.cross(vector, np.array([1, 0, 0])), 0):
                    surface_axis = 0
            fig.add_scatter3d(x=verts[:, 0], y=verts[:, 1], z=verts[:, 2],
                              mode='lines',
                              surfaceaxis=surface_axis,
                              surfacecolor=colour,
                              opacity=0.5,
                              name=str([*facet.miller_indices.hkl]),
                              text=str([*facet.miller_indices.hkl]),
                              line=dict(width=5, color='black'), hoverinfo="name", showlegend=False)

            fig.add_scatter3d(x=[facet.centre_of_geometry[0]],
                              y=[facet.centre_of_geometry[1]],
                              z=[facet.centre_of_geometry[2]],
                              text=str(facet.miller_indices.hkl), textfont=dict(size=16),
                              mode="text", showlegend=False)

        fig.update_scenes(bgcolor='white',
                          yaxis=dict(visible=False),
                          xaxis=dict(visible=False),
                          zaxis=dict(visible=False))

        return fig


class PlotlyZingg:
    def __init__(self, zone_opacity: float = 0.1):
        self.fig = None
        self.zone_opacity = zone_opacity
        self.construct_base_plot()

    @staticmethod
    def elongation_line() -> Tuple[List[float], List[float]]:
        """
        Following the definition of elongation give in eq. 5 of
        https://www.sciencedirect.com/science/article/pii/S0032591021009785.

        This function calculates the x and y values needed to draw a line which corresponds to an elongation value of
        0.2 (compare with Figure 5 of said paper). In this case x = c/b and y = b/a

        a, b, c are the major, median and minor sides of the bounding box.

        To solve, we set b=1 and solve for a.
        """

        # set b = 1, then solve for a using a solver
        # ac/(ac+1)-c/(a+c) = 0.2

        def solve_func(c):
            # these are the solutions of the equation to get a as a function of c
            # we have the negative roots
            sol1 = (c ** 2 + 1 - np.sqrt(c ** 4 + 98 * c ** 2 + 1)) / (8 * c)
            # and the positive roots
            sol2 = (c ** 2 + 1 + np.sqrt(c ** 4 + 98 * c ** 2 + 1)) / (8 * c)
            return sol1, sol2

        # get the "elongation" line
        # since we have set b=1, x = c/b = c
        x_values = np.linspace(0.0001, 1, 100)
        # Since a is one of the box dimensions, negative values don't have meaning, and we can ignore them
        # we keep positive solutions only
        y_values = 1 / solve_func(x_values)[1]

        # x_values are the c/b ratios
        # y_values are the b/a ratios
        return x_values, y_values

    @staticmethod
    def flatness_line() -> Tuple[List[float], List[float]]:
        """
        Following the definition of flatness give in eq. 5 of
        https://www.sciencedirect.com/science/article/pii/S0032591021009785.

        This function calculates the x and y values needed to draw a line which corresponds to a flatness value of
        0.2 (compare with Figure 5 of said paper). In this case x = c/b and y = b/a

        a, b, c are the major, median and minor sides of the bounding box.

        To solve, we set a=1 and solve for c.
        """

        # set a = 1, then solve flatness = 0.2 using a solver
        # b^2/(c+b^2) - c/(c+1) = 0.2

        def solve_func_b(b):
            # these are the solutions of the equation to get c as a function of b
            # we get two roots
            # Since c is one of the box dimensions, negative values don't have meaning, and we can ignore them
            # we have the positive roots
            sol1 = (- b ** 2 - 1 + np.sqrt(b ** 4 + 98 * b ** 2 + 1)) / 12
            # and we have the negative roots
            sol2 = - (b ** 2 + 1 + np.sqrt(b ** 4 + 98 * b ** 2 + 1)) / 12
            return sol1, sol2

        # get the "flatness" line
        # since we have set a = 1, y = b/a = b
        y_values = np.linspace(0.0001, 1, 1000)
        # Since c is one of the box dimensions, negative values don't have meaning, and we can ignore them
        # we keep positive solutions only
        x_values = solve_func_b(y_values)[0] / y_values

        # x_values are the c/b ratios
        # y_values are the b/a ratios
        return x_values, y_values

    def construct_base_plot(self):
        el_line = self.elongation_line()
        # line where flatness = 0.2
        fl_line = self.flatness_line()
        self.fig = go.Figure()

        self.fig.add_annotation(x=0.3, y=0.8, text="Plate", showarrow=False)
        self.fig.add_annotation(x=0.8, y=0.3, text="Needle", showarrow=False)
        self.fig.add_annotation(x=0.8, y=0.8, text="Block", showarrow=False)
        self.fig.add_annotation(x=0.3, y=0.3, text="Lath", showarrow=False)

        el_line_trace = go.Scatter(x=el_line[0],
                                   y=el_line[1],
                                   mode="lines",
                                   fill='tonexty', hoverinfo='skip',
                                   fillcolor=f'rgba(240,225,168,{self.zone_opacity})',
                                   showlegend=False)

        blank_line = go.Scatter(x=[0, 1],
                                y=[1, 1],
                                mode='lines',
                                fill='tonexty',
                                fillcolor=f'rgba(240,145,234,{self.zone_opacity})',
                                showlegend=False, hoverinfo='skip',
                                line=dict(color='rgba(0,0,0,0)'))

        fl_line_trace = go.Scatter(x=fl_line[0],
                                   y=fl_line[1],
                                   mode="lines", hoverinfo='skip',
                                   showlegend=False)

        blank_line2 = go.Scatter(x=[1, 1],
                                 y=[0, 1],
                                 mode='lines',
                                 fill='tonextx', hoverinfo='skip',
                                 fillcolor=f'rgba(120,240,226,{self.zone_opacity})',
                                 line=dict(color='rgba(0,0,0,0)'),
                                 showlegend=False)

        self.fig.add_vline(x=0.666, line_width=1, line_dash="dash", line_color="rgba(0,0,0,0.5)")
        self.fig.add_hline(y=0.666, line_width=1, line_dash="dash", line_color="rgba(0,0,0,0.5)")
        self.fig.add_traces([el_line_trace, blank_line, fl_line_trace, blank_line2])

    def plot_shape(self, minor_length: float, medium_length: float, major_length: float,
                   name: str = 'Structure') -> None:
        marker = go.Scatter(x=[minor_length / medium_length], y=[medium_length / major_length],
                            mode='markers', name=name, marker=dict(color='rgba(125,0,0,0.5)', size=16,
                                                                   line=dict(color='black', width=1.5)),
                            hovertemplate=f'<b>S/M</b>: {round(minor_length / medium_length, 2)} <br>'
                                          f'<b>M/L</b>: {round(medium_length / major_length, 2)} <br>'
                                          f'<b>Volume</b>:{round(medium_length * major_length * minor_length, 2)} A<sup>3<extra></extra>')
        self.fig.add_trace(marker)

        self.fig.update_layout(yaxis=dict(title='M / L', range=[0, 1]),
                               xaxis=dict(title='S / M', range=[0, 1]),
                               template='simple_white',
                               font_family="Courier New",
                               font_size=20,
                               scene_aspectmode="data")

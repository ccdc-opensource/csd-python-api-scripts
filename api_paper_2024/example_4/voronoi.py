#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2023-07-05 Created by Chris Kingsbury, the Cambridge Crystallographic Data Centre
# ORCID 0000-0002-4694-5566
#
# ccdc interface to Voronoi polyhedra which encapsulate the atomic / metal domains
#
#

import argparse
import json
import sys
from os import getcwd
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.linalg import norm
from scipy.spatial import ConvexHull, Voronoi

try:
    import plotly.graph_objects as go
    import pyvoro
except ImportError:
    print("pyvoro and plotly import error - use 'pip install plotly pyvoro' in miniconda (powershell)")

import ccdc.io
from ccdc.utilities import ApplicationInterface

WEIGHTING = json.load(open(Path(__file__).parent / "atom_weighting.json", "r"))


def weighting_from_name(weight_name):
    return WEIGHTING.get(weight_name, WEIGHTING["equal"])


DEFAULT_SETTINGS = {
    "radius": 20.0,
    # Angstroms; radius for the substructure contacts
    "plot_fraction": 0.95,
    # The amount by which each individual cell is contracted
    "output": "html",
    # Whether the output will be the plotly html or a pkl for future compatibility
    "background": False,
    # whether to draw the grid backing grid
    "opacity": 0.8,
    # opacity of the voronoi cells in the plot
    "metal_only": True,
    # Whether to run the molecular or metal-only representation
    "weighting": None,
    # how to bias the generation of the voronoi representation - by "vdw", "equal", "empirical", "calculated", "Bondi", "Batsanov", "Alvarez"
    # https://doi.org/10.1063%2F1.1725697
    # https://doi.org/10.1063/1.1712084
    # https://doi.org/10.1021%2Fjp8111556
    # https://pubs.rsc.org/en/content/articlehtml/2013/dt/c3dt50599e
}

PLOTLY_HEADER = """ <head><script src="https://cdn.plot.ly/plotly-2.20.0.min.js" charset="utf-8"></script>
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <style>.column {  float: left;   width: 45%;   padding: 10px;}
    .row:after { content: ""; display: table; clear: both;}  </style></head>"""

ATOMS_COLOR_DICT = {
    "C": "black",
    "H": "floralwhite",
    "N": "cornflowerblue",
    "S": "yellow",
    "O": "red",
    "Cs": "grey",
    "Li": "violet",
    "Cl": "green",
    "F": "lime",
    "Fe": "firebrick",
    "Mn": "magenta",
    "P": "orange",
    "Br": "maroon",
    "I": "darkviolet",
    "Na": "purple",
    "K": "violet",
    "B": "pink",
    "Cu": "lightblue",
    "Eu": "lightgreen",
}


def generate_superlist(pairs, chain):
    for i, pair in enumerate(pairs):
        for j, at in enumerate(chain):
            if at == pair[0]:
                chain.insert(j, pair[0])
                chain.insert(j + 1, pair[1])
                pairs.pop(i)
                return (pairs, chain)
    return (pairs, chain)


# the formula for determining the area of a triangle of three points in 3 dimensions
def area_tri(a, b, c):
    return np.linalg.norm(np.cross((b - a), (c - a))) / 2


# the area of a polygon using the above formula
def area_poly(arr):
    return sum([area_tri(arr[0], arr[i + 1], arr[i + 2]) for i in range(len(arr) - 2)])


class CrystalVoronoi:
    # initialises the class - this takes the general crystal object and uses atom positions (all or a subset) to determine a
    # Voronoi object. This has some uses in analysis of materials, in understanding how molecules interact, as well in the visualisation
    # of interacting molecular moieties.
    def __init__(
            self,
            ccdc_crystal,
            metal_only=False,
            weighting=None,
            radius=15,
    ):
        self.metal_only = metal_only
        self.crystal = ccdc_crystal
        if self.crystal.molecule.all_atoms_have_sites is False:
            raise ValueError("Some atoms do not have sites, please edit structure and re-calculate.")
        self.radius = radius
        self._shell = None

        if not weighting:
            self.weighting = weighting_from_name("equal")
        elif isinstance(weighting, str):
            self.weighting = weighting_from_name(weighting)
        else:
            self.weighting = weighting

    @property
    def metal_atoms(self):
        if self.shell:
            return [x for x in self.shell if x.is_metal]
        return None

    @property
    def shell(self):
        if self._shell:
            return self._shell
        self._shell = [
            x for x in
            self.crystal.molecular_shell(distance_range=(0.0, float(self.radius))).atoms + self.molec.atoms]
        return self._shell

    @shell.setter
    def shell(self, value):
        self._shell = value

    # This calculates the voronoi division of the space using pyvoro
    def calculate_voronoi(self):
        # The entire molecule for the atom plotting
        self.molec = self.crystal.molecule
        # all of the atoms that will be involved as the external faces of the
        # Voronoi polyhedra; ensures that these are of finite size
        if self.metal_only:
            self.shell = self.metal_atoms
        # prevents duplication of atoms (with multiple overlapping shells)
        self.calculate_unique_atoms()
        all_labels = [x.label for x in self.shell]
        self.uq_labels = [all_labels[y] for y in self.uq_indices]
        self.uq_atoms = [self.shell[y] for y in self.uq_indices]

        # Voronoi object - want to replace with something that will do
        # weighted voronoi by atom pair

        all_weights = [self.weighting.get(x.atomic_symbol, 100) for x in self.shell]
        self.uq_weights = [all_weights[y] for y in self.uq_indices]

        lims = [[min(x) - 1, max(x) + 1] for x in self.uq_coords.T]
        pyv = pyvoro.compute_voronoi(self.uq_coords, lims, 12, radii=self.uq_weights)
        self.pyv = pyv

        self.contacts = []
        for comp_index, comp in enumerate(self.molec.components):
            plot_coords = np.atleast_2d(
                np.array([np.array(x.coordinates) for x in comp.atoms if (x.coordinates is not None)]))

            vor_coords = plot_coords
            if self.metal_only:
                vor_coords = np.atleast_2d(plot_coords[[x.is_metal for x in comp.atoms]])

            if vor_coords.shape == (0, 3):
                print(f"component {str(comp_index)} contains no polyhedra")
                continue
            comp_indices = np.array(
                [i for i, x in enumerate(self.uq_coords)
                 if min(np.linalg.norm(x - vor_coords, axis=1)) < 0.1])

            for at_ix in comp_indices:
                at1 = self.uq_atoms[at_ix]
                verts = np.vstack(pyv[at_ix]["vertices"])
                for face in pyv[at_ix]["faces"]:
                    at2 = self.uq_atoms[face["adjacent_cell"]]
                    distance = np.linalg.norm(np.array(at1.coordinates) - np.array(at2.coordinates))
                    face_size = area_poly(verts[face["vertices"]])
                    if face["adjacent_cell"] > 0:
                        self.contacts.append([at1.atomic_symbol, at2.atomic_symbol, distance, face_size])

    def calculate_unique_atoms(self):
        self.uq_coords, self.uq_indices = np.unique(
            np.atleast_2d(np.array([np.array(x.coordinates) for x in self.shell])),
            axis=0,
            return_index=True,
        )

    # as above, using scipy.spatial - this has the limitation of only allowing equal weighting
    def calculate_voronoi_scipy(self):
        # The entire molecule for the atom plotting
        self.molec = self.crystal.molecule
        # all of the atoms that will be involved as the external faces of the
        # Voronoi polyhedra; ensures that these are of finite size

        if self.metal_only:
            self.shell = self.metal_atoms
        # prevents duplication of atoms (with multiple overlapping shells)
        self.uq_coords, self.uq_indices = np.unique(
            np.atleast_2d(np.array([np.array(x.coordinates) for x in self.shell])),
            axis=0,
            return_index=True,
        )
        self.tv = Voronoi(self.uq_coords)
        self.calculate_unique_atoms()
        all_labels = [x.label for x in self.shell]
        self.uq_labels = [all_labels[y] for y in self.uq_indices]
        self.uq_atoms = [self.shell[y] for y in self.uq_indices]

    # generates a plot from calculated voronoi divisions
    def generate_plot(
            self, plot_fraction=0.95, opacity=1, background=True, voronoi_type="pyvoro"
    ):
        if self.metal_only:
            cmin, cmax = 2, 12
        else:
            cmin, cmax = 0.5, 3

        self.figs, self.polys = [], []
        self.testing = []

        for comp_index, comp in enumerate(self.molec.components):
            plot_coords = np.atleast_2d(
                np.array([np.array(x.coordinates)
                          for x in comp.atoms
                          if (x.coordinates is not None)]))
            if self.metal_only:
                vor_coords = np.atleast_2d(plot_coords[[x.is_metal for x in comp.atoms]])
            else:
                vor_coords = plot_coords

            self.testing.append("met " + str(vor_coords))

            if vor_coords.shape == (0, 3):
                print(f"component {str(comp_index)} contains no polyhedra")
                continue
            comp_indices = np.array(
                [i
                 for i, x in enumerate(self.uq_coords)
                 if min(norm(x - vor_coords, axis=1)) < 0.1])

            self.testing.append("ind " + str(comp_indices))
            all_poly = []
            points, simplices, equations, atixes = [], [], [], []
            off = 0

            fig = go.Figure()
            labels = [y.label for y in comp.atoms]
            pairs = [[labels.index(x.label) for x in z.atoms] for z in comp.bonds]
            superlist = [0]

            for _ in range(len(pairs)):
                pairs, superlist = generate_superlist(pairs, superlist)
                pairs, superlist = generate_superlist(
                    [x[::-1] for x in pairs], superlist
                )

            xs, ys, zs = plot_coords[superlist].T
            asys = np.array([x.atomic_symbol for x in comp.atoms])[superlist]
            ats_labels = [labels[x] for x in superlist]
            ats_colours = [ATOMS_COLOR_DICT.get(x, "gold") for x in asys]

            fig.add_trace(
                go.Scatter3d(
                    x=xs,
                    y=ys,
                    z=zs,
                    hovertext=ats_labels,
                    mode="lines+markers",
                    line={"color": "black"},
                    marker={"color": ats_colours},
                    name="atoms",
                )
            )

            for at_ix in comp_indices:
                if voronoi_type == "scipy":
                    poly = self.tv.vertices[
                        self.tv.regions[self.tv.point_region[at_ix]]]
                else:
                    poly = np.vstack(self.pyv[at_ix].get("vertices"))

                conv = ConvexHull(poly)
                points.append(conv.points)
                simplices.append(conv.simplices + off)
                equations.append(conv.equations)
                [atixes.append(at_ix) for _ in range(len(conv.equations))]
                off += len(conv.points)
                all_poly.append(poly)

                x2, y2, z2 = (
                        plot_fraction * conv.points
                        + (1 - plot_fraction) * self.uq_coords[at_ix]
                ).T
                i2, j2, k2 = conv.simplices.T
                colours = abs(
                    np.dot(conv.equations[:, :3], self.uq_coords[at_ix])
                    + conv.equations[:, 3])

                if voronoi_type == "scipy":
                    colours = (colours * 2)
                else:
                    colours = (colours - self.uq_weights[at_ix])

                hover_info = self.uq_labels[at_ix]
                trace_name = f"{self.uq_labels[at_ix]}"
                if self.metal_only:
                    cmin, cmax = (0.0, 12.0)
                else:
                    cmin, cmax = (-0.5, 1.0)

                fig.add_trace(
                    go.Mesh3d(
                        x=x2,
                        y=y2,
                        z=z2,
                        i=i2,
                        j=j2,
                        k=k2,
                        intensity=colours,
                        intensitymode="cell",
                        colorscale="Spectral",
                        cmin=cmin,
                        cmax=cmax,
                        showlegend=True,
                        showscale=True,
                        opacity=opacity,
                        flatshading=True,
                        hovertext=hover_info,
                        text=np.round(colours, 2),
                        hoverinfo="text+name",
                        name=trace_name,
                    )
                )

            fig.update_layout(
                title=self.crystal.identifier,
                height=1000,
                legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
            )
            fig.layout.scene.camera.projection.type = "orthographic"
            if background:
                fig.update_scenes(
                    xaxis_visible=False, yaxis_visible=False, zaxis_visible=False
                )

            self.figs.append(fig)
            self.polys.append(all_poly)
        return self.figs, self.polys


def run_crystal_voronoi(settings=DEFAULT_SETTINGS):
    interface = ApplicationInterface(parse_commandline=False)
    interface.parse_commandline()
    entry = interface.current_entry
    crystal = entry.crystal
    reps = {"header": PLOTLY_HEADER}

    cv = CrystalVoronoi(
        crystal,
        metal_only=settings["metal_only"],
        weighting=settings["weighting"],
        radius=float(settings["radius"]),
    )
    cv.calculate_voronoi()
    cv.generate_plot(
        plot_fraction=settings["plot_fraction"],
        background=settings["background"],
        opacity=settings["opacity"],
    )

    reps.update(
        {str(i): x.to_html(include_plotlyjs="none") for i, x in enumerate(cv.figs)}
    )

    reps[
        "caption"
    ] = f""" Voronoi with settings: {str(settings)} <br>
    Dr. Chris Kingsbury (ckingsbury@ccdc.cam.ac.uk) (2023)"""

    with interface.html_report(
            title=f"Ligand Report for {interface.identifier}"
    ) as report:
        report.write("<br>".join(reps.values()))
        return report


def run_voronoi_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "idents",
        type=str,
        nargs="+",
        help="the files/refcodes that will be analysed (or .gcd)",
    )
    parser.add_argument(
        "-m", "--maxhits", type=int, default=10000, help="number of hits"
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=f"{getcwd()}/voronoi.pkl",
        help="output location / filetype (.pkl only)",
    )

    parser.add_argument("-v", "--verbose", action="store_true", help="print as you go")

    parser.add_argument(
        "-w",
        "--weighting",
        type=str,
        default="equal",
        help="weighting scheme (vdw, equal, empirical, calculated)",
    )

    parser.add_argument(
        "-c",
        "--columns",
        type=str,
        default="short",
        help="output columns (short, all)",
    )
    parser.add_argument(
        "-r",
        "--radius",
        type=float,
        default=10.0,
        help="maximum radius for interactions",
    )
    parser.add_argument(
        "-mo",
        "--metal_only",
        action="store_true",
        help="metal-type interactions (for magnets etc.)",
    )

    args = parser.parse_args()

    if "." not in "".join(args.idents):
        csd_reader = ccdc.io.CrystalReader("CSD")

    if args.idents[0] == "cwd":
        directory = Path(getcwd())
        args.idents = [
            str(directory / filenm)
            for filenm in list(directory.iterdir())
            if filenm.suffix in [".cif", ".mol2"]
        ]

    new_weighting = weighting_from_name(args.weighting)

    if args.verbose:
        print(
            f"weighting: {args.weighting if (args.weighting in list(WEIGHTING.keys())) else 'equal'}"
        )

    if args.idents[0].endswith(".gcd"):
        csd_reader = ccdc.io.CrystalReader("CSD")
        refcodes = open(args.idents[0], "r").readlines()
        args.idents = [
            x.rstrip("\n") for x in refcodes if len(x.rstrip("\n")) in [6, 8]
        ]

    polyhedra = []
    contacts = []
    cvs = []

    for index, ident in enumerate(args.idents):
        if args.verbose:
            print(f"{ident} ({index + 1} / {len(args.idents)})")
        if "." in ident:
            source = ccdc.io.CrystalReader(ident)
        else:
            source = csd_reader.crystal(ident.upper())
            try:
                [
                    x.assign_bond_types(which="unknown")
                    for x in source.molecule.components
                ]
            except RuntimeError:
                print("error with bond type assigner")

        try:
            cv = CrystalVoronoi(
                source,
                radius=args.radius,
                weighting=new_weighting,
                metal_only=args.metal_only,
            )
            cv.calculate_voronoi()
            polyhedra.append(cv.pyv)
            contacts.append(cv.contacts)
            cvs.append(cv)
        except (ValueError, IndexError):
            if args.verbose:
                print(f"skipping {ident} (likely atoms without coordinates)")

    df = pd.DataFrame(
        np.vstack(contacts), columns=["atom1", "atom2", "distance", "area"]
    )
    df[["distance", "area"]] = df[["distance", "area"]].astype(float)
    df["pair"] = [
        f"{x}-{y}" if x > y else f"{y}-{x}" for x, y in df[["atom1", "atom2"]].values
    ]

    if args.output.endswith(".pkl"):
        df.to_pickle(args.output)

    # elif args.output.endswith(".html"):
    #    [cv.generate_plot() for cv in cvs]


if __name__ == "__main__":
    # For running script from within Mercury
    if len(sys.argv) > 3 and sys.argv[3].endswith(".m2a"):
        import warnings

        warnings.filterwarnings("ignore", category=DeprecationWarning)
        run_crystal_voronoi(DEFAULT_SETTINGS)
    else:
        run_voronoi_cl()

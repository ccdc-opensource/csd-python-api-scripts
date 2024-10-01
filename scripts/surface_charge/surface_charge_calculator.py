#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# The following line states a licence feature that is required to show this script in Mercury and Hermes script menus.
# Created 18/08/2024 by Alex Moldovan
# ccdc-licence-features not-applicable
import os
import warnings
from typing import List, Tuple, Dict

import numpy as np
from ccdc.io import CrystalReader
from ccdc.particle import Surface
from ccdc.utilities import HTMLReport


class SurfaceCharge:
    def __init__(self, crystal, use_existing_charges: bool = False):
        if use_existing_charges:
            # if crystal.molecule.assign_partial_charges() is False:
            #     warnings.warn(f"Gasteiger charges could not be assigned to molecule: {crystal.identifier}",
            #                   RuntimeWarning)
            raise NotImplementedError("Use existing charges instead. Current implementation only supports Gasteiger.")
        self.crystal = crystal
        self.surface = None
        self._surface_charge = None
        self.surface_atom_charge = np.nan
        self.node_projected_charge = np.nan
        self.node_representative_charge = np.nan

    @staticmethod
    def sum_atom_charge(atoms: List[object]) -> float:
        return np.round(np.sum([atom.partial_charge for atom in atoms]), 3)

    def get_average_node_charge(self):
        self.node_charge_dictionary = {}
        node_list = list(self.surface.topology.nodes)
        for node, atoms in self.surface.surface_node_atom_contacts.items():
            node_index = node_list.index(node)
            average_node_charge = 0
            if len(atoms) > 0:
                average_node_charge = self.sum_atom_charge(atoms)
            self.node_charge_dictionary[node_index] = average_node_charge

    def calculate_triangles_properties(self,
                                       tri_index: List[Tuple[int, int, int]]) -> Dict[
        Tuple[int, int, int], Dict[str, float]]:
        surface_area = self.surface.descriptors.surface_area
        self.triangles_properties = {}
        triangle_areas = self.calculate_area_of_triangles(list(self.surface.topology.triangles))
        total_triangle_area = sum(triangle_areas)
        for node_index, triangle_area in zip(tri_index, triangle_areas):
            average_triangle_charge = np.mean([self.node_charge_dictionary[i] for i in node_index])
            triangle_representation = triangle_area / surface_area
            projected_charge = 0
            if np.isclose(triangle_area, 0.0) is False:
                projected_charge = average_triangle_charge / triangle_area
            self.triangles_properties[tuple(node_index)] = {'Average Charge': average_triangle_charge,
                                                            'Triangle Area': triangle_area,
                                                            'Percentage Area': triangle_representation,
                                                            'Node Representative Charge': average_triangle_charge * triangle_representation,
                                                            'Node Projected Surface Charge': projected_charge}

    def calculate_node_charges(self):
        tri_index = self.calculated_node_index_values(list(self.surface.topology.nodes),
                                                      list(self.surface.topology.triangles))
        self.get_average_node_charge()
        self.calculate_triangles_properties(tri_index)
        self.representative_charge = np.sum(
            [triangle['Node Representative Charge'] for triangle in self.triangles_properties.values()])
        self.node_charges = np.sum([triangle['Average Charge'] for triangle in self.triangles_properties.values()])
        return self.representative_charge

    @staticmethod
    def calculate_length(origin: np.ndarray, target: np.ndarray) -> float:
        """Returns distance between target and origin"""
        if not isinstance(origin, np.ndarray) or not isinstance(target, np.ndarray):
            raise TypeError("Please supply numpy arrays for lengths.")
        return np.linalg.norm(target - origin)

    @staticmethod
    def compute_triangle_area(a: float, b: float, c: float) -> float:
        """Calculates area of triangle using Heron's formula"""
        s = (a + b + c) / 2
        return np.sqrt(s * (s - a) * (s - b) * (s - c))

    def calculate_area_of_triangles(self, triangles: List) -> List:
        """ Calculates area of individual triangles from node positions using Heron's formula"""
        triangle_areas = []
        for triangle in triangles:
            pos_0, pos_1, pos_2 = np.array(triangle[0]), np.array(triangle[1]), np.array(triangle[2]),
            a_dist = self.calculate_length(pos_0, pos_1)
            b_dist = self.calculate_length(pos_0, pos_2)
            c_dist = self.calculate_length(pos_1, pos_2)
            triangle_areas.append(self.compute_triangle_area(a_dist, b_dist, c_dist))

        return triangle_areas

    @staticmethod
    def calculated_node_index_values(nodes: List, triangles: List) -> List:
        """
         Calculate index of all triangles for nodes

        :param nodes: list of nodes [x,y,z]
        :param triangles: list of triangles that contain nodes index values
        """
        search_dictionary = {e: i for i, e in enumerate(nodes)}
        return [[search_dictionary[n] for n in tri] for tri in triangles]

    def calculate_surface_charge(self, hkl: Tuple[int, int, int], offset: float):
        self.surface = Surface(self.crystal, miller_indices=hkl, offset=offset)
        if self.surface.slab.assign_partial_charges():
            self.surface_atom_charge = self.sum_atom_charge(atoms=self.surface.surface_atoms)
            self.node_representative_charge = self.calculate_node_charges()
            return
        warnings.warn(f"Gasteiger charges could not be assigned to molecule: {self.crystal.identifier}",
                      RuntimeWarning)


class SurfaceChargeController:
    def __init__(self, structure: str, hkl_and_offsets: List[Tuple[int, int, int, float]],
                 output_directory: str = None, use_existing_charges: bool = False):

        self.surface_node_charges = []
        self.surface_areas = []
        self.surface_node_representative_charge = []
        self.surface_atom_charges = []
        self.surface_charges_per_area = []
        self.surface_charges = []
        self.projected_area = []
        self.crystal = None
        if output_directory is None:
            output_directory = os.getcwd()
        self.output_directory = output_directory
        self.structure = structure
        self._surfaces = None
        self.get_structure()
        self.identifier = self.crystal.identifier
        self._surfaces = hkl_and_offsets
        self.use_existing_charges = use_existing_charges

    def get_structure(self):
        if "." not in self.structure:
            self.crystal = CrystalReader('CSD').crystal(self.structure)
        elif ".mol2" in self.structure:
            self.crystal = CrystalReader(self.structure)[0]
        else:
            raise IOError(" \n ERROR : Please supply refcode or  mol2")

    @property
    def surfaces(self):
        if self._surfaces:
            return self._surfaces

    def calculate_surface_charge(self):
        for surface in self.surfaces:
            charges = SurfaceCharge(crystal=self.crystal, use_existing_charges=self.use_existing_charges)
            charges.calculate_surface_charge(hkl=surface[:3], offset=surface[3])
            self.surface_atom_charges.append(charges.surface_atom_charge)
            self.surface_node_charges.append(charges.node_charges)
            self.surface_node_representative_charge.append(charges.node_representative_charge)

            self.projected_area.append(charges.surface.descriptors.projected_area)
            self.surface_areas.append(charges.surface.descriptors.surface_area)

    def make_report(self):
        html_content = self.generate_html_table()
        output_file = os.path.join(self.output_directory, self.identifier + "_surface_charge.html")
        with HTMLReport(file_name=output_file,
                        ccdc_header=True, append=False, embed_images=False,
                        report_title=f'Surface Charge Calculations - {self.identifier}') as self.report:
            self.report.write(html_content)

            print(f"Results saved to {output_file}")

    def generate_html_table(self):
        # HTML Table Header
        html = """
        <html>
        <head>
            <title>Calculation Results</title>
            <style>
                table { width: 100%; border-collapse: collapse; }
                th, td { border: 1px solid black; padding: 8px; text-align: center; }
                th { background-color: #f2f2f2; }
            </style>
        </head>
        <body>
            <h2>Calculation Results</h2>
            <table>
                <tr>
                    <th width="80px">hkl</th>
                    <th>Offset</th>
                    <th>P<sub>Area</sub>-Projected Area &#8491;<sup>2</sup></th>
                    <th>S<sub>Area</sub>-Surface Area &#8491;<sup>2</sup></th>
                    <th>Total Atom Surface Charge</th>
                    <th>Total Atom Surface Charge/P<sub>A</sub></th>
                    <th>Topological Surface Charge/ S<sub>Area</sub></th>
                </tr>
        """

        # HTML Table Rows
        for i, (h, k, l, offset) in enumerate(self.surfaces):
            hkl = "{" + f"{h}, {k}, {l}" + "}"
            html += f"""
                <tr>
                    <td>{hkl}</td>
                    <td>{offset:.2f}</td>
                    <td>{self.projected_area[i]:.3f}</td>
                    <td>{self.surface_areas[i]:.3f}</td>
                    <td>{self.surface_atom_charges[i]:.3f}</td>
                    <td>{self.surface_atom_charges[i] / self.projected_area[i]:.4f}</td>
                    <td>{self.surface_node_representative_charge[i]:.4f}</td>
                </tr>
            """

        # HTML Table Footer
        html += """
            </table>
            <p><i> *-Surface charge is based on gasteiger partial charges <a href="https://www.sciencedirect.com/science/article/pii/S0040403901949779?via%3Dihub">10.1016/S0040-4039(01)94977-9</a></i> </p>
        </body>
        </html>
        """
        return html

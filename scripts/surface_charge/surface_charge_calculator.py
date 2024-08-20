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
from typing import List, Tuple

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

    def calculate_surface_charge(self, hkl: Tuple[int, int, int], offset: float):
        self.surface = Surface(self.crystal, miller_indices=hkl, offset=offset)
        if self.surface.slab.assign_partial_charges():
            self._surface_charge = np.round(np.sum([atom.partial_charge for atom in self.surface.surface_atoms]), 3)
            return
        warnings.warn(f"Gasteiger charges could not be assigned to molecule: {self.crystal.identifier}",
                      RuntimeWarning)
        self._surface_charge = np.nan

    @property
    def surface_charge(self):
        if self._surface_charge is None:
            raise ValueError("Surface charge calculation not yet calculated, run calculate_surface_charge()")
        return self._surface_charge

    @property
    def surface_charge_per_area(self):
        return self.surface_charge / self.surface.descriptors.projected_area


class SurfaceChargeController:
    def __init__(self, structure: str, hkl_and_offsets: List[Tuple[int, int, int, float]],
                 output_directory: str = None, use_existing_charges: bool = False):
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
            self.surface_charges.append(charges.surface_charge)
            self.projected_area.append(charges.surface.descriptors.projected_area)
            self.surface_charges_per_area.append(charges.surface_charge_per_area)

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
                    <th>hkl</th>
                    <th>Offset</th>
                    <th>Projected Area</th>
                    <th>Surface Charge*</th>
                    <th>Surface Charge per Projected Area</th>
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
                    <td>{self.surface_charges[i]:.3f}</td>
                    <td>{self.surface_charges_per_area[i]:.4f}</td>
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

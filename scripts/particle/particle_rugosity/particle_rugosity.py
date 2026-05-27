#!/usr/bin/env python
#
# This script can be used for any purpose without limitation subject to the
# conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#

# script for calculating the particle rugosity for a given crystal structure (.cif in local working directory, or refcode of a CSD entry)
# requires a CCDC licence, and loading of the CSD-python miniconda environment if run from Linux command line

import argparse
from os.path import isfile

from ccdc.io import CrystalReader
from ccdc.particle import Surface
from ccdc.morphology import BFDHMorphology


def particle_rugosity(crystal_id):
    # crystal is a string, either refcode (eg: 'AABHTZ') or filename (eg: 'AABHTZ.cif')
    # load structure, either (1) crystal structure file in working directory, or
    # (2) CSD refcode
    if isfile(crystal_id):
        try:
            # Crystal structure file
            crystal_1 = CrystalReader(crystal_id)
            crystal = crystal_1[0]
        except:
            raise RuntimeError("Error in reading crystal structure input.\nPlease enter a CSD refcode or provide the path to your crystal structure file")
    else:
        try:
        # CSD refcode
            crystal = CrystalReader('CSD').crystal(crystal_id)
        except:
            raise RuntimeError("Error in reading crystal structure input.\nPlease enter a CSD refcode or provide the path to your crystal structure file")

    print(f"The crystal structure {crystal_id} has been loaded")
    print("Calculating particle rugosity ...")

    # run BFDH morphology to determine predicted particle surface area = hkl planes, and d-spacing of the planes
    morphology = BFDHMorphology(crystal)
    facets = morphology.facets
    all_hkl = [f.miller_indices.hkl for f in facets]
    all_mi = [f.miller_indices for f in facets]
    facets_relA = [morphology.relative_area(mi) for mi in all_mi]
    all_d_hkl = [f.miller_indices.d_spacing for f in facets]
    num_face = len(all_d_hkl)

    # generate surfaces and record rugosity, weigh by particle surface area
    w_part_rug = []
    for i in range(num_face):
        hkl = all_hkl[i]
        rug_list = []

        # check rugosity of crystal plane at 0, 0.25, 0.5, and 0.75 offset from origin 
        for split in range(1, 5):
            os = all_d_hkl[i] / split
            surface = Surface(crystal, hkl, offset=os)
            rug_list.append(surface.descriptors.rugosity)

        # only take the minimum rugosity value (from all offset calculations) for each hkl
        rug = min(rug_list)
        w_surf_rug = rug * facets_relA[i]
        w_part_rug.append(w_surf_rug)

    # sum weighted rugosity values and print total particle rugosity
    print(f'{crystal_id} particle rugosity: {round(sum(w_part_rug), 3)}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("crystal_id")
    args = parser.parse_args()
    crystal = args.crystal_id

    particle_rugosity(crystal)

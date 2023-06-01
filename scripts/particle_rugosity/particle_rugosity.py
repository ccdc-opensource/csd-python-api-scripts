# script for calculating the particle rugosity for a given crystal structure (.cif in local working directory, or refcode of a CSD entry)
# requires a CCDC license, and loading of the CSD-python miniconda environment if run from Linux command line

# import ccdc modules
import ccdc
from ccdc.io import CrystalReader
from ccdc import particle
from ccdc.particle import Surface
from ccdc import morphology
from ccdc.morphology import BFDHMorphology

def particle_rugosity(crystal):
# load structure, either (1) crystal structure file in working directory, or 
# (2) CSD refcode
    head, sep, tail = crystal.partition('.')
    if sep == '.':
    # agrument is crystal structure file 
        xname = crystal
        xtal1 = CrystalReader(xname)
        xtal = xtal1[0]
        xtal1.close
    else:
    # agrument is a CSD refcode 
        xname = crystal
        xtal = CrystalReader('CSD').crystal(xname)

    print(f"The crystal structure {crystal} has been loaded")
    print("Calculating particle rugosity ...")
# run BFDH morphology to determine predicted particle surface area = hkl planes, and d-spacing of the planes
    morphology = BFDHMorphology(xtal)
    facets = morphology.facets
    f = facets[0]
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
        for split in range(1,5):
            os = all_d_hkl[i] / split
            surface = Surface(xtal, hkl, offset=os)
            rug_list.append(surface.descriptors.rugosity)

        # only take the minimum rugosity value (from all offset calculations) for each hkl
        rug = min(rug_list)
        w_surf_rug = rug * facets_relA[i]
        w_part_rug.append(w_surf_rug)

# sum weighted rugosity values and print total particle rugosity
    tot_rug = sum(w_part_rug)
    print(f'{crystal} particle rugosity: ',tot_rug)

# ================================
# argument parser
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input_x")
args = parser.parse_args()
crystal = args.input_x

# run code
particle_rugosity(crystal)


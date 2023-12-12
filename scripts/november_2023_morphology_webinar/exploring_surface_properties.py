#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#

"""
This script shows a few uses of the CSD Particle tools available.
"""
import pandas as pd
from ccdc.io import EntryReader
from ccdc.particle import Surface
from ccdc.morphology import VisualHabit

csd = EntryReader('CSD')
ibuprofen = csd.crystal("IBPRAC18")

visualhabit_settings = VisualHabit.Settings()
visualhabit_settings.potential = 'dreidingII'
vh_results = VisualHabit().calculate(ibuprofen)
vh_morphology = vh_results.morphology

# Surface topology
ibuprofen_100 = Surface(ibuprofen, (1, 0, 0))
print()
print("Ibuprofen (100) topology descriptors (no offset):")
print(ibuprofen_100.descriptors.rugosity)

ibuprofen_100 = Surface(ibuprofen, (1, 0, 0), offset=4.308)
surface_topology = [[(1, 0, 0),
                     ibuprofen_100.descriptors.rugosity,
                     ibuprofen_100.descriptors.rmsd,
                     ibuprofen_100.descriptors.skewness,
                     ibuprofen_100.descriptors.kurtosis]]

print()
print("Ibuprofen (100) topology descriptors:")
print(pd.DataFrame(surface_topology, columns=["Facet", "Rugosity", "RMSD", "Skewness", "Kurtosis"]))

# Surface chemistry
surface_chemistry = [[(1, 0, 0),
                      ibuprofen_100.descriptors.hb_donors,
                      ibuprofen_100.descriptors.hb_acceptors,
                      ibuprofen_100.descriptors.aromatic_bonds]]

print()
print("Ibuprofen (100) surface chemistry descriptors:")
print(pd.DataFrame(surface_chemistry,
                   columns=["Facet", "Density, HB Donors", "Density, HB Acceptors", "Density, Aromatic Bonds"]))

# Average surface properties
facet_rugosity_data = []
for facet in vh_morphology.facets:
    hkl = facet.miller_indices.hkl
    relative_area = round(vh_morphology.relative_area(facet.miller_indices), 3)
    rugosity = Surface(ibuprofen, facet.miller_indices.hkl).descriptors.rugosity
    weighted_rugosity = relative_area * rugosity
    facet_rugosity_data.append([hkl, relative_area, rugosity, weighted_rugosity])

facet_rugosities = pd.DataFrame(facet_rugosity_data,
                                columns=["Facet", "Relative Area", "Rugosity", "Weighted Rugosity"])
print()
print("Facet rugosities:")
print(facet_rugosities)

# Particle rugosity
import morphology_plot as plotter

print()
print("Weighted rugosity:")
print(round(facet_rugosities["Weighted Rugosity"].sum(), 3))
plotter.generate_morphology_plot(vh_morphology, labels=list(facet_rugosities["Rugosity"]))

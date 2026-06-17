#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#

"""
This script runs a morphology calculation and outputs the facet properties.
"""
import pandas as pd
from ccdc.io import EntryReader
from ccdc.morphology import BFDHMorphology

csd = EntryReader('CSD')
ibuprofen = csd.crystal("IBPRAC18")
bfdh_morphology = BFDHMorphology(ibuprofen)

bfdh_facet_data = [[f.miller_indices.hkl, round(f.area, 3)] for f in bfdh_morphology.facets]

print()
print("BFDH Data:")
print(pd.DataFrame(bfdh_facet_data, columns=["Miller Index", "Facet Area"]))

from ccdc.morphology import VisualHabit

visualhabit_settings = VisualHabit.Settings()
visualhabit_settings.potential = 'dreidingII'
vh_results = VisualHabit().calculate(ibuprofen)
vh_morphology = vh_results.morphology

vh_facet_data = [[f.miller_indices.hkl, f.area, f.attachment_energy] for f in vh_morphology.facets]

print()
print("VisualHabit Data:")
print(pd.DataFrame(vh_facet_data, columns=["Miller Index", "Facet Area", "Attachment Energy"]))

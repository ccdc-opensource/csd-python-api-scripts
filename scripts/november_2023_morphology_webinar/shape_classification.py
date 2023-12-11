#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#

"""
Script to

"""

import morphology_shape_classification as classifier
from ccdc.io import EntryReader
from ccdc.morphology import BFDHMorphology
from ccdc.morphology import VisualHabit

csd = EntryReader('CSD')
visualhabit_settings = VisualHabit.Settings()
visualhabit_settings.potential = 'dreidingII'

succinic_acid_bfdh = BFDHMorphology(csd.crystal("SUCACB02"))
succinic_acid_vh_dreiding = VisualHabit().calculate(csd.crystal("SUCACB02")).morphology

bfdh_shape = classifier.ShapeClassifier(succinic_acid_bfdh)
vh_shape = classifier.ShapeClassifier(succinic_acid_vh_dreiding)
print(f"BFDH morphology is classified as a {bfdh_shape.shape_description}")
print(f"VisualHabit morphology is classified as a {vh_shape.shape_description}")

classifier.plot_multiple_shapes([bfdh_shape, vh_shape])

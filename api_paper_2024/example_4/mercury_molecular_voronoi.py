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
# Mercury interface to the metal Voronoi polyhedra which encapsulate atoms
# with weighting schemes available
#

import voronoi

new_settings = voronoi.DEFAULT_SETTINGS.copy()
new_settings.update({
    "radius": 10.0,  # Angstroms; radius for the substructure contacts
    "opacity": 1.0,
    "metal_only": False,
    "weighting": voronoi.WEIGHTING['vdw'],
})

if __name__ == "__main__":
    voronoi.run_crystal_voronoi(new_settings)

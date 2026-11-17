#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2024-07-26 Created by Chris Kingsbury, the Cambridge Crystallographic Data Centre
# ORCID 0000-0002-4694-5566
# 2024-10-30 Modified by Alex Moldovan, the Cambridge Crystallographic Data Centre
# flake8: noqa
import sys
from simplifier import MercuryInterface

import numpy as np

if __name__ == "__main__":
    if len(sys.argv) > 3 and sys.argv[3].endswith(".m2a"):
        mercury_interface = MercuryInterface()
        if not mercury_interface.write_readme_if_none_selected():
            atom_clicking_order = np.argsort(np.argsort(mercury_interface.interface.selected_indices))
            selected_atoms_in_order = [mercury_interface.interface.selected_atoms[x] for x in atom_clicking_order]
            mercury_interface.simplifier.define_simplification(crystal=mercury_interface.crystal,
                                                               user_selected_atoms=selected_atoms_in_order)
            mercury_interface.simplifier.use_simplification(mercury_interface.crystal)
            mercury_interface.write_output()
    else:
        raise NotImplementedError("Implementation is for Mercury - CSD Python API menu only")

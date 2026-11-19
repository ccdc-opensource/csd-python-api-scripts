#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2023-08-18 Created by Chris Kingsbury, the Cambridge Crystallographic Data Centre
# ORCID 0000-0002-4694-5566
# flake8: noqa
import sys

from simplifier import MercuryInterface

if __name__ == "__main__":
    if len(sys.argv) > 3 and sys.argv[3].endswith(".m2a"):
        mercury_interface = MercuryInterface(leave_unmatched_molecules=False)
        mercury_interface.simplifier.use_simplification(crystal=mercury_interface.crystal)
        mercury_interface.write_output()
    else:
        raise NotImplementedError("Implementation is for Mercury - CSD Python API menu only")

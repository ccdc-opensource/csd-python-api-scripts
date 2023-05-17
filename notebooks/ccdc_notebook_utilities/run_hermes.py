#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2022-05-17: created by Jason Cole, The Cambridge Crystallographic Data Centre
#

from ccdc.io import csd_directory
from pathlib import Path
from platform import platform
from subprocess import Popen

"""
From inside a notebook, run CCDC Hermes utility or fail if we cant find it. Assumes the
software is installed alongside the data folder currently.
"""
def run_hermes(*filenames):

    try:
        hermes_dir = Path(csd_directory()) / '..' / '..' / 'ccdc-software' / 'hermes'
        hermes_exe = (hermes_dir / 'hermes.exe' if platform().startswith('Windows') else hermes_dir / 'hermes').as_posix()
        _ = Popen([hermes_exe, *filenames], creationflags=0x00000008)
    except Exception as e:
        print( f"Couldnt run Hermes {e}")

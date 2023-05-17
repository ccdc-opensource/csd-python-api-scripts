#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2022-05-17: created by Jason Cole, The Cambridge Crystallographic Data Centre
#
import logging
from platform import platform
import sys
import ccdc.io
import os


def create_logger(verbose=True):
    """
    From inside a notebook, create a logger and log starting information
    """
    
    logger = logging.getLogger(__name__)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('[%(asctime)s %(levelname)-7s] %(message)s', datefmt='%y-%m-%d %H:%M:%S'))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    if verbose:
        logger.info(f"""
Platform:                     {platform()}

Python exe:                   {sys.executable}
Python version:               {'.'.join(str(x) for x in sys.version_info[:3])}

CSD version:                  {ccdc.io.csd_version()}
CSD directory:                {ccdc.io.csd_directory()}
API version:                  {ccdc.__version__}

CSDHOME:                      {os.environ.get('CSDHOME', 'Not set')}
CCDC_LICENSING_CONFIGURATION: {os.environ.get('CCDC_LICENSING_CONFIGURATION', 'Not set')}
""")
    return logger

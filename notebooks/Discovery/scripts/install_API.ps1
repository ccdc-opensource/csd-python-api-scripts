##############################################################################################################
# 
# This document describes how to set up the CSD Python API for use with the Discovery Notebooks.
# 
# This version assumes Windows PowerShell is being used; a Linux/MacOS version is available on request. 
# 
# For the official installation instructions, see:
#
#     https://downloads.ccdc.cam.ac.uk/documentation/API/installation_notes.html
# 
# As noted in that document, there are two ways of accessing the CSD Python API. One can either use the Miniconda 
# distribution included with the CSD release (e.g. in 'C:\Program Files\CCDC\Python_API_2020\miniconda\')
# or the API package can be installed into a separate Python distribution. I recommend the second approach,
# as the release distribution is intended for the use of CSD software and there can be version compatibility
# issues when installing extra modules into an environment derived from the release distribution.
# 
# In the instructions below, a full Anaconda distribution could be substituted for Miniconda. I prefer
# building up an install from Miniconda, however, as it gives a better picture of what is actually necessary.
# 
# Python distributions other than Miniconda/Anaconda can also be used as the API is also available as a
# pip package. However, as I have always found the Miniconda-based solution to be preferable in practice I
# don't have much familiarity with the alternatives. If you wish to use a different approach, please refer
# to the instructions document linked above.
# 
# If you do not have Miniconda (or Anaconda, if preferred) installed already, please download and install it
# before attempting to follow the instructions below. Note that the CCDC now only supports Python 3, so please
# ensure you choose the right version.
# Download the latest version of Miniconda from here: https://repo.anaconda.com/miniconda
# Follow the instructions here: https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html
#
# The CSD Python API conda channel is available via the downloads site, using your normal mechanism:
# 
#     https://www.ccdc.cam.ac.uk/support-and-resources/csdsdownloads/
# 
# If you are a CCDC Research Partner, sign in to the downloads site, click on the link 'CSD Python API'
# (under 'CSD-System') and you will see a list of the packages you can download. If, instead, you have a
# Customer Number and Activation Key, go to the site and request a set of links to be emailed to you.
# Either way, you will need the main CSD Python API conda package. This must be the appropriate one
# for your Operating System, i.e. Windows, Linux or macOS. Note that the name might differ slightly
# depending on which source you get them from, but they will have a common prefix. For example, on
# the RP site, for Windows it will be ‘Python 3.7 CSD Python API 3.0.3 conda package (Windows 64-bit)’;
# in a downloads email it would be ‘Python 3.7 CSD Python API 3.0.3 Windows 64-bit conda installer’.
# 
##############################################################################################################

# Name for the environment (must not exist)...

$env_name = "ccdc"

# Miniconda location (must exist)...
# This was the default name and location offered by the installer, and may need altering if you used
# an alternative location or if the defaults have changed. Please check it exists before proceeding.

$miniconda = "${HOME}\miniconda"  # "${HOME}\miniconda3"

# Path to CSD Python API conda channel archive (must exist)...
# This assumes you have manually saved the archive into your personal downloads folder.
# In addition, the filename may need updating to that of the particular version you have downloaded.

$conda_channel_archive = "${HOME}\Downloads\csd-python-api-3.0.4-win32-64-py3.7-conda.zip"

##############################################################################################################

# Check the Miniconda dir exists...

If ( -Not (Test-Path $miniconda) ) {  

    Write-Warning -WarningAction Stop -Message "Warning! Miniconda dir '${miniconda}' not found!" 
}

# Check conda channel archive file exists...

If ( -Not (Test-Path $conda_channel_archive) ) {  

    Write-Warning -WarningAction Stop -Message "Warning! Conda channel archive '${conda_channel_archive}' not found!" 
}

##############################################################################################################

# At this stage, conda might not be accessible if the base environment isn't activated... 

Try {
    
    Get-Command conda -ErrorAction Stop
}
Catch {  # Activate the base environment...
    
    . ${miniconda}\shell\condabin\conda-hook.ps1
    
    conda activate    
}

# The base conda environemnt should now be activated...

Get-Command conda -ErrorAction Stop

######

# Check whether the environment already exists...

If ( $(conda env list | Where { $_.Split()[0] -Eq $env_name }) -Ne $null ) {

    Write-Warning -WarningAction Stop -Message "Warning! Conda channel '${env_name}' already exists!" 

}

#######

# Deactivate any currently active environment...

conda deactivate

Clear-Host

##############################################################################################################

# Create and activate a new conda environment (note that the CCDC now only supports Python 3.7)...

conda create --yes --name $env_name python=3.7

conda activate $env_name

######

# Unzip the conda channel archive into the currect directory, such that directory .\ccdc_conda_channel is created.


7z x $conda_channel_archive

If ( -Not (Test-Path .\ccdc_conda_channel) ) {  # Check local conda channel exists

    Write-Warning -WarningAction Stop -Message "Warning! Expected directory 'ccdc_conda_channel' not found!" 
}

######

# Install the CSD Python API into the active environment...

conda install --yes --channel $(Resolve-Path .\ccdc_conda_channel) csd-python-api  # N.B. Conda requires an absolute path

# Check that the CSD API was installed correctly...

python -c 'import ccdc; print(ccdc.__version__)'

If ($LastExitCode -ne 0) { 

    Write-Warning -WarningAction Stop -Message "Warning! The CSD Python API module does not seem to have been installed properly." 
}

# Clean up...

Remove-Item -Recurse .\ccdc_conda_channel

######

# Install RDKit (N.B. this will bring in Pandas as a dependency)...
# N.B. If you have the RDBase environment variable set, please unset it as it can cause problems here.

conda install --yes --channel conda-forge rdkit

# Check RDKit was installed properly...

python -c 'import rdkit; print(rdkit.__version__)'

If ($LastExitCode -ne 0) { 

    Write-Warning -WarningAction Stop -Message "Warning! The RDKit module does not seem to have been installed properly." 
}

######

# Install Jupyter-lab...

conda install --yes jupyterlab 

#######

# At time of writing I have been experiencing two separate problems with Jupyter that necessitate downgrading two packages.
# 
# The first involves a DeprecationWarning being issued whenever a cell was executed. It is described in the following SO link:
# https://stackoverflow.com/questions/63413807/deprecation-warning-from-jupyter-should-run-async-will-not-call-transform-c

conda install --yes ipython=7.10.1

# The second issue involved autocompletion being broken. It is describe in the following GitHub link:
# https://github.com/ipython/ipython/issues/12740

conda install --yes jedi=0.17.2

######

# Install other required packages...

conda install --yes --channel conda-forge altair #  nbstripout

######

# Install pymol...
# N.B. this conda version of PyMOL can coexist with other installations and is convenient for the notebooks

conda install --yes --channel schrodinger pymol

############

# Start Jupyter Lab...

# Start-Process jupyter-lab

##############################################################################################################
# End
##############################################################################################################
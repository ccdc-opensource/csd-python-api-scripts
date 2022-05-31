##########################################################################################################################
# 
# visualise.py
# 
# This is a helper script descigned to be used alongside the Jupyter notebook Interaction_maps.ipynb
# 
##########################################################################################################################

from pymol import cmd
from pathlib import Path

##########################################################################################################################

cmd.do("reinitialize everything")

# Settings...

cmd.set('stick_radius', 0.125)
cmd.set('transparency', 0.5)
cmd.set('valence', 1)
cmd.set('sphere_scale', 0.1)

######

# Not sure how to pass parameters to a PyMOL script run this way, so do it via a file...

with Path('results_dir.txt').open('r') as fp:

    results_dir = Path(fp.read().strip())

######

# Load protein, if present...

protein_file = results_dir / '2UW7.mol2'

if protein_file.exists():

    cmd.load(protein_file, 'protein')
    
    # util.cbag('protein')

# Load ligand...

ligand_file = results_dir / 'A_GVP1351.mol2'

cmd.load(ligand_file, 'ligand')

cmd.show('sticks', 'ligand') 

# util.cbaw('ligand')

######

cmd.hide('everything', '(h. and (e. c extend 1))', ) # Hide non-polar Hs

##########################################################################################################################

probe_names = ['Uncharged NH Nitrogen', 'Carbonyl Oxygen', 'Aromatic CH Carbon']
colours     = ['blue',                  'red',             'yellow']

for probe_name, colour in zip(probe_names, colours):

    name = probe_name.replace(' ', '_')
    
    # Load grid...
    
    grid_file = results_dir / (name + '.grd')
    
    cmd.load(grid_file, name)
    
    surface_name = f'{name}_surface'
    
    cmd.isosurface(surface_name, name, 1.5)
    
    cmd.color(colour, surface_name)
    
    ######
    
    # Load hotspots (map maxima or 'peaks')...
    
    peaks_file  = results_dir / (name + '.peaks.mol2')
    
    peaks_name = f'{name}_peaks'
    
    cmd.load(peaks_file, peaks_name)
    
    cmd.hide('everything', peaks_name)
    
    cmd.show('sphere', peaks_name)
    
    cmd.label(peaks_name, "'%.2f' % partial_charge")
    
    cmd.color(colour, peaks_name)

    cmd.disable(peaks_name)  # Hide by default

##########################################################################################################################

cmd.orient('ligand')

cmd.save(str(results_dir / 'visualisation.pse'))

##########################################################################################################################
# End
##########################################################################################################################
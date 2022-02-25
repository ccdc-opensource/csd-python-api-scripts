#find_binding_conformation.py 

We know that most pharmaceutically relevant compounds bind to their targets in a relaxed conformation. The challenge in discovery is to figure out rapidly which conformations are readily accessible for the molecules we are considering. There is now a new solution to address this based on statistical, rather than just energetic approaches.

Driven by the wealth and diversity of bond, angle and torsion information in the Cambridge Structural Database (CSD), the CSD Conformer Generator produces realistic ensembles of low energy ligand structures. These are ready to be exploited for drug design in the presence and also in the absence of detailed knowledge about the three-dimensional structure of the protein active site.

Starting from a list of PDB-codes, this script generates idealized conformers
for ligands and evaluates their RMSD to the conformation in the PDB.

#Author
'Brandl, Giangreco, Higueruelo, Schaerfer and Sykes'

#Requirements
This script uses the CSD Python API and needs full installation of CSD. 
This script also uses PDBe's and RCSB's API to obtain PDB related information.

#Usage

```Python
python find_binding_conformation.py pdb_example.txt
```

#Output

Subdirectories for each PDB entry with the conformers generated for each ligand, and a spreadsheet (.csv) with the results of the comparison.


For more details about conformer generator see this blog https://www.ccdc.cam.ac.uk/Community/blog/post-61/

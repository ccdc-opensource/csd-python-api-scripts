The target here is SYK tyrosine kinase (5LMA).

The ligands in 'input.mol2' are built from the SMILES in 'input.csv'. If the name is a PDBe code, it means the SMILES
corresponds to the crystallographic ligand from that structure (with conventional ionissation states assigned).

Note that not all ligands can be cross-docked, as there is an induced-fit effect in SYK that GOLD cannot reproduce.
If the name has a suffix, the SMILES is a manually-generated analogue. The suffix '_bad' means the structure is designed
to fail, even in the parent crystal structure (this is to e.g. illustrate constrain-violation penalties). 
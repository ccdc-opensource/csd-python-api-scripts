The notebooks in this directory illustrate the use of the CSD Docking API.

The material in gold_multi.zip illustrates the use of the standard Python multiprocessing module
and the Docking API to parallelize GOLD docking at a multicore workstation level. 

For very large-scale docking, please enquire about our GOLD Cloud or GOLD Cluster tools.

Material about covalent docking with GOLD is also available on request.

Note on the input files provided
--------------------------------

The example target provided here is SYK tyrosine kinase (5LMA).

A small number of example ligands are provided as SMILES in input_files/input.csv. If the name is a PDB code, it means
the SMILES correspondes to the crystallographic ligand from that structure (with conventional ionization states assigned).
If the name has a suffix, the SMILES is a manually-generated analogue. Note that not all ligands can be cross-docked into
5LMA, as there are induced-fit effects in SYK that GOLD cannot reproduce.
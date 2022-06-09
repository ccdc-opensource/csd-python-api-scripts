# Discovery Notebooks

The notebooks in this directory will be of interest primarily to scientists working in Drug Discovery. They show how the API may be used to complement the applications in the [CSD-Discovery](https://www.ccdc.cam.ac.uk/Solutions/csd-discovery/) suite.

## Contents

Directory                                                | Contents
-----                                                    | -----
[00_Background](00_Background)                           | Basic introduction to API concepts
[01_CSD_Search](01_CSD_Search)                           | Searching the CSD and extracting geometrical parameters (as in ConQuest)
[02_Protein_Ligand](02_Protein_Ligand)                   | Searching the CrossMiner database of protein binding sites and extracting geometrical parameters
[03_Molecular_geometries](03_Molecular_geometries)       | Indentify unusal intramolecular geometries uisng CSD data (as in Mogul)
[04_Conformer_generation](04_Conformer_generation)       | Generating conformers using parameter distributions derived from the CSD 
[05_Molecular_interactions](05_Molecular_interactions)   | Investigate intermolecular interactions (as in IsoStar)
[06_Interaction_maps](06_Interaction_maps)               | Investigate intermolecular interactions (as in SuperStar)
[07_Cavities](07_Cavities)                               | Search and compare protein cavities
[08_Docking](08_Docking)                                 | Programmatic GOLD docking
[09_Covalent_Docking](09_Covalent_Docking)               | Covalent docking with GOLD
[10_Editing_molecules](10_Editing_molecules)             | Creating molecules from scratch

## Requirements

Beyond the API, the only Python modules required are [Pandas](https://pandas.pydata.org/) and [Plotly](https://plotly.com/); for the 
Covalent Docking notebooks, [RDKit](https://rdkit.org/) is also required. These may be installed from [conda-forge](https://conda-forge.org/):

```
conda install --yes --channel=conda-forge pandas plotly rdkit
```
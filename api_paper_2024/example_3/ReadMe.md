# Similarity-driven docking

A key feature used in this script is GOLD in interactive docking mode;
this allows a user to set up a GOLD daemon object on a socket and send individual molecules to the socket for docking.
This avoids repeated initialisation when used to receive ligands from external data sources.
The ChEMBL Python API (Davies et al., 2015) is used to facilitate similarity searching in the ChEMBL database.

## Dependencies

Optional:

- chembl_webresource_client

## To Run
```bash
python similarity_docking.py Example/unconstrained.conf
```

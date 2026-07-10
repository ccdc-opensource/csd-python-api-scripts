# Disorder Workflow Demo

## Summary

Exports a CSD entry to a manually assembled CIF that preserves atom-site occupancy,
disorder assembly/group annotations, anisotropic displacement parameters, and bond
connectivity.

This example demonstrates a workflow that is not covered by the built-in CIF writers:
reading a disordered entry from the CSD Python API and reconstructing a CIF with
explicit `_atom_site_occupancy`, `_atom_site_disorder_assembly`, and
`_atom_site_disorder_group` fields.

Example output:

```text
Licence requirement: CSD-Materials or another licence tier with CSD Python API database access
Exported ABACIR -> output/ABACIR.cif (atoms=148, partial_occ=16, bonds=154)
Exported ABABUB -> output/ABABUB.cif (atoms=84, partial_occ=8, bonds=88)
```

## Requirements

- CSD Python API with database access
- Python 3.10+

## Licensing Requirements

CSD-Materials or another licence tier that provides access to the CSD Python API and CSD database.

## Instructions on Running

```bash
python demo_csd_disorder_workflow.py --refcodes ABACIR ABABUB --output-dir ./output
```

Arguments:

- `--refcodes`: one or more CSD refcodes
- `--output-dir`: destination directory for exported CIF files

## Author

Ming-Yu Guo (2026-06-29)

> For feedback or to report any issues please contact [support@ccdc.cam.ac.uk](mailto:support@ccdc.cam.ac.uk)

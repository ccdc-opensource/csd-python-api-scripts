## Conformer Demo

This is a short script to generate conformers with some rudimentary analysis for a single molecule.
There are also options to overlay the results to view in Hermes.


### Example output showing what the user can expect to see:
```
Reading file: AZD9291.mol2 ... done.
Generating conformers, maximum of 20 ... done, generated 20 conformers.
Sampling limit reached? No.
How many rotamers had no observations? 0.
Normalised score of most probable conformer: 0.0.
Most probable conformer RMSD wrt input: 3.276; wrt minimised: 3.202.
Scores of top 10 conformers: 0.000, 0.000, 0.000, 0.027, 0.027, 0.027, 0.027, 0.027, 0.029, 0.029.
Overlaying conformers ... done.
Writing file superimposed ... done.
```

CCDC Python API Licence required, minimum version: 3.0.15


There is an accompanying mol2 file with this script, but users may use any small molecule provided in a file format readable by our API (e.g. mol, mol2, sdf, etc)


Author: Chris Ringrose - 22/11/24


    For feedback or to report any issues please contact support@ccdc.cam.ac.uk

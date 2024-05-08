# Voronoi

## Visualising metal-metal interactions using Voronoi Tessellation

CCDC interface to Voronoi polyhedra which encapsulate the atomic / metal domains

## Dependencies

- `plotly`
- `pyvoro`

## To Run

For general use, it is best to run mercury_molecular_voronoi.py through Mercury -
this script will output an HTML with the Voronoi graph. More complicated queries can
be achieved with the notebook (see notebooks/voronoi) or command line as below.

```cmd
usage: voronoi.py [-h] [-m MAXHITS] [-o OUTPUT] [-v] [-w WEIGHTING] [-c COLUMNS]
[-r RADIUS] [-mo]
idents [idents ...]

positional arguments:
idents the files/refcodes that will be analysed (or .gcd)

optional arguments:
-h, --help show this help message and exit
-m MAXHITS, --maxhits MAXHITS
number of hits
-o OUTPUT, --output OUTPUT
output location / filetype (.pkl only)
-v, --verbose print as you go
-w WEIGHTING, --weighting WEIGHTING
weighting scheme (vdw, equal, empirical, calculated)
-c COLUMNS, --columns COLUMNS
output columns (short, all)
-r RADIUS, --radius RADIUS
maximum radius for interactions
-mo, --metal_only metal-type interactions (for magnets etc.)
```

## Author

Chris Kingsbury (2024)
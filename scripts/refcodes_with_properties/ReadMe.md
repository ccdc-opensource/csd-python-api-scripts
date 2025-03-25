# Refcode List Generator

## Summary

A script that allows you to create refcode lists (or CSV files of properties for a refcode list) for simple properties.
The advantage of the script is that the control is via an easy to read file so you can keep an interprettable record of
how a test set was generated in research. You can also then reproduce the list, or indeed run it on a new database and
update it with the same conditions.

### Relevance

We want research to be FAIR (Findable, Attributable, Interoperable and Reproducible) - this script means we can create a
simple description of the test set used that any researcher could then reproduce from the script and the description.

## Requirements

- Tested with CSD Python API version 3.9 on Linux and Windows
- ccdc.io
- ccdc.search

## Licensing Requirements

- CSD-Core

## Instructions on Running

### Linux command line

- load the CSD Python API Miniconda environment
- create a text control file with the various control lines specified
- call Python to read the script and specify necessary arguments

~~~bash
python refcodes_with_properties.py --help
~~~

The above will print an extended help message that describes the registered

You can run the script with an Example file. Results are printed by default and can be redirected to be saved in an
output file, e.g.

~~~
python refcodes_with_properties.py -c example_control_file.txt -o mylist.gcd
~~~

This will generate a GCD file that can be used in other work.

### Windows CSD Python API

- launch a CMD window
- Use the installed version of the CSD Python API, for example C:\Users\<YOUR WINDOWS USERNAME>
  \CCDC\ccdc-software\csd-python-api assuming the CCDC tools are installed in the ususal place do this

~~~bat
C:\Users\<YOUR WINDOWS USERNAME>\CCDC\ccdc-software\csd-python-api\run_python_api.bat refcodes_with_properties.py --help
~~~

## Author

_Jason C.Cole_ 2025

> For feedback or to report any issues please contact [support@ccdc.cam.ac.uk](mailto:support@ccdc.cam.ac.uk)

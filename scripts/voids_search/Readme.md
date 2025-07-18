# Void Search

Script for searching the structures in the CSD for void properties pre-calculated with PoreAnalyser. Can be run on the command line with:

```python
python void_search.py
```

or by launching from the CSD Python API dropdown in Mercury. If Launching from Mercury, put the script, with the accompanying search_result_template.html file, in your /mercury/scripts folder,
or use the options dialog to "Add Location" to where you have saved the script. When the script is run, it will launch a dailog box with properties to search on, and spaces to fill in queries.
These query boxes accept values like "200-300" or ">5" to constrain the search. A .tsv file is written with the search results, and the table is displayed in the Data Analysis tab; non-Ascii
characters in chemical names will be deleted.
A .gcd file of refcodes in the search results is also written, and read into Mercury

## Author

Chris Kingsbury <kingsbury@ccdc.cam.ac.uk>

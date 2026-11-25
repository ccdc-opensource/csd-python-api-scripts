#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2015-04-02: created by Clare Macrae, the Cambridge Crystallographic Data Centre
# 2015-06-18: made available by the Cambridge Crystallographic Data Centre
# 2025-04-17: Modified by Chris Kingsbury
# 2025-07-02: Modified from SMARTS_QUERY to query pore / void params
# 2025-09-15: Modified to search Topology information.py
# 2026-02-23: Added SearchDialog2 for hover-over hints and more info
# 2026-11-19: adapted from above for searching topology data fields

"""
search_topology_csd.py - return results from the calculated topology using CSD internals
See ../topology for more information
"""
from pathlib import Path
import search_tools
import sys

import ccdc.utilities
import ccdc.io
import search_tools
import sqlite3
from pathlib import Path
import pandas as pd

DATABASE_PATHS = sorted(Path(ccdc.io.csd_directory()).glob("*.sqlite"))

def string_to_query(variable: str, table: str, input: str):
    """
    converts a string input to a valid SQL subquery
    """
    if not input:
        return ""
    if variable == "MaxDimensionality":
        spec = search_tools.query_string_to_tuple(input, 3)
        return f"{table}.{variable} BETWEEN {spec[0]} AND {spec[1]}"
    return f"{table}.{variable} LIKE '%{input}%'"

def perform_csd_topology_search(
        query_dict: dict = {},
        interface: ccdc.utilities.ApplicationInterface = None,
        database_paths: list = None
        ):
    if not database_paths:
        database_paths = DATABASE_PATHS
    columns = ", ".join([
        "IDENTIFIER.EntryIdentifier",
        "TOPOLOGY_PROPERTIES.VertexSymbol",
        "TOPOLOGY_PROPERTIES.UniqueIdentifier",
        "TOPOLOGY_PROPERTIES.RcsrAssignment",
        "TOPOLOGY_PROPERTIES.MaxDimensionality"])
    table = "TOPOLOGY_PROPERTIES INNER JOIN IDENTIFIER ON IDENTIFIER.Id == TOPOLOGY_PROPERTIES.EntryId"
    subqueries = [string_to_query(variable, "TOPOLOGY_PROPERTIES", input) for variable, input in query_dict.items() if input]
    query_where = " AND ".join(subqueries)
    query = f"""SELECT {columns} FROM {table} WHERE {query_where}"""
    results = []

    for path in DATABASE_PATHS:
        conn = sqlite3.connect(path)
        cur = conn.execute("PRAGMA table_list")
        table_names = [x[1] for x in cur.fetchall()]
        has_topology_information = "TOPOLOGY_PROPERTIES" in table_names
        if not has_topology_information:
            continue
        cur = conn.execute(query)
        results.append(cur.fetchall())

    tsv_columns = ["Identifier", "Vertex Symbol", "Unique Identifier", "RCSR Assignment", "Max Dimensionality"]
    df = pd.DataFrame([a for b in results for a in b], columns = tsv_columns)
    if interface:
        search_tools.write_tsv_from_dataframe(df, interface.output_tsv_file)
        search_tools.write_search_result_html(
        title="Search",
        data=f"""
        Query: {query_dict} ({query}) <br>
        Result:{len(df)} hits in {len(df["Identifier"].unique())} structures in {len(results)} databases
        """,
        path=interface.output_html_file)
    return df

test_query_dict = {"MaxDimensionality": "3"}

perform_csd_topology_search(test_query_dict)

CSD_QUERY_DICT = {
      "VertexSymbol":{
            "label": "Vertex Symbol",
            "header": "VertexSymbol",
            "hovertext": """
            A representation of the network as the vertex symbol, in
 HTML format - e.g. '4(4<sup>4</sup>.8<sub>2</sub><sup>2</sup>)'
 is the layers/sql network"""
            },
      "RcsrAssignment":{
            "label": "RCSR Assignment",
            "header": "RcsrAssignment",
            "hovertext": "The exact three letter code from the RCSR list \n (where available)"
      },
      "MaxDimensionality":{
            "label": "Dimensionality",
            "header": "MaxDimensionality",
            "hovertext": "The number of dimensions that are spanned by the \n connected network"
      },
      "UniqueIdentifier":{
            "label": "Unique Identifier",
            "header": "UniqueIdentifier",
            "hovertext": "a generated 8 letter code from SHA1 hash of the \n vertex and cs10 representations"
      }
      }

def main():
    search_info = search_tools.get_query_text_with_hovertext(CSD_QUERY_DICT)
    if len(sys.argv) > 3 and sys.argv[3].endswith(".m2a"):
        perform_csd_topology_search(search_info, interface = ccdc.utilities.ApplicationInterface())
    else:
        perform_csd_topology_search(search_info)

if __name__ == "__main__":
    main()

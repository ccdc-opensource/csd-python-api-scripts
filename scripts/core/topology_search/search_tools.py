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
#
# ccdc-licence-features not-applicable


"""search_tools.py - simple way to make a search UI"""

import sys
import os

if sys.platform.startswith("darwin"):
    os.environ["CCDC_PYTHON_API_NO_QAPPLICATION"] = "1"

import ccdc.utilities

from jinja2 import Template
from pathlib import Path
import sqlite3
import pandas as pd

import sys

import os
import tkinter
import tkinter.font
from tkinter import ttk

if sys.platform.startswith("darwin"):
    os.environ["CCDC_PYTHON_API_NO_QAPPLICATION"] = "1"

RESULT_TEMPLATE_PATH = Path(__file__).parent / "templates" / "search_result_template.html"

STANDARD_SEARCH_ERROR_MESSAGE = """Nothing was given to search for.
    Please enter either one number (to search +/- 5 %) or two numbers
    separated by a "-" in numeric search boxes, or a string into the string boxes"""


def write_search_result_html(title, data, path):
    with open(path, "w") as report:
        tl = Template(
            open(
                RESULT_TEMPLATE_PATH,
                "r",
            ).read()
        )
        report.write(
            tl.render(
                title=title,
                data=data
            )
        )


def write_tsv_from_dataframe(
        df: pd.DataFrame,
        path: str):
    """
    uses ccdc.utilities.CSVWriter to write a tsv from a pandas dataframe
    """
    with ccdc.utilities.CSVWriter(
        path, header=df.columns, delimiter="\t"
    ) as tsv_file:
        for row in df.values:
            tsv_file.write_row(
                [
                    str(value) if (str(value).lower() not in ("nan")) else "0"
                    for value in row
                ]
            )


def run_sqlite_search(query_dict: dict,
                      sqlite_path: Path,
                      interface: ccdc.utilities.ApplicationInterface = None,
                      settings: dict = {}):
    """Search an sqlite file for descriptors set out in the query_dict.
    This will generate a list of structures as a .tsv and a summary as
    .html

    query_dict: dict
        dictionary describing a search to be performed, as in topology_query_dict.json
        with keys describing the different search fields (e.g. for dimensionality)
        and the values denoting the search type (numeric, etc) the table name and the
        sql code to run this search
    sqlite_path: Path, list[Path]
        the path to the sqlite file or files to be queried
    interface: ccdc.utilities.ApplicationInterface
        an interface to print results to Mercury
    settings: dict = {}
        contains title, search type and caption

    returns : None
    """

    if not interface:
        # Create the ApplicationInterface without immediately parsing command line parameters
        # so that we can add a custom parameter to pass the search query in directly.
        interface = ccdc.utilities.ApplicationInterface(parse_commandline=False)

        # Add a custom command-line parameter for the compound name.
        for k, v in query_dict.items():
            interface.commandline_parser.add_argument(
                f"--{k}",
                required=False,
                default="",
                help=f"{v} to search for, bypassing the GUI dialog.",
            )

        # Parse the commandline including checking for the compound parameter.
        interface.parse_commandline()

    filtered_query_dict = {k: v for k, v in query_dict.items() if v.get("show_in_search", "True") == "True"}
    # open a GUI dialog asking for the search query.

    if settings.get("hovertext", False):
        search_dict = get_query_text_with_hovertext(filtered_query_dict)
    else:
        search_dict = get_query_text(filtered_query_dict)
    # Error out if the compound name given is empty
    if len("".join(search_dict.values())) == 0:
        interface.show_script_error(
            STANDARD_SEARCH_ERROR_MESSAGE
        )
        return

    # Set up and run the CSD search
    interface.update_progress("Running search for %s..." % str(search_dict))

    sql_query_string, search_report = yield_sql_query(
        search_dict, query_dict=query_dict
    )
    if isinstance(sqlite_path, (Path, str)):
        sql_connection = sqlite3.connect(sqlite_path)
        additional_databases = []
    elif isinstance(sqlite_path, (list, tuple)):
        sql_connection = sqlite3.connect(sqlite_path[0])
        additional_databases = sqlite_path[1:]

    if "first_row_name" in settings.keys():
        first_row_name = settings["first_row_name"]
    elif "NAME" in settings.get("output_columns"):
        first_row_name = "NAME"
    elif "refcode" in settings.get("output_columns"):
        first_row_name = "refcode"
    elif "ident" in settings.get("output_columns"):
        first_row_name = "ident"
    else:
        table_name = list(query_dict.values())[0].get('table')
        row_name_gen = sql_connection.execute(f"PRAGMA table_info({table_name});")
        print(next(row_name_gen))
        first_row_name = next(row_name_gen)[1]

    length_query = sql_query_string.replace("SELECT * ", f"SELECT COUNT({first_row_name}) ")
    unique_structures_query = sql_query_string.replace("SELECT * ", f"SELECT COUNT(DISTINCT({first_row_name})) ")

    maxhits = search_dict.get("max_hits_override", settings.get("max_hits", None))
    if isinstance(maxhits, int):
        sql_query_string += f" LIMIT {maxhits}"

    result_df = pd.read_sql_query(sql_query_string, sql_connection)
    unique_structures = next(sql_connection.execute(unique_structures_query))[0]
    total_length = next(sql_connection.execute(length_query))[0]

    sql_connection.close()

    for sql_path in additional_databases:
        sql_connection = sqlite3.connect(sql_path)

        total_length += next(sql_connection.execute(length_query))[0]
        unique_structures += next(sql_connection.execute(unique_structures_query))[0]
        new_rows = pd.read_sql_query(sql_query_string, sql_connection)
        sql_connection.close()

        if len(new_rows) > 0:
            result_df = pd.concat([result_df, new_rows])

    result_df["NAME"] = result_df[first_row_name]
    result_df["Identifier"] = [
        f"{interface.identifier}|{c}|{b}" for b, c in enumerate(result_df[first_row_name])
    ]

    # Write hits to the output GCD and TSV files
    interface.update_progress("Writing results to output files...")

    refcodes = result_df[first_row_name].unique()

    with ccdc.utilities.output_file(interface.output_gcd_file) as gcd_file:
        gcd_file.write("\n".join([x.partition("_")[0] for x in refcodes]))

    write_tsv_from_dataframe(
        result_df[settings.get("output_columns", result_df.columns)],
        interface.output_tsv_file)

    write_search_result_html(
        title=settings.get("title", "Search"),
        data=f"""
        Query: {search_report} ({sql_query_string}) <br>
        Result:{total_length} hits in {unique_structures} structures (multiple databases) <br>
        Display:{len(result_df)} hits in {len(refcodes)} unique structures <br>
        {settings.get("caption", "")}
        """,
        path=interface.output_html_file)


class SearchDialog2(ttk.Frame):
    def __init__(self, query_dict: dict):
        super().__init__()
        self.query_dict = query_dict
        self.entries = list()
        self.outputs = [""] * len(self.query_dict)
        self.custom_font = tkinter.font.Font(family="Arial", size=14)
        self.initUI()

    def generate_entry(self, output, parent_frame=None, input_dict: dict = {}, index=0):
        text = input_dict.get("label", "")
        hovertext = input_dict.get("hovertext", "")

        def on_enter(enter):
            self.infobox.configure(text=hovertext)

        if parent_frame:
            frame = ttk.Frame(parent_frame)
        else:
            frame = ttk.Frame(self)

        label = ttk.Label(frame, text=text, width=20)
        entry = ttk.Entry(frame, textvariable=output)

        label.pack(side=tkinter.LEFT, expand=True, fill=tkinter.BOTH, padx=1, pady=5)
        entry.pack(side=tkinter.LEFT, expand=True, fill=tkinter.BOTH, padx=1, pady=5)

        frame.bind("<Enter>", on_enter)
        frame.pack(side=tkinter.TOP)
        self.entries.append(entry)

    def initUI(self):

        self.master.title("Search for topology")

        rframe = ttk.Frame(self, width=400, height=1000, border=1)

        self.infobox = tkinter.Label(rframe,
                                     text="Hover for description",
                                     width=400,
                                     height=150,
                                     background="ivory1",
                                     border=1,
                                     wraplength=500)

        self.imagebox = tkinter.Label(rframe,
                                      text="image here",
                                      width=400,
                                      height=150,
                                      background="lightblue", border=1)

        self.infobox.pack(side=tkinter.TOP, expand=True)
        self.imagebox.pack(side=tkinter.TOP, expand=True)

        lframe = ttk.Frame(self, width=300, height=1000, border=1)
        [
            self.generate_entry(output=output, parent_frame=lframe, input_dict=input_dict, index=index)
            for index, (input_dict, output) in
            enumerate(zip(self.query_dict.values(), self.outputs))
        ]

        btn = tkinter.Button(lframe, text="Submit", command=self.onSubmit)
        btn.pack(side=tkinter.TOP, padx=1, pady=5)

        lframe.pack(side=tkinter.LEFT, fill=tkinter.Y)
        rframe.pack(side=tkinter.LEFT, fill=tkinter.Y)
        self.pack()

    def onSubmit(self):
        self.outputs = [x.get() for x in self.entries]
        self.quit()


class SearchDialog(ttk.Frame):

    def __init__(self, query_dict):
        super().__init__()
        self.query_dict = query_dict
        self.entries = list()
        self.outputs = [""] * len(self.query_dict)
        self.initUI()

    def generate_entry(self, text, output):
        frame = ttk.Frame(self)
        frame.pack(fill=tkinter.X)
        label = ttk.Label(frame, text=text, width=30)
        label.pack(side=tkinter.LEFT, padx=5, pady=9)
        entry = ttk.Entry(frame, textvariable=output)
        entry.pack(fill=tkinter.X, padx=5, expand=True)
        self.entries.append(entry)

    def initUI(self):

        self.master.title("Simple Dialog")
        self.pack(fill=tkinter.BOTH, expand=True)

        [
            self.generate_entry(text, output)
            for text, output in zip(
                [x["label"] for x in self.query_dict.values()], self.outputs
            )
        ]

        frame_button = ttk.Frame(self)
        frame_button.pack(fill=tkinter.X)

        # Command tells the form what to do when the button is clicked
        btn = ttk.Button(frame_button, text="Submit", command=self.onSubmit)
        btn.pack(padx=5, pady=10)

    def onSubmit(self):
        self.outputs = [x.get() for x in self.entries]
        self.quit()


def get_query_text_with_hovertext(query_dict):
    """Open a GUI dialog to ask for a search term.

    :returns: (:obj:`str`) The search query entered in the dialog.
    """
    # Create a Tk GUI instance and sets the size
    height = len(query_dict.keys()) * 40
    root = tkinter.Tk()
    root.geometry(f"1000x{height}+300+300")
    app = SearchDialog2(query_dict)

    root.mainloop()

    # Vars are where it retrieves the app outputs (the values you entered) into the query
    vars = {k: v for k, v in zip(query_dict.keys(), app.outputs)}

    try:
        root.destroy()
    except tkinter.TclError:
        pass

    return vars


def get_query_text(query_dict):
    """Open a GUI dialog to ask for a search term.

    :returns: (:obj:`str`) The search query entered in the dialog.
    """
    # Create a Tk GUI instance and sets the size
    height = len(query_dict.keys()) * 50
    root = tkinter.Tk()
    root.geometry(f"400x{height}+300+300")
    app = SearchDialog(query_dict)

    root.mainloop()

    # Vars are where it retrieves the app outputs (the values you entered) into the query
    vars = {k: v for k, v in zip(query_dict.keys(), app.outputs)}

    try:
        root.destroy()
    except tkinter.TclError:
        pass

    return vars


def query_string_to_tuple(query_string, default_upper_limit=9999.9):
    """this converts a string to a tuple of values representing the upper and lower bounds.
    Compatible with a single value (x, <x, >x) or a range (x-y)"""
    if "-" in query_string:
        ll, _, ul = query_string.partition("-")
    elif "<" in query_string:
        ll, _, ul = query_string.partition("<")
        ll = ll.rstrip("=")
        ul = ul.lstrip("=")
    elif ">" in query_string:
        ul, _, ll = query_string.partition(">")
        ll = ll.lstrip("=")
        ul = ul.rstrip("=")
    elif "=" in query_string:
        _, _, ll = query_string.partition("=")
        ul = ll
    elif query_string in ("0", "0.0"):
        ll, ul = 0.0, 0.1
    else:
        ll, ul = float(query_string) * 0.95, float(query_string) * 1.05

    if ll != "":
        lower_limit = float(ll)
    else:
        lower_limit = 0.0
    if ul != "":
        upper_limit = float(ul)
    else:
        upper_limit = default_upper_limit
    return (lower_limit, upper_limit)


ATOMIC_NUMBERS = {
    at: ix + 1
    for ix, at in enumerate(
        ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al",
         "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
         "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
         "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
         "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
         "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
         "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
         "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
         "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
         "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]
    )
}


def atomic_symbols_between(lano, uano):
    # lower atomic number, upper
    return [
        k
        for k, v in ATOMIC_NUMBERS.items()
        if all((((v >= lano) or (lano < 0)),
                ((v < uano) or (uano < 0))))
    ]


def query_atomic_symbols_split(query_string):
    """this converts a string representing multiple atomic symbols to a list"""

    if "-" in query_string:
        lat, _, uat = query_string.partition("-")
        lano = ATOMIC_NUMBERS.get(lat.strip(), -2)
        uano = ATOMIC_NUMBERS.get(uat.strip(), -2) + 1
        atomic_symbols = atomic_symbols_between(lano, uano)

    elif "<=" in query_string:
        lat, _, uat = query_string.partition("<=")
        lano = ATOMIC_NUMBERS.get(lat.strip(), -2)
        uano = ATOMIC_NUMBERS.get(uat.strip(), -2)
        atomic_symbols = atomic_symbols_between(lano, uano)

    elif "<" in query_string:
        lat, _, uat = query_string.partition("<")
        lano = ATOMIC_NUMBERS.get(lat.strip(), -2) + 1
        uano = ATOMIC_NUMBERS.get(uat.strip(), -2) - 1
        atomic_symbols = atomic_symbols_between(lano, uano)

    elif ">=" in query_string:
        uat, _, lat = query_string.partition(">=")
        lano = ATOMIC_NUMBERS.get(lat.strip(), -2)
        uano = ATOMIC_NUMBERS.get(uat.strip(), -2)
        atomic_symbols = atomic_symbols_between(lano, uano)

    elif ">" in query_string:
        uat, _, lat = query_string.partition(">")
        lano = ATOMIC_NUMBERS.get(lat.strip(), -2) + 1
        uano = ATOMIC_NUMBERS.get(uat.strip(), -2) - 1
        atomic_symbols = atomic_symbols_between(lano, uano)

    elif "," in query_string:
        atomic_symbols = [
            x.strip()
            for x in query_string.split(",")
            if (x.strip() in ATOMIC_NUMBERS.keys())
        ]

    elif " " in query_string:
        atomic_symbols = [
            x.strip()
            for x in query_string.split()
            if (x.strip() in ATOMIC_NUMBERS.keys())
        ]

    elif query_string.strip() in ATOMIC_NUMBERS.keys():
        atomic_symbols = [query_string.strip()]
    else:
        atomic_symbols = list()
    return ",".join([f"'{x}'" for x in atomic_symbols])


def yield_sql_query(search_dict, query_dict):
    query_report = str()
    subqueries = list()

    for key, value in search_dict.items():
        if not value:
            continue

        key_info = query_dict[key]

        if key_info["process"] == "query_atomic_symbols_split":
            val = query_atomic_symbols_split(value)
            query_report += key_info["report_query"].format(val=val)
            subqueries.append(key_info["subquery"].format(val=val))

        elif key_info["process"] == "query_string_to_tuple":
            vals = query_string_to_tuple(value)
            query_report += key_info["report_query"].format(val1=vals[0], val2=vals[1])
            subqueries.append(key_info["subquery"].format(val1=vals[0], val2=vals[1]))

        elif key_info["process"] == "query_string_lookup":
            val = key_info["lookup_dict"].get(value.lower(), None)
            query_report += key_info["report_query"].format(val=val)
            subqueries.append(key_info["subquery"].format(val=val))

        elif key_info["process"] == "query_string_refcode_splitter":
            import re

            vals = re.findall(
                "[A-Z]{6}[0-9]{0,2}",
                value.upper(),
            )
            vals_str = "(" + ",".join(["'" + x + "'" for x in vals if x]) + ")"
            query_report += key_info["report_query"].format(val=vals_str)
            subqueries.append(key_info["subquery"].format(val=vals_str))

        elif key_info["process"] == "query_int_as_substring":
            import re
            import numpy as np
            vals = re.findall("[0-9][0-9]|[0-9]", value)
            vals_as_uq = dict(zip(*np.unique([str(int(x)) for x in vals], return_counts=True)))

            vals_in_perc_signs = ["%" + ",".join([f"({x})"] * y) + "%"
                                  for x, y in vals_as_uq.items()]

            vals_str = "('" + "'&'".join(vals_in_perc_signs) + "')"
            subq = [
                key_info["subquery"].format(val=f"'{x}'")
                for x in vals_in_perc_signs
            ]
            query_report += key_info["report_query"].format(val=vals_str)
            [subqueries.append(x) for x in subq]

        elif key_info["process"] == "none":
            query_report += key_info["report_query"].format(val=value)
            subqueries.append(key_info["subquery"].format(val=value))

        elif key_info["process"] == "add_quotes":
            query_report += key_info["report_query"].format(val=f"'{value}'")
            subqueries.append(key_info["subquery"].format(val=f"'{value}'"))

    query = f"SELECT * FROM {key_info['table']} WHERE {' AND '.join(subqueries)}"

    return query, query_report

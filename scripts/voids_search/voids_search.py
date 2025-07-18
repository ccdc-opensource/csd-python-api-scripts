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

"""voids_search.py - return results from the calculated void properties of the
CSD for use in Mercury.
"""
import tkinter as Tkinter
import unicodedata

# Tell the CSD Python API not to create a QApplication object on macOS.
# This script doesn't require any QApplication-dependent functionality and on
# macOS, an Apple bug will cause Python crashes when using Tk if a QApplication
# object exists in the same Python script.
import sys
import os
if sys.platform.startswith("darwin"):
    os.environ['CCDC_PYTHON_API_NO_QAPPLICATION'] = '1'

from ccdc.search import TextNumericSearch
import ccdc.utilities
from jinja2 import Template
from pathlib import Path

from tkinter import BOTH, X, LEFT
from tkinter.ttk import Frame, Label, Entry, Button

pore_query_dict = {"total_surface_area": "Pore Total Surface Area (Å^2)",
                   "total_geometric_volume": "Pore Total Geometric Volume (Å^3)",
                   "pore_limiting_diameter": "Pore Limiting Diameter (Å)",
                   "max_pore_diameter": "Pore Maximum Diameter (Å)",
                   "num_percolated_dimensions": "Number of Percolated Dimensions"}


class VoidDialog(Frame):

    def __init__(self):
        super().__init__()
        self.entries = list()
        self.outputs = [""]*5
        self.initUI()

    def generate_entry(self, text, output):
        frame = Frame(self)
        frame.pack(fill=X)
        label = Label(frame, text=text, width=30)
        label.pack(side=LEFT, padx=5, pady=9)
        entry = Entry(frame, textvariable=output)
        entry.pack(fill=X, padx=5, expand=True)
        self.entries.append(entry)

    def initUI(self):

        self.master.title("Simple Dialog")
        self.pack(fill=BOTH, expand=True)

        [self.generate_entry(text, output) for text, output in zip(pore_query_dict.values(), self.outputs)]

        frame_button = Frame(self)
        frame_button.pack(fill=X)

        # Command tells the form what to do when the button is clicked
        btn = Button(frame_button, text="Submit", command=self.onSubmit)
        btn.pack(padx=5, pady=10)
        # btn.bind('<Return>', command=self.onSubmit)

    def onSubmit(self):
        self.outputs = [x.get() for x in self.entries]
        self.quit()


def get_query_text():
    """Open a GUI dialog to ask for a search term.

    :returns: (:obj:`str`) The search query entered in the dialog.
    """
    # Create a Tk GUI instance and sets the size
    root = Tkinter.Tk()
    root.geometry("400x300+300+300")
    app = VoidDialog()

    root.mainloop()

    # Vars are where it retrieves the app outputs (the values you entered) into the query
    values = {k: v for k, v in zip(pore_query_dict.keys(), app.outputs)}

    if root.winfo_exists():
        root.destroy()

    return values


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
    elif query_string in ("0", "0.0"):
        ll, ul = 0.0, 0.1
    else:
        ll, ul = float(query_string)*0.95, float(query_string)*1.05

    if ll != "":
        lower_limit = float(ll)
    else:
        lower_limit = 0.0
    if ul != "":
        upper_limit = float(ul)
    else:
        upper_limit = default_upper_limit
    return (lower_limit, upper_limit)


def convert_to_ascii(text):
    return unicodedata.normalize('NFKD', text).encode('ascii', 'ignore').decode('ascii')


def run_search(interface=None):
    """Search the CSD for a set of void descriptors. This will generate a list of structures

    :param interface: (:obj:`ccdc.utilities.ApplicationInterface`) An ApplicationInterface instance.
    """
    if not interface:
        # Create the ApplicationInterface without immediately parsing command line parameters
        # so that we can add a custom parameter to pass the search query in directly.
        interface = ccdc.utilities.ApplicationInterface(parse_commandline=False)

        # Add a custom command-line parameter for the compound name.
        interface.commandline_parser.add_argument(
            '--total_surface_area', required=False, default='',
            help='total_surface_area to search for, bypassing the GUI dialog.')
        interface.commandline_parser.add_argument(
            '--total_geometric_volume', required=False, default='',
            help='total_geometric_volume to search for, bypassing the GUI dialog.')
        interface.commandline_parser.add_argument(
            '--pore_limiting_diameter', required=False, default='',
            help='pore_limiting_diameter to search for, bypassing the GUI dialog.')
        interface.commandline_parser.add_argument(
            '--max_pore_diameter', required=False, default='',
            help='max_pore_diameter to search for, bypassing the GUI dialog.')
        interface.commandline_parser.add_argument(
            '--num_percolated_dimensions', required=False, default='',
            help='num_percolated_dimensions to search for, bypassing the GUI dialog.')

        # Parse the commandline including checking for the compound parameter.
        interface.parse_commandline()

    # open a GUI dialog asking for the search query.
    void_search_dict = get_query_text()

    query_report = str()
    # Error out if the compound name given is empty
    if len("".join(void_search_dict.values())) == 0:
        interface.show_script_error('''Nothing was given to search for. Please enter either one number
    (to search +/- 5 %) or two numbers separated by a "-" in one of the available boxes''')
        return

    # Set up and run the CSD search
    interface.update_progress('Running search for pore %s ...' % str(void_search_dict))

    TNS = TextNumericSearch()
    if void_search_dict.get("total_surface_area", False):
        tsa = query_string_to_tuple(void_search_dict["total_surface_area"], 10000.0)
        query_report += f"Pore Total Surface Area (&#8491;^2) {tsa[0]}-{tsa[1]} <br>"
        TNS.add_pore_analysis_total_surface_area(tsa)

    if void_search_dict.get("total_geometric_volume", False):
        tgv = query_string_to_tuple(void_search_dict["total_geometric_volume"], 100000.0)
        query_report += f"Pore Total Geometric Volume (&#8491;^3): {tgv[0]}-{tgv[1]} <br>"
        TNS.add_pore_analysis_total_geometric_volume(tgv)

    if void_search_dict.get("pore_limiting_diameter", False):
        pld = query_string_to_tuple(void_search_dict["pore_limiting_diameter"], 200.0)
        query_report += f"Pore Limiting Diameter (&#8491;): {pld[0]}-{pld[1]} <br>"
        TNS.add_pore_analysis_pore_limiting_diameter(pld)

    if void_search_dict.get("max_pore_diameter", False):
        pmd = query_string_to_tuple(void_search_dict["max_pore_diameter"], 200.0)
        query_report += f"Pore Maximum Diameter (&#8491;): {pmd[0]}-{pmd[1]} <br>"
        TNS.add_pore_analysis_max_pore_diameter(pmd)

    if void_search_dict.get("num_percolated_dimensions", False):
        npd = query_string_to_tuple(void_search_dict["num_percolated_dimensions"], 3.1)
        query_report += f"Number of Percolated Dimensions: {npd[0]}-{npd[1]} <br>"
        TNS.add_pore_analysis_max_pore_diameter(npd)

    hits = TNS.search()
    refcode_list = list()

    # Write hits to the output GCD and TSV files
    interface.update_progress('Writing results to output files...')
    with ccdc.utilities.output_file(interface.output_gcd_file) as gcd_file, \
            ccdc.utilities.CSVWriter(interface.output_tsv_file,
                                     header=["Identifier",
                                             "Pore Total Surface Area (\u212b^2)",
                                             "Pore Total Geometric Volume (\u212b^3)",
                                             "Pore Limiting Diameter (\u212b)",
                                             "Pore Maximum Diameter (\u212b)",
                                             "Pore Number of Percolated Dimensions",
                                             'Chemical Name(s)',
                                             ],
                                     delimiter='\t') as tsv_file:

        for h in hits:
            names = h.entry.chemical_name if h.entry.chemical_name else ''
            if h.entry.synonyms:
                names += ';' + ';'.join(h.entry.synonyms)
            names = convert_to_ascii(names)

            # Write identifier to GCD file and identifier and names to TSV
            print(h.identifier, file=gcd_file)
            pa = h.entry.calculated_properties.pore_analyser
            refcode_list.append(h.identifier)
            tsv_file.write_row([h.identifier,
                                pa.total_surface_area,
                                pa.total_geometric_volume,
                                pa.pore_limiting_diameter,
                                pa.max_pore_diameter,
                                pa.num_percolated_dimensions,
                                names])

    with open(interface.output_html_file, "w") as report:
        tl = Template(
            open(
                Path(__file__).parent / "search_result_template.html",
                "r",
            ).read()
        )
        report.write(
            tl.render(
                title="Voids_search",
                data=f"""
                Query: {query_report} <br>
                Result:{len(refcode_list)} hits in {len(set(refcode_list))} structures <br>
                More information on the pore calculations is available here:
                <a href="https://downloads.ccdc.cam.ac.uk/documentation/API/descriptive_docs/calculated_properties.html">
                link
                </a>
                """,
            )
        )


if __name__ == '__main__':
    run_search()

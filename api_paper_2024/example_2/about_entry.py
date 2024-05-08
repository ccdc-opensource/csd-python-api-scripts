# -*- coding: utf-8 -*-
#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
"""
Write out a simple HTML report that contains the CrossRef (see https://crossref.org/ ) of the publication associated with a CSD entry and
mines OpenAlex (see https://openalex.org ) for references and citations.

To run this you will need to install habanero - see https://pypi.org/project/habanero/ along side the CCDC python API
"""

# For now, to prevent certain packages tripping things up when running
# through mercury.
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)
import os
import time

from ccdc.utilities import html_table, Logger
from ccdc.search import TextNumericSearch, CombinedSearch
from ccdc.io import EntryReader
from ccdc.utilities import ApplicationInterface
from references import ReferencesAndCitations, ReferencesAndCitationsError

logger = Logger()

#
# Imports - you will need to install habanero. Wordcloud and nltk are optional

have_habenero = True
have_wordcloud = True
try:
    import crossref
except ImportError as e:
    have_habanero = False
    print(e)

# This is optional, but if it can be imported, the report will have a word cloud created

try:
    import language_processing
except ImportError as e:
    have_wordcloud = False
    print(e)


def wordcloud(list_of_strings, output_directory_path, identifier):
    """ Use the language processing code to create a wordcloud

    :param list_of_string: a list of strings to add to the input text (e.g. could be a list of titles)
    :param output_directory_path: path to a directory to write the word cloud to
    :param identifier: some unique identifier to include in the output file name
    """
    if have_wordcloud:
        text = '\n'.join([x for x in list_of_strings if x is not None])
        if len(text) > 0:
            fname = os.path.join(output_directory_path, f"WordCloud_{identifier}.png")
            language_processing.word_frequency_analysis(1, text, fname)
            return fname
    return None


def get_entries_for_doi(doi):
    """
    Search the Cambridge Structural Database for any entries associated with a single DOI
    :param doi: The document object identifier to search for
    """
    start = time.time()
    if doi is None:
        return []

    searcher = TextNumericSearch()
    searcher.add_doi(doi)
    x = searcher.search()
    logger.debug(f"Search time for doi {doi}: {time.time() - start}")
    return x


def get_entries_for_dois(dois):
    """
    Search the Cambridge Structural Database for any entries associated with a list of DOIs
    :param dois: The document object identifiers to search for
    """
    start = time.time()
    if len(dois) == 0:
        return []

    if len(dois) == 1:
        return get_entries_for_doi(dois[0])

    # As we are searching for more than one DOI lets build a combination search 'ORing' each
    # case
    combo = TextNumericSearch()
    combo.add_doi(dois[0])

    for d in dois[1:]:
        ts = TextNumericSearch()
        ts.add_doi(d)
        combo |= ts

    # We can then run the search using the combined search method
    searcher = CombinedSearch(combo)
    x = searcher.search()
    logger.debug(f"Search time for multiple dois : {time.time() - start}")
    return x


def linked_csd_entries(doi, root_identifier):
    hits = get_entries_for_doi(doi)
    return [hit.identifier for hit in hits if hit.identifier != root_identifier]


def tabulate_crossref_record(publication_doi, logger, email_address=None):
    """
    Search in CrossRef tabulate information returned as HTML
    """
    crossref_info = crossref.CrossrefSearch(logger, email_address).crossref_record(publication_doi)
    data = [[f"<b>{k}</b>", f"{v}"] for k, v in crossref_info.items()]
    data.append(["<b>Link to publication</b>", f"<a href=https:/www.doi.org/{publication_doi}>Click here</a>"])
    return data


def search_citation_data(doi, root_identifier, email_address=None, datasource="openalex"):
    """
    Search in either OpenAlex or OpenCitations for citations and references
    """

    if doi is None:
        logger.debug("No doi to search")
        return []

    # Run the search in the external data sources for this document object identifier

    metadata_searcher = ReferencesAndCitations(logger, email_address, datasource)
    metadata_searcher.run(doi)

    # Go and search the CSD for citations and references with the associated DOIs

    hits = get_entries_for_dois(metadata_searcher.citation_dois)
    citing_identifiers = sorted([hit.identifier for hit in hits if hit.identifier != root_identifier])

    hits = get_entries_for_dois(metadata_searcher.reference_dois)
    referenced_identifiers = sorted([hit.identifier for hit in hits if hit.identifier != root_identifier])

    return [citing_identifiers, referenced_identifiers,
            metadata_searcher.reference_titles + metadata_searcher.citation_titles]


def write_entry_list(out_file, identifiers):
    """
    For a given set of CSD identifiers (refcodes) look up data in the CSD and write
    a tab-seperated file that can be sent to Mercury

    :param out_file: The output file to write to
    :param identifiers: A list of identifiers
    """
    with EntryReader('CSD') as er, open(out_file, 'w', encoding="utf-8") as out_tsv:
        out_tsv.write('refcode\tcompound name\tsynomnyms\n')
        for id in sorted(identifiers):
            entry = er.entry(id)
            out_tsv.write(f'{entry.identifier}\t"{entry.chemical_name}"\t"{" ".join(entry.synonyms)}"\n')


def write_report(interface: ApplicationInterface, doi: str, email_address: str, identifier: str, source: str):
    with interface.html_report(
            title='Citation and Reference Information for Paper DOI associated with %s' % doi) as report:

        # Write the section header for Entry Details
        report.write_section_header('Crossref Record')
        # Assemble a dictionary of information and labels to go into the "Entry Details" table
        # Generate a HTML table from the entry details and write it to the report

        data = tabulate_crossref_record(doi, logger, email_address=email_address)

        identifiers = linked_csd_entries(doi, identifier)
        if len(identifiers) > 0:
            data.append(["<b>Other entries from the same paper</b>", '\n'.join(identifiers)])
        else:
            data.append(["<b>No other entries from the same paper</b>"])
        try:
            citation_data = search_citation_data(doi, identifier, email_address=email_address, datasource=source)

        except ReferencesAndCitationsError as e:
            citation_data = None
            data.append(["<b>Citation lookup failed</b>", str(e)])

        if citation_data:
            data.append(["<b>Number of cited references found</b>", str(len(citation_data[1]))])
            data.append(["<b>Number of papers that cite this paper</b>", str(len(citation_data[0]))])

        report.write(html_table(data=data, table_id='crossref_details'))
        if citation_data:
            file_name = wordcloud(citation_data[2], interface.output_directory_path, identifier)
            all_ids = citation_data[0] + citation_data[1]
            if len(all_ids) > 0:
                write_entry_list(interface.output_tsv_file, all_ids)
            else:
                logger.info("No associated CSD entries found in citing articles or referenced articles")

            if file_name and os.path.exists(file_name):
                report.write_section_header('Data Analysis')
                report.write_figure(file_name, caption="Wordcloud of related article titles")


############################################################################################################
def main(interface=None):
    if interface is None:
        # Create the ApplicationInterface without immediately parsing command line parameters
        # so that we can add a custom parameter to pass the search query in directly.
        interface = ApplicationInterface(parse_commandline=False)

        # Add a custom command-line parameter for the doi or refcode.
        interface.commandline_parser.add_argument(
            '-d', '--doi-or-refcode', required=False, default=None,
            help='Refcode or doi to search for, allowing command line use')

        # Add a custom command-line parameter for the doi or refcode.
        interface.commandline_parser.add_argument(
            '-e', '--email-address', required=False, default=None,
            help='Email to use in APIs (to be polite)')

        # Add a custom command-line parameter for the doi or refcode.
        interface.commandline_parser.add_argument(
            '-s', '--source', required=False, default="openalex",
            choices=["openalex", "opencitations", "all"],
            help='Where to search for reference information')

        choice_map = {"debug": Logger.DEBUG, "info": Logger.INFO, "warning": Logger.WARNING, "error": Logger.ERROR,
                      "critical": Logger.CRITICAL}
        interface.commandline_parser.add_argument('-l', '--log-level',
                                                  default="error",
                                                  choices=choice_map.keys(),
                                                  help="Log level"
                                                  )

        # Parse the commandline including checking for the compound parameter.
        interface.parse_commandline()

        logger.set_log_level(choice_map[interface.log_level])

        entry_id = interface.doi_or_refcode
        source = interface.source
        email_address = interface.email_address
    else:
        logger.set_log_level(Logger.ERROR)
        source = "openalex"
        entry_id = None
        email_address = None

    if have_habenero is False:
        interface.exit_with_error("Your python distribution doesnt have habanero installed.\n\n"
                                  "See https://habanero.readthedocs.io/en/latest/intro/install.html "
                                  "for more information on installing this package")
    if entry_id is None:
        entry = interface.current_entry
        doi = entry.publication.doi
        identifier = entry.identifier
    else:
        reader = EntryReader('CSD')
        try:
            entry = reader.entry(entry_id)
            doi = entry.publication.doi
            identifier = entry.identifier
        except RuntimeError:
            # No entry with this name: treat it as a raw doi
            doi = entry_id
            identifier = entry_id.replace('/', '_').replace('(', '_').replace(')', '_')
    write_report(interface, doi, identifier, email_address, source)


if __name__ == "__main__":
    main()

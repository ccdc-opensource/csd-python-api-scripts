#!/usr/bin/env python
#
# This script can be used for any purpose without limitation subject to the
# conditions at https://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2024-05-02: created by the Cambridge Crystallographic Data Centre
#

"""
Written by Jason C. Cole

Credit: to Unai Saralegui
Much of this code was inspired by code in the opencitingpy project.

See https://pypi.org/project/opencitingpy/
"""

from http import HTTPStatus
# When we move to Python 3.12 we can import batched from itertools directly
#
# from itertools import batched
from itertools import islice

import url_requesting


def batched(iterable, n):
    # batched('ABCDEFG', 3) --> ABC DEF G
    if n < 1:
        raise ValueError('n must be at least one')
    it = iter(iterable)
    while batch := tuple(islice(it, n)):
        yield batch


OPEN_CITATIONS_API_URL = 'https://w3id.org/oc/index/api/v1'


class OpenCitationsError(RuntimeError):
    def __init__(self, msg):
        super().__init__(msg)


class OpenCitationsMetadataResult:
    """
    Credit:
    This is adapted from opencitingpy (the Metadata class) by Unai Saralegui
    for parsing the opencitation meta data object returned as json
    """

    def __init__(self, data):
        """
        Custom class for OpenCitation Metadata type objects.

        :param data: The dictionary (json) with the information for each item returned by the OpenCitations API
        :type data: dict
        """

        def make_list(input_string):
            if input_string is None:
                return None
            return [item.strip() for item in input_string.split("; ")]

        self.author = make_list(data.get('author', None))
        self.year = data.get('year', None)
        self.title = data.get('title', None)
        self.source_title = data.get('source_title', None)
        self.source_id = data.get('source_id', None)
        self.volume = data.get('volume', None)
        self.issue = data.get('issue', None)
        self.page = data.get('page', None)
        self.doi = data.get('doi', None)
        self.reference = make_list(data.get('reference', None))
        self.citation = make_list(data.get('citation', None))
        self.citation_count = data.get('citation_count', None)
        self.oa_link = data.get('oa_link', None)
        self.data = data

    def __str__(self):
        d = dict(self.data)
        d['author'] = self.author
        d['reference'] = self.reference
        d['citation'] = self.citation
        return str(d)


class OpenCitationsSearcher:

    def __init__(self, logger=None, email=None):
        self.logger = logger
        self.email = email

    def run(self, dois):

        def run_batched(subset, size):
            # Run max of 10 dois at a time to avoid bad request errors
            batches = list(batched(subset, size))
            failed = []
            for batch in batches:
                try:
                    failed += self.run_subset(batch)
                except url_requesting.URLRequestError as exc:
                    code = exc.status_code
                    if code == HTTPStatus.BAD_REQUEST or code == HTTPStatus.NOT_FOUND:
                        failed += batch
                    else:
                        raise
            return failed

        if isinstance(dois, str):
            dois = [dois]

        self._raw_results = []

        # If there is a bad DOI one gets a bad request response
        # so divide into ever smaller batches
        failed = dois
        size = len(dois)
        if size > 100:
            size = 100
        while len(failed) > 0:
            failed = run_batched(failed, size)
            if size == 1:
                break
            new_size = 1
            if size > 2:
                new_size = min(int(len(failed) / 2), 1 + int(size / 2))
            size = new_size

        # Now parse the results
        self.results = [OpenCitationsMetadataResult(data) for data in self._raw_results]
        self.not_found = failed

    def run_subset(self, dois):

        doi_string = '__'.join(dois)
        url = f'{OPEN_CITATIONS_API_URL}/metadata/{doi_string}'

        res = url_requesting.URLRequest(self.logger).run(url)

        failed = []
        if len(res) != len(dois):
            hits = set([v['doi'] for v in res])
            failed = [d for d in dois if d not in hits]

        self._raw_results += res
        return failed


class WorkCitationsAndReferences:
    def __init__(self, logger=None, email=None):
        """
        Get doi for all citations and references to a root doi
        """
        self._searcher = OpenCitationsSearcher(logger, email)
        self.citation_info = []
        self.reference_info = []

    def run(self, doi):
        self._searcher.run(doi)
        metadata = self._searcher.results
        # harvest any citing dois and reference dois
        if len(metadata) > 0:
            citation_dois = metadata[0].citation
            reference_dois = metadata[0].reference
            self._searcher.run(citation_dois)
            self.citation_info = [(d.doi, d.title) for d in self._searcher.results]
            self._searcher.run(reference_dois)
            self.reference_info = [(d.doi, d.title) for d in self._searcher.results]


def test_dois():
    test_dois_links = ['10.1021/jp103212z', '10.1016/j.jorganchem.2012.03.020',
                       '10.1016/j.poly.2018.06.030', '10.1002/bip.360231115',
                       '10.1007/978-1-4684-1354-0_4', '10.1002/chin.197846077',
                       '10.1038/292055a0']

    import ccdc.utilities

    logger = ccdc.utilities.Logger()
    logger.set_log_level(ccdc.utilities.Logger.DEBUG)
    x = OpenCitationsSearcher(logger)
    x.run(test_dois_links)
    print(x.results)
    print(x.not_found)

    test_dois_links.append('10.9999/test_doi')
    x.run(test_dois_links)
    print(x.results)
    print(x.not_found)

# -*- coding: utf-8 -*-
#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
"""
Wrapper classes to unfiy citation data sources
"""

import openalex
import opencitations
import url_requesting


class ReferencesAndCitationsError(RuntimeError):
    def __init__(self, msg):
        super().__init__(msg)


class ReferencesAndCitations:
    """
        Wrapper class to allow different data sources for harvesting citation and reference information
    """

    def __init__(self, logger=None, email_address=None, datasource="openalex"):
        self._logger = logger
        self._email_address = email_address
        self._datasource = datasource
        self._citation_info = []
        self._reference_info = []

    def _openalex(self, doi):
        try:
            metadata_searcher = openalex.WorkCitationsAndReferences(self._logger, self._email_address)
            metadata_searcher.run(doi)
        except url_requesting.URLRequestError as e:
            raise ReferencesAndCitationsError(
                f"Unable to access meta-data from openalex for related DOIs due to exception {e}")
        self._citation_info += metadata_searcher.citation_info
        self._reference_info += metadata_searcher.reference_info

    def _opencitations(self, doi):
        try:
            metadata_searcher = opencitations.WorkCitationsAndReferences(self._logger, self._email_address)
            metadata_searcher.run(doi)
        except url_requesting.URLRequestError as e:
            raise ReferencesAndCitationsError(
                f"Unable to access meta-data from opencitations for related DOIs due to exception {e}")
        self._citation_info += metadata_searcher.citation_info
        self._reference_info += metadata_searcher.reference_info

    def run(self, doi):
        if self._datasource == "openalex":
            self._openalex(doi)
        elif self._datasource == "opencitations":
            self._opencitations(doi)
        else:
            self._openalex(doi)
            self._opencitations(doi)

    @property
    def reference_titles(self):
        return list(set([t[1] for t in self._reference_info if t[1] is not None]))

    @property
    def citation_titles(self):
        return list(set([t[1] for t in self._citation_info if t[1] is not None]))

    @property
    def reference_dois(self):
        return list(set([t[0] for t in self._reference_info if t[0] is not None]))

    @property
    def citation_dois(self):
        return list(set([t[0] for t in self._citation_info if t[0] is not None]))

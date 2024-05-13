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

import time

from habanero import Crossref, WorksContainer, RequestError
from habanero import counts
from requests.exceptions import HTTPError, RequestException

TZERO = time.time()


class CrossrefSearch:
    def __init__(self, logger, email_address=None):
        """
        :param logger: a CCDC logger to log outcomes
        :param email_address: An email address. If not None and a valid address, CrossRef will use the polite API
        """
        self._email_address = email_address
        self._logger = logger
        self._repeat_pause = 0.03
        self._not_available_message = "Not Available in CrossRef"

    def debug_log(self, msg):
        if self._logger is not None:
            self._logger.debug(msg)

    def is_crossref_doi(self, publication_doi):
        """
         Look up if a given doi is a CrossRef doi

        :param publication_doi: the DOI
        :param email_address: An email address. If not None and a valid address, CrossRef will use the polite API
        """
        if publication_doi is not None and len(publication_doi) > 0:
            time.sleep(self._repeat_pause)  # Be nice to crossref - dont hit it too hard.
            cr = Crossref(mailto=self._email_address)
            try:
                ra = cr.registration_agency(publication_doi)
                return ra and len(ra) > 0
            except KeyError:
                pass
            except HTTPError:
                pass

        return False

    def crossref_record(self, publication_doi):
        """
        Generate a dictionary based on information in CrossRef for a publication DOI
        :param publication_doi: a publication object possibly extracted from the CSD `ccdc.entry.Entry`
        """

        self.debug_log(f"Start CROSSRef tabulation {time.time() - TZERO}")
        data = {}
        if self.is_crossref_doi(publication_doi) is False:
            data["Information"] = "Paper DOI not available"
            return data

        result = None
        try:
            cr = Crossref(mailto=self._email_address)

            cr_data = WorksContainer(cr.works(publication_doi))
            citation_count = counts.citation_count(publication_doi)
            self.debug_log(f"Time after looking up publication doi in CrossRef: {time.time() - TZERO}")
            titles = cr_data.title[0]
            # Appears to be a list of lists
            if len(titles) > 1:
                titles = "\n".join(titles)
            elif len(titles) == 0:
                titles = "Not available"
            else:
                titles = titles[0]
            data["Titles"] = titles
            subjects = self._not_available_message
            if cr_data.subject[0]:
                subjects = "; ".join(cr_data.subject[0])
            data["Scopus Subject Keywords"] = subjects

            try:
                first_author_and_institution = "None listed"
                if 'author' in cr_data.works[0]:
                    first_author_list = self.get_author_list(authors=cr_data.works[0]['author'])
                    first_author_and_institution = "; ".join(first_author_list)
            except AttributeError:
                first_author_and_institution = self._not_available_message

            data["First Author(s) & Affiliation(s)"] = first_author_and_institution

            data["Paper Citation Count"] = str(citation_count)

            try:
                abstract = cr_data.abstract[0].replace("jats:", "")
            except AttributeError:
                abstract = self._not_available_message
            data["Abstract"] = abstract

            try:
                funder = "None listed"
                if 'funder' in cr_data.works[0]:
                    funder = "; ".join([x['name'] for x in cr_data.works[0]['funder']])
            except AttributeError:
                funder = self._not_available_message

            data["Funder(s)"] = funder

        except RequestError as re:
            result = f"Look up failed with an exception {re}"
            if re.status_code == 404:
                result = "DOI not found in CrossRef"
        except RequestException as e:
            result = f"Look up failed with an exception {e}"

        if result is not None:
            data["Result"] = result
        return data

    @staticmethod
    def get_author_list(authors):
        first_author_list = []
        for author in authors:
            if author['sequence'] == 'first':
                v = f"{author['given']} {author['family']}"
                if 'affiliation' in author and len(author['affiliation']) > 0:
                    affiliations = ", ".join([a['name'] for a in author['affiliation']])
                    v += f" ({affiliations})"
                first_author_list.append(v)
        return first_author_list

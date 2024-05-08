"""
Credit: by Unai Saralegui

Much of this code was inspired by code in the opencitingpy project.
"""

import url_requesting


class OpenAlexSearcher:
    def __init__(self, logger=None, email=None):
        self.logger = logger
        self.email = email

    def run(self, url):
        return url_requesting.URLRequest(self.logger).run(url)

    def add_parameters(self, url, select_parameters):
        if self.email and url.find("email") == -1:
            url += f"?email={self.email}"
            if len(select_parameters) > 0:
                url += '&'

        if len(select_parameters) > 0:
            if url[-1] != '&':
                if url.find('?') == -1:
                    url += '?'
                else:
                    url += '&'
            url += f"select={','.join(select_parameters)}"
        return url

    def convert_to_api_url(self, raw_url):
        return raw_url.replace("https://openalex.org/", "https://api.openalex.org/works/")

    def works(self, doi, select_parameters):
        ref = f'https://api.openalex.org/works/https://doi.org/{doi}'
        self.add_parameters(ref, select_parameters)
        return self.run(ref)


class WorkCitationsAndReferences:

    def __init__(self, logger=None, email=None):
        """
        Get doi for all citations and references to a root doi
        """
        self._searcher = OpenAlexSearcher(logger, email)
        self.reference_info = []
        self.citation_info = []
        self._logger = logger

    def run(self, doi):

        raw_data = self._searcher.works(doi, ['cited_by_api_url', 'referenced_works'])
        ref = self._searcher.add_parameters(raw_data['cited_by_api_url'], ['doi', 'title'])
        citations = self._searcher.run(ref)

        self.citation_info = [(x['doi'].replace("https://doi.org/", "").encode("utf-8").decode('utf-8'),
                               x['title'].encode("utf-8").decode('utf-8')) for x in citations['results'] if
                              x['doi'] is not None]

        # Alas have to do magic with the references
        self.reference_info = []
        for ref in raw_data['referenced_works']:
            ref = self._searcher.convert_to_api_url(ref)
            ref = self._searcher.add_parameters(ref, ['title', 'doi'])
            ref_data = self._searcher.run(ref)

            doi = None
            if ref_data['doi'] is not None:
                doi = ref_data['doi'].replace("https://doi.org/", "")

            self.reference_info.append((doi, ref_data['title']))


def test_work_citation():
    x = WorkCitationsAndReferences(email=None)
    x.run('10.1021/jp103212z')
    print(x.citation_info)
    print(x.reference_info)

    try:
        x.run('10.8989/junk')
        print("Broken - should throw with an unknown doi")
    except url_requesting.URLRequestError as e:
        print(f"Threw as expected {e}")


if __name__ == "__main__":
    test_work_citation()

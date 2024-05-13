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
Credit - this code is adapted from  The impact of the Cambridge Structural Database and the small molecule crystal structures it contains: a bibliographic and literature study
Peter Willett, Jason C. Cole   and  Ian J. Bruno

@Article{D0CE00045K,
author ="Willett, Peter and Cole, Jason C. and Bruno, Ian J.",
title  ="The impact of the Cambridge Structural Database and the small molecule crystal structures it contains: a bibliographic and literature study",
journal  ="CrystEngComm",
year  ="2020",
volume  ="22",
issue  ="43",
pages  ="7233-7241",
publisher  ="The Royal Society of Chemistry",
doi  ="10.1039/D0CE00045K",
url  ="http://dx.doi.org/10.1039/D0CE00045K",
"""
# For now, to prevent deprecation warnings tripping things up
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)  # noqa

# Flake doesnt like imports after the line above, but we have to have the line above to
# prevent deprecation warnings from these third party packages sending spurious output
# which we dont want, as it is 'seen' as an error when using this code through Mercury.

from wordcloud import WordCloud # noqa
import nltk # noqa
import nltk.collocations # noqa
from nltk.stem.wordnet import WordNetLemmatizer # noqa
import pandas as pd # noqa
from csv import writer # noqa


class FrequencyCalculator:
    """
    A class to calculate the frequency of words in a text corpus that build on nltk
    """

    def __init__(self, text):

        def update_word(word):
            if word == "bonding":
                word = "bond"
            return word

        def keep(word, stop_words):
            if len(word) <= 3 and word.lower() != "tin":
                return False
            if word.startswith('/'):
                return False
            if word in stop_words:
                return False
            return True

        self.stop = set(nltk.corpus.stopwords.words('english'))
        self.bigrams_to_ignore = set()
        self.bigram_words_to_ignore = set()

        lem = WordNetLemmatizer()
        self._words = [update_word(lem.lemmatize(word)) for word in nltk.word_tokenize(text.lower()) if
                       keep(word, self.stop)]

        self._fdist = nltk.FreqDist(self._words)

    def most_common_words(self, how_many=None):
        n_elements = len(self._words)
        if how_many is None:
            rowdata = self._fdist.most_common()
        else:
            rowdata = self._fdist.most_common(how_many)
        return n_elements, rowdata

    def write_word_frequency_table(self, how_many, filename='single_word_frequencies.csv'):

        n_elements, rowdata = self.most_common_words(how_many)

        with open(filename, 'w', newline='', encoding='utf-8') as csvfile:
            w = writer(csvfile)
            w.writerow(['word', 'count', 'frequency'])
            for row in rowdata:
                towrite = [row[0], row[1], row[1] * 100.0 / n_elements]
                w.writerow(towrite)

    """
        Based on
        https://medium.com/@nicharuch/collocations-identifying-phrases-that-act-like-individual-words-in-nlp-f58a93a2f84a
    """

    def right_bigram_type(self, ngram):
        if '-pron-' in ngram or 't' in ngram:
            return False
        for word in ngram:
            if word in self.stop or word.isspace() or word in self.bigram_words_to_ignore:
                return False
        acceptable_types = ('JJ', 'JJR', 'JJS', 'NN', 'NNS', 'NNP', 'NNPS')
        second_type = ('NN', 'NNS', 'NNP', 'NNPS')
        tags = nltk.pos_tag(ngram)
        if tags[0][1] in acceptable_types and tags[1][1] in second_type:
            # Further checks
            return ngram not in self.bigrams_to_ignore
        else:
            return False

    def calculate_bigram_frequencies(self):

        finder = nltk.collocations.BigramCollocationFinder.from_words(self._words)
        bigram_freq = finder.ngram_fd.items()
        frequency_table = pd.DataFrame(list(bigram_freq), columns=['bigram', 'freq']).sort_values(by='freq',
                                                                                                  ascending=False)

        return frequency_table[frequency_table.bigram.map(lambda x: self.right_bigram_type(x))]


def make_word_cloud_from_text(text, fname):
    # lower max_font_size
    wc = WordCloud(max_font_size=40).generate(text)
    wc.to_file(fname)


def make_word_cloud_from_frequencies(dataframe, min_count, fname):
    d = {}
    for a, x in dataframe.values:
        if x >= min_count:
            d[" ".join(a)] = x

    wc = WordCloud(width=800, height=600, background_color='white', prefer_horizontal=0.5)
    wc.generate_from_frequencies(frequencies=d)
    wc.to_file(fname)


def word_frequency_analysis(min_bigrams, text_to_process, fname):
    c = FrequencyCalculator(text_to_process)
    make_word_cloud_from_frequencies(c.calculate_bigram_frequencies(), min_bigrams, fname)

import attr
import numpy as np

from pycortex.utils import revcomp

EDGE_SET_LENGTH = 8
EDGE_SET_LETTER_LOOKUP = {
    'a': 0, 'c': 1, 'g': 2, 't': 3,
    'A': 4, 'C': 5, 'G': 6, 'T': 7,
}
EDGE_IDX_TO_LETTER = ['A', 'C', 'G', 'T']


def lexlo(kmer_string):
    """Returns the lexicographically lowest version of the kmer string"""
    alt_kmer_string = revcomp(kmer_string)
    if alt_kmer_string < kmer_string:
        return alt_kmer_string
    return kmer_string


@attr.s(slots=True, cmp=False)
class EdgeSet(object):
    """Adds methods for accessing an edge set array (data)"""
    data = attr.ib()

    @data.validator
    def check(self, _, value):
        assert len(value) == EDGE_SET_LENGTH

    def is_edge(self, letter):
        return self.data[EDGE_SET_LETTER_LOOKUP[letter]]

    def add_edge(self, letter):
        self.data[EDGE_SET_LETTER_LOOKUP[letter]] = 1

    def remove_edge(self, letter):
        self.data[EDGE_SET_LETTER_LOOKUP[letter]] = 0

    def __getitem__(self, item):
        return self.data[item]

    def __eq__(self, other):
        return np.array_equal(self.data, other.data)

    @property
    def outgoing(self):
        return self.data[EDGE_SET_LENGTH // 2:]

    @property
    def incoming(self):
        return self.data[:EDGE_SET_LENGTH // 2]

    def get_incoming_kmers(self, kmer_string):
        kmers = []
        kmer_suffix = kmer_string[:-1]
        for edge_idx, edge in enumerate(self.incoming):
            if edge:
                kmers.append(lexlo(EDGE_IDX_TO_LETTER[edge_idx] + kmer_suffix))
        return kmers

    def get_outgoing_kmers(self, kmer_string):
        kmers = []
        kmer_prefix = kmer_string[1:]
        for edge_idx, edge in enumerate(self.outgoing):
            if edge:
                kmers.append(lexlo(kmer_prefix + EDGE_IDX_TO_LETTER[edge_idx]))
        return kmers

    def to_str(self, *, as_revcomp=False):
        es_letters = []
        for letter, index in EDGE_SET_LETTER_LOOKUP.items():
            if self.is_edge(letter):
                es_letters.append(letter)
            else:
                es_letters.append('.')
        es_letters = ''.join(es_letters)
        if as_revcomp:
            es_letters = revcomp(es_letters)
        return es_letters


def empty():
    return EdgeSet(np.zeros(EDGE_SET_LENGTH, dtype=np.uint8))

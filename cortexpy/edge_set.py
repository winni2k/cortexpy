import attr
import numpy as np

from cortexpy.utils import revcomp, lexlo

EDGE_SET_LENGTH = 8
EDGE_SET_ORDERED_LETTERS = 'acgtACGT'
EDGE_SET_LETTER_LOOKUP = {letter: idx for idx, letter in enumerate(EDGE_SET_ORDERED_LETTERS)}

EDGE_IDX_TO_LETTER = ['A', 'C', 'G', 'T']


@attr.s(slots=True, cmp=False)
class EdgeSet(object):
    """Adds methods for accessing an edge set array (data)"""
    data = attr.ib(attr.Factory(lambda: np.zeros(EDGE_SET_LENGTH, dtype=np.uint8)))

    @data.validator
    def check(self, _, value):  # noqa
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
        for letter in EDGE_SET_ORDERED_LETTERS:
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

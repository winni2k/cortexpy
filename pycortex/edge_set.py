import attr
import numpy as np

EDGE_SET_LENGTH = 8
EDGE_SET_LETTER_LOOKUP = {
    'a': 0, 'c': 1, 'g': 2, 't': 3,
    'A': 4, 'C': 5, 'G': 6, 'T': 7,
}
EDGE_IDX_TO_LETTER = ['A', 'C', 'G', 'T']


@attr.s(slots=True, cmp=False)
class EdgeSet(object):
    """Adds methods for accessing an edge set array (data)"""
    data = attr.ib()

    @data.validator
    def check(self, _, value):
        assert len(value) == EDGE_SET_LENGTH

    def is_edge(self, letter):
        return self.data[EDGE_SET_LETTER_LOOKUP[letter]]

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
                kmers.append(EDGE_IDX_TO_LETTER[edge_idx] + kmer_suffix)
        return kmers

    def get_outgoing_kmers(self, kmer_string):
        kmers = []
        kmer_prefix = kmer_string[1:]
        for edge_idx, edge in enumerate(self.outgoing):
            if edge:
                kmers.append(kmer_prefix + EDGE_IDX_TO_LETTER[edge_idx])
        return kmers

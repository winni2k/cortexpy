import attr
import numpy as np

EDGE_SET_LETTER_LOOKUP = {
    'a': 0, 'c': 1, 'g': 2, 't': 3,
    'A': 4, 'C': 5, 'G': 6, 'T': 7,
}
EDGE_IDX_TO_LETTER = ['A', 'C', 'G', 'T']


@attr.s(slots=True, cmp=False)
class EdgeSet(object):
    edges = attr.ib()

    def __attrs_post_init__(self):
        if isinstance(self.edges, str):
            self.edges = edge_set_string_to_array(self.edges)

    @edges.validator
    def check(self, _, value):
        assert len(value) == 8

    def is_edge(self, letter):
        return self.edges[EDGE_SET_LETTER_LOOKUP[letter]]

    def __getitem__(self, item):
        return self.edges[item]

    def __eq__(self, other):
        return np.array_equal(self.edges, other.edges)

    @property
    def outgoing(self):
        return self.edges[4:]

    @property
    def incoming(self):
        return self.edges[:4]

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


def edge_set_string_to_array(edge_set_string):
    assert len(edge_set_string) == 8
    edge_set = []
    for edge in edge_set_string:
        if edge == '.':
            edge_set.append(0)
        else:
            edge_set.append(1)
    return np.array(edge_set, dtype=np.uint8)

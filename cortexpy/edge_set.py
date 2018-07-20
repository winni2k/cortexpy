import attr
import numpy as np

from cortexpy.constants import EdgeTraversalOrientation
from cortexpy.utils import revcomp, lexlo

EDGE_SET_LENGTH = 8
HALF_EDGE_SET_LENGTH = 4
EDGE_SET_DATA_LETTER_ORDER = 'acgtTGCA'
EDGE_SET_REPR_LETTER_ORDER = 'acgtACGT'
EDGE_SET_LETTER_LOOKUP = {letter: idx for idx, letter in enumerate(EDGE_SET_DATA_LETTER_ORDER)}

EDGE_IDX_TO_LETTER = ['A', 'C', 'G', 'T', 'T', 'G', 'C', 'A']


@attr.s(slots=True, cmp=False)
class EdgeSet(object):
    """Adds methods for accessing an edge set array (data)"""
    data = attr.ib()

    @data.validator
    def check(self, _, value):  # noqa
        assert len(value) == EDGE_SET_LENGTH

    def is_edge(self, letter):
        return self.data[EDGE_SET_LETTER_LOOKUP[letter]] == 1

    def add_edge(self, letter):
        data = list(self.data)
        data[EDGE_SET_LETTER_LOOKUP[letter]] = 1
        self.data = tuple(data)

    def remove_edge(self, letter):
        data = list(self.data)
        data[EDGE_SET_LETTER_LOOKUP[letter]] = 0
        self.data = tuple(data)

    def __getitem__(self, item):
        return self.data[item]

    def __eq__(self, other):
        return all(s == o for s, o in zip(self.data, other.data))

    def __str__(self):
        return self.to_str(as_revcomp=False)

    @property
    def outgoing(self):
        return self.data[HALF_EDGE_SET_LENGTH:]

    @property
    def incoming(self):
        return self.data[:HALF_EDGE_SET_LENGTH]

    def num_outgoing(self):
        return sum(self.outgoing)

    def num_incoming(self):
        return sum(self.incoming)

    def _get_kmer_strings(self, sub_kmer_string, is_incoming, is_lexlo):
        if is_incoming:
            if is_lexlo:
                edges = self.incoming
            else:
                edges = self.outgoing
            for edge_idx, edge in enumerate(edges):
                if edge:
                    yield EDGE_IDX_TO_LETTER[edge_idx] + sub_kmer_string
        else:
            if is_lexlo:
                edges = self.outgoing
            else:
                edges = self.incoming
            for edge_idx, edge in enumerate(edges):
                if edge:
                    yield sub_kmer_string + EDGE_IDX_TO_LETTER[edge_idx + 4]

    def get_incoming_kmer_strings(self, kmer_string, is_lexlo=None):
        if is_lexlo is None:
            is_lexlo = bool(kmer_string == lexlo(kmer_string))
        return self._get_kmer_strings(kmer_string[:-1], True, is_lexlo)

    def get_outgoing_kmer_strings(self, kmer_string, is_lexlo=None):
        if is_lexlo is None:
            is_lexlo = bool(kmer_string == lexlo(kmer_string))
        return self._get_kmer_strings(kmer_string[1:], False, is_lexlo)

    def get_incoming_kmers(self, kmer_string):
        lexlo_string = lexlo(kmer_string)
        assert lexlo_string == kmer_string
        return [lexlo(kmer_string) for kmer_string in
                self.get_incoming_kmer_strings(kmer_string, is_lexlo=True)]

    def get_outgoing_kmers(self, kmer_string):
        lexlo_string = lexlo(kmer_string)
        assert lexlo_string == kmer_string
        return [lexlo(kmer_string) for kmer_string in
                self.get_outgoing_kmer_strings(kmer_string, is_lexlo=True)]

    def to_str(self, *, as_revcomp=False):
        es_letters = []
        for letter in EDGE_SET_REPR_LETTER_ORDER:
            if self.is_edge(letter):
                es_letters.append(letter)
            else:
                es_letters.append('.')
        es_letters = ''.join(es_letters)
        if as_revcomp:
            es_letters = revcomp(es_letters)
        return es_letters

    def oriented(self, orientation):
        return OrientedEdgeSet(self, orientation)

    def dump(self, buffer):
        binary = np.packbits(self.data)
        assert 1 == len(binary)
        buffer.write(binary[0])


def empty():
    return EdgeSet(np.zeros(EDGE_SET_LENGTH, dtype=np.uint8))


@attr.s(slots=True)
class OrientedEdgeSet(object):
    edge_set = attr.ib()
    orientation = attr.ib()
    neighbor_kmers = attr.ib(init=False)
    neighbor = attr.ib(init=False)
    neighbor_kmer_strings = attr.ib(init=False)
    _num_neighbor = attr.ib(init=False)

    def __attrs_post_init__(self):
        if self.orientation == EdgeTraversalOrientation.original:
            self.neighbor_kmers = self.edge_set.get_outgoing_kmers
            self.neighbor = self.edge_set.outgoing
            self._num_neighbor = self.edge_set.num_outgoing
            self.neighbor_kmer_strings = self.edge_set.get_outgoing_kmer_strings
        else:
            self.neighbor_kmers = self.edge_set.get_incoming_kmers
            self.neighbor = self.edge_set.incoming
            self._num_neighbor = self.edge_set.num_incoming
            self.neighbor_kmer_strings = self.edge_set.get_incoming_kmer_strings

    def other_orientation(self):
        return OrientedEdgeSet(self.edge_set, EdgeTraversalOrientation.other(self.orientation))

    def num_neighbor(self, kmer_string):
        lexlo_kmer_string = lexlo(kmer_string)
        if lexlo_kmer_string != kmer_string:
            return self.other_orientation().num_neighbor(lexlo(kmer_string))
        else:
            return self._num_neighbor()

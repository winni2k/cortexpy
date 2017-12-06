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

    def num_outgoing(self):
        return sum(self.outgoing)

    def num_incoming(self):
        return sum(self.incoming)

    def _get_kmer_strings(self, sub_kmer_string, is_incoming, is_lexlo):
        kmers = []
        edge_funcs = [self.incoming, self.outgoing]
        if is_incoming:
            edge_funcs = list(reversed(edge_funcs))
        if not is_lexlo:
            edges = edge_funcs[0][::-1]
        else:
            edges = edge_funcs[1]
        for edge_idx, edge in enumerate(edges):
            if edge:
                if is_incoming:
                    new_kmer_string = EDGE_IDX_TO_LETTER[edge_idx] + sub_kmer_string
                else:
                    new_kmer_string = sub_kmer_string + EDGE_IDX_TO_LETTER[edge_idx]
                kmers.append(new_kmer_string)
        return kmers

    def get_incoming_kmer_strings(self, kmer_string, is_lexlo=None):
        if is_lexlo is None:
            is_lexlo = bool(kmer_string == lexlo(kmer_string))
        sub_kmer_string = kmer_string[:-1]
        return self._get_kmer_strings(sub_kmer_string, True, is_lexlo)

    def get_outgoing_kmer_strings(self, kmer_string, is_lexlo=None):
        if is_lexlo is None:
            is_lexlo = bool(kmer_string == lexlo(kmer_string))
        sub_kmer_string = kmer_string[1:]
        return self._get_kmer_strings(sub_kmer_string, False, is_lexlo)

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

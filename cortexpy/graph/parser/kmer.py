from struct import unpack
import attr
import numpy as np
import math

import cortexpy.edge_set
from cortexpy.edge_set import EdgeSet
from cortexpy.graph.parser.constants import (
    NUM_TO_LETTER, UINT64_T, UINT32_T, LETTER_TO_NUM,
    LETTERS_PER_BYTE,
)
from cortexpy.utils import revcomp, lexlo


def check_kmer_string(kmer_string):
    if len(kmer_string) % 2 == 0:
        raise ValueError('kmer_string needs to be odd length')
    if len(kmer_string) < 3:
        raise ValueError('kmer_string needs to length 3 or more')


@attr.s(slots=True)
class EmptyKmerBuilder(object):
    num_colors = attr.ib(0)
    _seen_kmers = attr.ib(attr.Factory(dict))

    @num_colors.validator  # noqa
    def not_less_than_zero(self, attribute, value):  # noqa
        if value < 0:
            raise ValueError('value less than zero')

    def _build_from_lexlo(self, kmer_string, is_lexlo=False):
        """Build empty kmer from a lexicographically lowest string"""
        return Kmer(EmptyKmer(kmer=kmer_string,
                              coverage=np.zeros(self.num_colors,
                                                dtype=np.uint8),
                              kmer_size=len(kmer_string),
                              num_colors=self.num_colors))

    def build(self, kmer_string):
        """Build empty kmer from a kmer string"""
        check_kmer_string(kmer_string)
        kmer_string_to_use = lexlo(kmer_string)
        return self._build_from_lexlo(kmer_string_to_use)

    def build_or_get(self, kmer_string):
        """Build empty kmer or return a cached kmer for a kmer string"""
        check_kmer_string(kmer_string)
        kmer_string_to_use = lexlo(kmer_string)
        if kmer_string_to_use in self._seen_kmers.keys():
            return self._seen_kmers[kmer_string_to_use]
        kmer = self._build_from_lexlo(kmer_string_to_use)
        self._seen_kmers[kmer_string_to_use] = kmer
        return kmer


def flip_kmer_string_to_match(flip, ref, flip_is_after_reference_kmer=True):
    flip_revcomp = revcomp(flip)
    if flip_is_after_reference_kmer:
        ref_core = ref[1:]
        flip_core = flip[:-1]
        flip_revcomp_core = flip_revcomp[:-1]
    else:
        ref_core = ref[:-1]
        flip_core = flip[1:]
        flip_revcomp_core = flip_revcomp[1:]

    if ref_core == flip_core:
        return flip, False
    elif ref_core == flip_revcomp_core:
        return flip_revcomp, True
    else:
        raise ValueError


def connect_kmers(first, second, color, identical_kmer_check=True):
    if identical_kmer_check and first == second and first is not second:
        raise ValueError('Kmers are equal, but not the same object')
    connection_is_set = False
    for flip_is_after_reference_kmer in [True, False]:
        for reverse_first_second in [True, False]:
            if reverse_first_second:
                flip_kmer, ref_kmer = second, first
            else:
                flip_kmer, ref_kmer = first, second
            try:
                flipped_string, is_flipped = flip_kmer_string_to_match(
                    flip_kmer.kmer, ref_kmer.kmer,
                    flip_is_after_reference_kmer=flip_is_after_reference_kmer
                )
            except ValueError:
                pass
            else:
                connection_is_set = True
                if flip_is_after_reference_kmer:
                    ref_letter = flipped_string[-1]
                    flip_letter = ref_kmer.kmer[0].lower()
                else:
                    ref_letter = flipped_string[0].lower()
                    flip_letter = ref_kmer.kmer[-1]

                if is_flipped:
                    flip_letter = revcomp(flip_letter).swapcase()
                ref_kmer.edges[color].add_edge(ref_letter)
                flip_kmer.edges[color].add_edge(flip_letter)

    if not connection_is_set:
        raise ValueError(
            'first kmer ({}) cannot be connected to second kmer ({})'.format(first.kmer,
                                                                             second.kmer)
        )


def kmer_eq(self, other):
    if any([self.kmer != other.kmer,
            self.kmer_size != other.kmer_size,
            self.num_colors != other.num_colors,
            not np.all(self.coverage == other.coverage),
            (self.edges is None) != (other.edges is None),
            (self.edges is not None and not np.all(self.edges == other.edges))]):
        return False
    return True


@attr.s(slots=True, cmp=False)
class EmptyKmer(object):
    kmer = attr.ib()
    coverage = attr.ib()
    kmer_size = attr.ib()
    num_colors = attr.ib()
    edges = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.edges = [cortexpy.edge_set.empty() for _ in range(len(self.coverage))]


@attr.s(slots=True)
class RawKmerConverter(object):
    kmer_size = attr.ib()

    def to_letters(self, raw_kmer):
        kmer_as_uint64ts = np.frombuffer(raw_kmer, dtype='<u8')
        big_endian_kmer = kmer_as_uint64ts.astype('>u8')
        kmer_as_bits = np.unpackbits(np.frombuffer(big_endian_kmer.tobytes(), dtype=np.uint8))
        kmer = (kmer_as_bits.reshape(-1, 2) * np.array([2, 1])).sum(1)
        return NUM_TO_LETTER[kmer[(len(kmer) - self.kmer_size):]]


NUM_TO_BITS = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])


@attr.s(slots=True)
class StringKmerConverter(object):
    kmer_size = attr.ib()
    _kmer_container_size_in_uint64ts = attr.ib(init=False)
    _padding_array = attr.ib(init=False)

    def __attrs_post_init__(self):
        num_letters_per_uint = UINT64_T * LETTERS_PER_BYTE
        self._kmer_container_size_in_uint64ts = math.ceil(
            self.kmer_size / num_letters_per_uint)
        padding_size = num_letters_per_uint - self.kmer_size % num_letters_per_uint
        if padding_size == num_letters_per_uint:
            padding_size = 0
        self._padding_array = np.zeros(padding_size, dtype=np.uint8)

    def to_uints(self, kmer_string):
        """Converts kmer_string to big-endian uint64 array"""
        assert isinstance(kmer_string, str)
        encoded_kmer_string = kmer_string.encode()
        translated_kmer_string = encoded_kmer_string.translate(LETTER_TO_NUM)
        letter_vals = np.frombuffer(translated_kmer_string, dtype=np.uint8)
        letter_vals = np.concatenate((letter_vals[:-self.kmer_size],
                                      self._padding_array,
                                      letter_vals[-self.kmer_size:]))
        letter_val_bits = NUM_TO_BITS[letter_vals]
        return np.packbits(letter_val_bits).view('uint64').newbyteorder()


@attr.s(slots=True, cmp=False)
class KmerData(object):
    _data = attr.ib()
    kmer_size = attr.ib()
    num_colors = attr.ib()
    _kmer_container_size_in_uint64ts = attr.ib(1)
    _kmer = attr.ib(None)
    _coverage = attr.ib(None)
    _edges = attr.ib(None)
    _incoming_kmers = attr.ib(None)
    _outgoing_kmers = attr.ib(None)
    _kmer_vals_to_delete = attr.ib(init=False)

    def __attrs_post_init__(self):
        n_vals_left_over = self.kmer_size % 4
        n_vals_to_remove = 4 - n_vals_left_over
        if n_vals_to_remove > 0:
            object.__setattr__(self, "_kmer_vals_to_delete",
                               np.arange(0, n_vals_to_remove) + self.kmer_size - n_vals_left_over)

    @property
    def kmer(self):
        if self._kmer is None:
            raw_kmer = self._data[:self._kmer_container_size_in_uint64ts * UINT64_T]
            kmer_letters = RawKmerConverter(self.kmer_size).to_letters(raw_kmer)
            object.__setattr__(self, "_kmer", kmer_letters.astype('|S1').tostring().decode('utf-8'))
        return self._kmer

    @property
    def coverage(self):
        if self._coverage is None:
            start = self._kmer_container_size_in_uint64ts * UINT64_T
            coverage_raw = self._data[start:(start + self.num_colors * UINT32_T)]
            fmt_string = ''.join(['I' for _ in range(self.num_colors)])
            object.__setattr__(self, "_coverage",
                               np.array(unpack(fmt_string, coverage_raw), dtype=np.uint32))
        return self._coverage

    @property
    def edges(self):
        if self._edges is None:
            start = (
                self._kmer_container_size_in_uint64ts * UINT64_T + self.num_colors * UINT32_T
            )
            edge_bytes = np.frombuffer(self._data[start:], dtype=np.uint8)
            edge_sets = np.unpackbits(edge_bytes)
            edge_sets = edge_sets.reshape(-1, 8)
            edge_sets = map(EdgeSet, edge_sets)
            edge_sets = list(edge_sets)
            object.__setattr__(self, "_edges", edge_sets)
        return self._edges


@attr.s(slots=True, cmp=False)
class Kmer(object):
    _kmer_data = attr.ib()
    _kmer = attr.ib(None)
    _coverage = attr.ib(None)
    _edges = attr.ib(None)
    kmer_size = attr.ib(init=False)
    num_colors = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.kmer_size = self._kmer_data.kmer_size
        self.num_colors = self._kmer_data.num_colors

    @property
    def kmer(self):
        if self._kmer is None:
            self._kmer = self._kmer_data.kmer
        return self._kmer

    @kmer.setter
    def kmer(self, val):
        self._kmer = val

    @property
    def coverage(self):
        if self._coverage is None:
            self._coverage = self._kmer_data.coverage
        return self._coverage

    @coverage.setter
    def coverage(self, val):
        self._coverage = val

    @property
    def edges(self):
        if self._edges is None:
            self._edges = self._kmer_data.edges
        return self._edges

    @edges.setter
    def edges(self, val):
        return self._edges

    def increment_color_coverage(self, color):
        """Increment the coverage of a color by one"""
        num_colors = max(self.num_colors, color + 1)
        if num_colors > self.num_colors:
            coverage_array = np.append(self.coverage, np.zeros(num_colors - self.num_colors,
                                                               dtype=self.coverage.dtype))
            self.coverage = coverage_array
            for edge_idx in range(self.num_colors, num_colors):
                assert len(self.edges) < edge_idx + 1
                self.edges.append(cortexpy.edge_set.empty())
        self.coverage[color] += 1
        self.num_colors = num_colors

    def __eq__(self, other):
        return kmer_eq(self, other)

    def __str__(self):
        string_parts = [self.kmer]
        string_parts += [str(c) for c in self.coverage]
        string_parts += [e.to_str() for e in self.edges]
        return ' '.join(string_parts)

    @property
    def colors(self):
        return range(self.num_colors)

    def has_outgoing_edge_to_kmer_in_color(self, other, color):
        if self.kmer[1:] == other.kmer[:-1]:
            other_kmer_letter = other.kmer[-1]
            this_kmer_letter = self.kmer[0].lower()
        else:
            other_kmer_revcomp = revcomp(other.kmer)
            if self.kmer[1:] == other_kmer_revcomp[:-1]:
                other_kmer_letter = other_kmer_revcomp[-1]
                this_kmer_letter = revcomp(self.kmer[0])
            else:
                raise ValueError('Kmers are not neighbors')
        edge_set = self.edges[color]
        if edge_set.is_edge(other_kmer_letter) != other.edges[color].is_edge(this_kmer_letter):
            raise ValueError('Kmers ({}) and ({}) do not agree on connection'.format(self, other))
        return edge_set.is_edge(other_kmer_letter)

    def has_incoming_edge_from_kmer_in_color(self, other, color):
        if self.kmer[:-1] == other.kmer[1:]:
            other_kmer_letter = other.kmer[0].lower()
            this_kmer_letter = self.kmer[-1]
        else:
            other_kmer_revcomp = revcomp(other.kmer)
            if self.kmer[:-1] == other_kmer_revcomp[1:]:
                other_kmer_letter = other_kmer_revcomp[0].lower()
                this_kmer_letter = revcomp(self.kmer[-1]).lower()
            else:
                raise ValueError('Kmers are not neighbors')
        edge_set = self.edges[color]
        if edge_set.is_edge(other_kmer_letter) != other.edges[color].is_edge(this_kmer_letter):
            raise ValueError('Kmers do not agree on connection')
        return edge_set.is_edge(other_kmer_letter)


@attr.s(slots=True, cmp=False)
class KmerByStringComparator(object):
    kmer = attr.ib(None)
    kmer_object = attr.ib(None)

    def __attrs_post_init__(self):
        if self.kmer is None:
            self.kmer = self.kmer_object.kmer

    def __eq__(self, other):
        return self.kmer == other.kmer

    def __lt__(self, other):
        return self.kmer < other.kmer


@attr.s(slots=True, cmp=False)
class KmerUintComparator(object):
    kmer_uints = attr.ib()

    def __lt__(self, other):
        for i, val in enumerate(self.kmer_uints):
            if val < other.kmer_uints[i]:
                return True
            elif val != other.kmer_uints[i]:
                return False
        return False

    def __eq__(self, other):
        return np.all(self.kmer_uints == other.kmer_uints)

    def __len__(self):
        return len(self.kmer_uints)

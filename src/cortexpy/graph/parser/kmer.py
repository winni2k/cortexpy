"""Cortex kmers
===============

This module provides classes and functions for working with Cortex kmers.
"""

import math
import struct
from itertools import repeat

import attr
import numpy as np

import cortexpy.edge_set
from cortexpy.utils import revcomp, lexlo
from .constants import (
    UINT64_T, UINT32_T, LETTER_TO_NUM,
    NUM_LETTERS_PER_UINT, NUM_TO_BITS,
)
from .kmer_ext import raw_kmer_to_string, raw_edges_to_list, raw_to_coverage


def check_kmer_string(kmer_string):
    if len(kmer_string) % 2 == 0:
        raise ValueError('kmer_string needs to be odd length')
    if len(kmer_string) < 3:
        raise ValueError('kmer_string needs to length 3 or more')


@attr.s(slots=True)
class EmptyKmerBuilder(object):
    num_colors = attr.ib(1)
    default_coverage = attr.ib(1)
    _seen_kmers = attr.ib(attr.Factory(dict))

    @num_colors.validator  # noqa
    def not_less_than_zero(self, attribute, value):  # noqa
        if value < 0:
            raise ValueError('value less than zero')

    def _build_from_lexlo(self, kmer_string):
        """Build empty kmer from a lexicographically lowest string"""
        return Kmer.from_kmer_data(EmptyKmer(kmer=kmer_string,
                                             coverage=tuple(
                                                 repeat(self.default_coverage, self.num_colors)),
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


def revcomp_target_to_match_ref(target, ref, rc_is_after_reference_kmer=True):
    target_revcomp = revcomp(target)
    if rc_is_after_reference_kmer:
        ref_core = ref[1:]
        flip_core = target[:-1]
        flip_revcomp_core = target_revcomp[:-1]
    else:
        ref_core = ref[:-1]
        flip_core = target[1:]
        flip_revcomp_core = target_revcomp[1:]

    if ref_core == flip_core:
        return target, False
    elif ref_core == flip_revcomp_core:
        return target_revcomp, True
    else:
        raise ValueError(
            "ref_core ({}) does not match"
            " flip_core ({}) or flip_revcomp_core ({})".format(ref_core,
                                                               flip_core,
                                                               flip_revcomp_core))


def find_neighbors(first, second, *, rc_is_after_reference_kmer, revcomp_second):
    if revcomp_second:
        revcomp_kmer, ref_kmer = second, first
    else:
        revcomp_kmer, ref_kmer = first, second
    try:
        revcomp_string, is_revcomp = revcomp_target_to_match_ref(
            revcomp_kmer.kmer, ref_kmer.kmer,
            rc_is_after_reference_kmer=rc_is_after_reference_kmer
        )
    except ValueError:
        pass
    else:
        if rc_is_after_reference_kmer:
            ref_letter = revcomp_string[-1]
            flip_letter = ref_kmer.kmer[0].lower()
        else:
            ref_letter = revcomp_string[0].lower()
            flip_letter = ref_kmer.kmer[-1]

        if is_revcomp:
            flip_letter = revcomp(flip_letter).swapcase()
        yield ref_kmer, revcomp_kmer, ref_letter, flip_letter


def find_all_neighbors(first, second):
    """Return kmers and letters to get from first kmer to second"""
    for rc_is_after_reference_kmer in [True, False]:
        for revcomp_second in [True, False]:
            yield from find_neighbors(first, second,
                                      rc_is_after_reference_kmer=rc_is_after_reference_kmer,
                                      revcomp_second=revcomp_second)


def connect_kmers(first, second, color, identical_kmer_check=True):
    """Connect two kmers"""
    if identical_kmer_check and first == second and first is not second:
        raise ValueError('Kmers are equal, but not the same object')
    are_neighbors = False
    for ref_kmer, flip_kmer, ref_letter, flip_letter in find_all_neighbors(first, second):
        are_neighbors = True
        ref_kmer.edges[color].add_edge(ref_letter)
        flip_kmer.edges[color].add_edge(flip_letter)
    if not are_neighbors:
        raise ValueError(
            'first kmer ({}) cannot be connected to second kmer ({})'.format(first.kmer,
                                                                             second.kmer)
        )


def disconnect_kmers(first, second, colors):
    """Disconnect two kmers"""
    are_neighbors = False
    for ref_kmer, flip_kmer, ref_letter, flip_letter in find_all_neighbors(first, second):
        are_neighbors = True
        for color in colors:
            ref_kmer.edges[color].remove_edge(ref_letter)
            flip_kmer.edges[color].remove_edge(flip_letter)
    if not are_neighbors:
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

    def get_raw_kmer(self):
        return StringKmerConverter(self.kmer_size).to_raw(self.kmer)


@attr.s(slots=True)
class RawKmerConverter(object):
    kmer_size = attr.ib()

    def to_string(self, raw_kmer):
        return raw_kmer_to_string(self.kmer_size, raw_kmer)


def calc_kmer_container_size(kmer_size):
    return math.ceil(kmer_size / NUM_LETTERS_PER_UINT)


@attr.s(slots=True)
class StringKmerConverter(object):
    """Converts kmer strings to various binary representations"""
    kmer_size = attr.ib()
    _kmer_container_size_in_uint64ts = attr.ib(init=False)
    _padding_array = attr.ib(init=False)

    def __attrs_post_init__(self):
        self._kmer_container_size_in_uint64ts = calc_kmer_container_size(self.kmer_size)
        padding_size = NUM_LETTERS_PER_UINT - self.kmer_size % NUM_LETTERS_PER_UINT
        if padding_size == NUM_LETTERS_PER_UINT:
            padding_size = 0
        self._padding_array = np.zeros(padding_size, dtype=np.uint8)

    def to_uints(self, kmer_string):
        """Converts kmer_string to big-endian uint64 array"""
        kmer_string = str(kmer_string)
        encoded_kmer_string = kmer_string.encode()
        translated_kmer_string = encoded_kmer_string.translate(LETTER_TO_NUM)
        letter_vals = np.frombuffer(translated_kmer_string, dtype=np.uint8)
        letter_vals = np.concatenate((letter_vals[:-self.kmer_size],
                                      self._padding_array,
                                      letter_vals[-self.kmer_size:]))
        letter_val_bits = NUM_TO_BITS[letter_vals]
        return np.packbits(letter_val_bits).view('uint64').newbyteorder()

    def to_raw(self, kmer_string):
        uints = self.to_uints(kmer_string)
        little_endian_uints = uints.astype('<u8')
        return little_endian_uints.tobytes()


@attr.s(slots=True, cmp=False)
class KmerData(object):
    _data = attr.ib()
    kmer_size = attr.ib()
    num_colors = attr.ib()
    _kmer = attr.ib(None, init=False)
    _coverage = attr.ib(None, init=False)
    _edges = attr.ib(None, init=False)

    def get_raw_kmer(self):
        return self._data[:self.kmer_container_size_in_uint64ts * UINT64_T]

    @property
    def kmer(self):
        if self._kmer is None:
            self._kmer = raw_kmer_to_string(self.kmer_size, self.get_raw_kmer())
        return self._kmer

    @property
    def coverage(self):
        if self._coverage is None:
            self._coverage = raw_to_coverage(self._data,
                                             self.kmer_container_size_in_uint64ts * UINT64_T,
                                             self.num_colors)
        return self._coverage

    @coverage.setter
    def coverage(self, val):
        self._coverage = val

    @property
    def edges(self):
        if self._edges is None:
            start = (
                self.kmer_container_size_in_uint64ts * UINT64_T + self.num_colors * UINT32_T
            )
            tuples = raw_edges_to_list(self._data[start:])
            self._edges = [cortexpy.edge_set.EdgeSet(t) for t in tuples]
        return self._edges

    @edges.setter
    def edges(self, val):
        self._edges = val

    @property
    def kmer_container_size_in_uint64ts(self):
        return calc_kmer_container_size(self.kmer_size)


@attr.s(slots=True, cmp=False)
class Kmer:
    """Represents a Cortex kmer

    This class wraps a kmer data object with attributes and methods for inspecting and manipulating
    the underlying kmer data object."""
    _kmer_data = attr.ib()
    num_colors = attr.ib()
    kmer_size = attr.ib()
    _revcomp = attr.ib(None)

    @num_colors.validator
    def check(self, attribute, value):
        if value < 1:
            raise ValueError("num_colors must be greater than 0")

    @classmethod
    def from_kmer_data(cls, kmer_data):
        return cls(kmer_data, kmer_size=kmer_data.kmer_size, num_colors=kmer_data.num_colors)

    @property
    def revcomp(self):
        if self._revcomp is None:
            self._revcomp = revcomp(self.kmer)
        return self._revcomp

    @property
    def kmer_container_size(self):
        return self._kmer_data.kmer_container_size_in_uint64ts

    @property
    def kmer(self):
        return self._kmer_data.kmer

    @property
    def coverage(self):
        return self._kmer_data.coverage

    @coverage.setter
    def coverage(self, value):
        assert isinstance(value, tuple)
        if len(value) != self.num_colors:
            raise ValueError(f'coverage length ({value}) must match num_colors ({self.num_colors})')
        assert any(c > 0 for c in value)
        self._kmer_data.coverage = value

    @property
    def edges(self):
        return self._kmer_data.edges

    @edges.setter
    def edges(self, val):
        self._kmer_data.edges = val

    def increment_color_coverage(self, color):
        """Increment the coverage of a color by one"""
        num_colors = max(self.num_colors, color + 1)
        coverage = list(self.coverage)
        if num_colors > self.num_colors:
            coverage += list(repeat(0, num_colors - self.num_colors))
            for edge_idx in range(self.num_colors, num_colors):
                assert len(self.edges) < edge_idx + 1
                self.edges.append(cortexpy.edge_set.empty())
        coverage[color] += 1
        self.num_colors = num_colors
        self.coverage = tuple(coverage)

    def __eq__(self, other):
        return kmer_eq(self, other)

    def __str__(self):
        string_parts = [self.kmer]
        string_parts += [str(c) for c in self.coverage]
        string_parts += [e.to_str() for e in self.edges]
        return ' '.join(string_parts)

    def __repr__(self):
        return str(self)

    @property
    def colors(self):
        return range(self.num_colors)

    def has_outgoing_edge_to_kmer_in_color(self, other, color):
        try:
            _, _, ref_letter, flip_letter = next(find_neighbors(self, other,
                                                                revcomp_second=True,
                                                                rc_is_after_reference_kmer=True))
        except StopIteration:
            raise ValueError('Kmers ({}) and ({}) do not agree on connection'.format(self, other))
        edge_set = self.edges[color]
        if edge_set.is_edge(ref_letter) != other.edges[color].is_edge(flip_letter):
            raise ValueError('Kmers ({}) and ({}) do not agree on connection'.format(self, other))
        return edge_set.is_edge(ref_letter)

    def has_incoming_edge_from_kmer_in_color(self, other, color):
        try:
            _, _, ref_letter, flip_letter = next(find_neighbors(self, other,
                                                                revcomp_second=True,
                                                                rc_is_after_reference_kmer=False))
        except StopIteration:
            raise ValueError('Kmers ({}) and ({}) do not agree on connection'.format(self, other))
        edge_set = self.edges[color]
        if edge_set.is_edge(ref_letter) != other.edges[color].is_edge(flip_letter):
            raise ValueError('Kmers do not agree on connection')
        return edge_set.is_edge(ref_letter)

    def get_raw_kmer(self):
        return self._kmer_data.get_raw_kmer()

    def dump(self, buffer):
        if self.kmer == self._kmer_data.kmer:
            buffer.write(self.get_raw_kmer())
        else:
            raise NotImplementedError
        for cov in self.coverage:
            buffer.write(struct.pack('I', cov))
        for edge_set in self.edges:
            edge_set.dump(buffer)

    def __getitem__(self, item):
        """Allow kmers to be used as node objects in networkx graphs"""
        if item == 'kmer':
            return self
        raise KeyError


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

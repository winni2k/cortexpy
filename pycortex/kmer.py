from struct import unpack

import attr
import numpy as np

from pycortex.graph.parser.constants import NUM_TO_LETTER, UINT64_T, UINT32_T


@attr.s(slots=True)
class Kmer(object):
    _raw_data = attr.ib()
    kmer_size = attr.ib()
    num_colors = attr.ib()
    kmer_container_size_in_uint64ts = attr.ib(1)
    _kmer = attr.ib(None)
    _coverage = attr.ib(None)
    _edges = attr.ib(None)
    _kmer_vals_to_delete = attr.ib(init=False)

    def __attrs_post_init__(self):
        n_vals_left_over = self.kmer_size % 4
        n_vals_to_remove = 4 - n_vals_left_over
        if n_vals_to_remove > 0:
            self._kmer_vals_to_delete = (
                np.arange(0, n_vals_to_remove) + self.kmer_size - n_vals_left_over
            )

    @property
    def kmer(self):
        if self._kmer is None:
            kmer_as_uint64ts = np.frombuffer(
                self._raw_data[:self.kmer_container_size_in_uint64ts * 8],
                dtype='<u8')
            kmer_as_uint64ts_be = kmer_as_uint64ts.byteswap().newbyteorder()  # change to big endian
            kmer_as_properly_ordered_bits_right_aligned = np.unpackbits(
                np.frombuffer(kmer_as_uint64ts_be.tobytes(), dtype=np.uint8)
            )
            kmer = (
                kmer_as_properly_ordered_bits_right_aligned.reshape(-1, 2) * np.array([2, 1])
            ).sum(1)
            self._kmer = ''.join(NUM_TO_LETTER[num] for num in kmer[(len(kmer) - self.kmer_size):])
        return self._kmer

    @property
    def coverage(self):
        if self._coverage is None:
            start = self.kmer_container_size_in_uint64ts * UINT64_T
            coverage_raw = self._raw_data[start:(start + self.num_colors * UINT32_T)]
            fmt_string = ''.join(['I' for _ in range(self.num_colors)])
            self._coverage = unpack(fmt_string, coverage_raw)
        return self._coverage

    @property
    def edges(self):
        if self._edges is None:
            start = (
                self.kmer_container_size_in_uint64ts * UINT64_T + self.num_colors * UINT32_T
            )
            edge_bytes = list(self._raw_data[start:])
            edge_sets = np.unpackbits(np.array(edge_bytes, dtype=np.uint8)).reshape(-1, 4)
            edge_sets[1::2] = np.fliplr(edge_sets[1::2])
            self._edges = tuple(map(tuple, edge_sets.reshape(-1, 8)))
        return self._edges


class KmerByStringComparator(object):
    def __init__(self, *, kmer=None, kmer_object=None):
        self.kmer = kmer
        self.kmer_object = kmer_object
        if self.kmer is None:
            self.kmer = self.kmer_object.kmer

    def __eq__(self, other):
        return self.kmer == other.kmer

    def __lt__(self, other):
        return self.kmer < other.kmer

    def __repr__(self):
        return self.kmer

from bisect import bisect_left
from collections import Sequence, Mapping
from functools import lru_cache
from io import SEEK_END

import attr
import numpy as np

import cortexpy.graph.parser.header
from cortexpy.graph.parser.constants import UINT64_T
from cortexpy.graph.parser.kmer import (
    Kmer, KmerData, KmerUintComparator,
    StringKmerConverter,
)
from cortexpy.graph.parser.streaming import kmer_generator_from_stream_and_header
from cortexpy.utils import lexlo


@attr.s(slots=True)
class RandomAccess(Mapping):
    """Provide fast k-mer access to Cortex graph in log(n) time (n = number of kmers in graph)"""
    graph_handle = attr.ib()
    kmer_cache_size = attr.ib(None)
    kmer_cache_size_binary_search = attr.ib(None)
    header = attr.ib(init=False)
    graph_sequence = attr.ib(init=False)
    graph_kmer_sequence = attr.ib(init=False)
    n_records = attr.ib(init=False)
    _cached_get_kmer_data_for_string = attr.ib(init=False)

    def __attrs_post_init__(self):
        assert self.graph_handle.seekable()
        self.graph_handle.seek(0)
        self.header = cortexpy.graph.parser.header.from_stream(self.graph_handle)
        body_start_stream_position = self.graph_handle.tell()

        self.graph_handle.seek(0, SEEK_END)
        body_size = self.graph_handle.tell() - body_start_stream_position
        if body_size % self.header.record_size != 0:
            raise ValueError(
                "Body size ({}) % Record size ({}) != 0".format(body_size,
                                                                self.header.record_size))
        self.n_records = body_size // self.header.record_size
        if self.kmer_cache_size is None:
            self.kmer_cache_size = self.n_records
        if self.kmer_cache_size_binary_search is None:
            self.kmer_cache_size_binary_search = self.n_records
        self.graph_sequence = KmerRecordSequence(graph_handle=self.graph_handle,
                                                 body_start=body_start_stream_position,
                                                 header=self.header,
                                                 n_records=self.n_records,
                                                 kmer_cache_size=self.kmer_cache_size)
        self.graph_kmer_sequence = KmerUintSequence(
            graph_handle=self.graph_handle,
            body_start=body_start_stream_position,
            header=self.header,
            n_records=self.n_records,
            kmer_cache_size=self.kmer_cache_size_binary_search
        )

        self._cached_get_kmer_data_for_string = lru_cache(maxsize=self.kmer_cache_size)(
            self._get_kmer_data_for_string)

    def _get_kmer_data_for_string(self, kmer_string):
        uints = self.graph_kmer_sequence.kmer_string_converter.to_uints(kmer_string)
        index = self.graph_kmer_sequence.index_uint_vector(uints)
        if index < self.n_records:
            if KmerUintComparator(uints) == self.graph_kmer_sequence[index]:
                kmer_data = self.graph_sequence[index]._kmer_data
                kmer_data._kmer = kmer_string
                return kmer_data
        raise KeyError('Could not retrieve kmer: ' + kmer_string)

    def __getitem__(self, kmer_string):
        return Kmer(self._cached_get_kmer_data_for_string(kmer_string))

    def __len__(self):
        return max(0, self.n_records)

    def __iter__(self):
        """Iterate over kmers in graph in order stored in graph"""
        self.graph_handle.seek(self.graph_sequence.body_start)
        return kmer_generator_from_stream_and_header(self.graph_handle, self.header)

    def get_kmer_for_string(self, string):
        """Will compute the revcomp of kmer string before getting a kmer"""
        return self[lexlo(string)]

    @property
    def num_colors(self):
        return self.header.num_colors

    @property
    def colors(self):
        return self.header.colors

    @property
    def sample_names(self):
        return self.header.sample_names

    @property
    def kmer_size(self):
        return self.header.kmer_size


@attr.s(slots=True)
class KmerRecordSequence(Sequence):
    graph_handle = attr.ib()
    header = attr.ib()
    body_start = attr.ib()
    n_records = attr.ib()
    kmer_cache_size = attr.ib(0)
    record_size = attr.ib(init=False)
    num_colors = attr.ib(init=False)
    kmer_size = attr.ib(init=False)
    kmer_container_size = attr.ib(init=False)
    _cached_get_kmer_data_for_item = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.record_size = self.header.record_size
        self.kmer_size = self.header.kmer_size
        self.num_colors = self.header.num_colors
        self.kmer_container_size = self.header.kmer_container_size
        self._cached_get_kmer_data_for_item = lru_cache(maxsize=self.kmer_cache_size)(
            self._get_kmer_data_for_item)

    def __getitem__(self, item):
        if not isinstance(item, int):
            raise TypeError("Index must be of type int")
        if item >= self.n_records or item < 0:
            raise IndexError("Index ({}) is out of range".format(item))
        kmer_data = self._cached_get_kmer_data_for_item(item)
        return Kmer(kmer_data)

    def __len__(self):
        return max(0, self.n_records)

    def _get_kmer_data_for_item(self, item):
        self.graph_handle.seek(self.body_start + self.record_size * item)
        kmer_bytes = self.graph_handle.read(self.record_size)
        return KmerData(
            kmer_bytes,
            kmer_size=self.kmer_size,
            num_colors=self.num_colors,
        )


@attr.s(slots=True)
class KmerUintSequence(Sequence):
    graph_handle = attr.ib()
    header = attr.ib()
    body_start = attr.ib()
    n_records = attr.ib()
    kmer_cache_size = attr.ib(0)
    record_size = attr.ib(init=False)
    kmer_container_size = attr.ib(init=False)
    _cached_get_kmer_data_for_item = attr.ib(init=False)
    kmer_string_converter = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.record_size = self.header.record_size
        self.kmer_container_size = self.header.kmer_container_size
        if self.kmer_cache_size == 0:
            self._cached_get_kmer_data_for_item = self._get_kmer_data_for_item
        else:
            self._cached_get_kmer_data_for_item = lru_cache(maxsize=self.kmer_cache_size)(
                self._get_kmer_data_for_item)
        self.kmer_string_converter = StringKmerConverter(self.header.kmer_size)

    def __getitem__(self, item):
        # if not isinstance(item, int):
        #     raise TypeError("Index must be of type int")
        if item >= self.n_records or item < 0:
            raise IndexError("Index ({}) is out of range".format(item))
        kmer_uints = self._cached_get_kmer_data_for_item(item)
        return KmerUintComparator(kmer_uints=kmer_uints)

    def __len__(self):
        return max(0, self.n_records)

    def _get_kmer_data_for_item(self, item):
        self.graph_handle.seek(self.body_start + self.record_size * item)
        kmer_bytes = self.graph_handle.read(self.kmer_container_size * UINT64_T)
        return np.frombuffer(kmer_bytes, dtype='<u8')

    def index_kmer_string(self, kmer_string):
        uints = self.kmer_string_converter.to_uints(kmer_string)
        return self.index_uint_vector(uints)

    def index_uint_vector(self, uints):
        return bisect_left(self, KmerUintComparator(uints))

from bisect import bisect_left
from collections import Sequence, Mapping
from io import SEEK_END

import attr

import cortexpy.graph.parser.header
from cortexpy.graph.parser.kmer import KmerByStringComparator, Kmer, KmerData
from cortexpy.graph.parser.streaming import kmer_generator_from_stream_and_header
from cortexpy.utils import lexlo


class RandomAccessError(KeyError):
    """Raise this if a random access parser could not find a kmer"""


@attr.s(slots=True)
class RandomAccess(Mapping):
    graph_handle = attr.ib()
    header = attr.ib(init=False)
    graph_sequence = attr.ib(init=False)
    n_records = attr.ib(init=False)

    def __attrs_post_init__(self):
        assert self.graph_handle.seekable()
        self.graph_handle.seek(0)
        self.header = cortexpy.graph.parser.header.from_stream(self.graph_handle)
        body_start_stream_position = self.graph_handle.tell()

        self.graph_handle.seek(0, SEEK_END)
        body_size = self.graph_handle.tell() - body_start_stream_position
        if body_size % self.header.record_size != 0:
            raise RuntimeError(
                "Body size ({}) % Record size ({}) != 0".format(body_size, self.header.record_size))
        self.n_records = body_size // self.header.record_size
        self.graph_sequence = KmerRecordSequence(graph_handle=self.graph_handle,
                                                 body_start=body_start_stream_position,
                                                 header=self.header,
                                                 n_records=self.n_records)

    def __getitem__(self, kmer_string):
        kmer = KmerByStringComparator(kmer=kmer_string)
        try:
            kmer_comparator = index_and_retrieve(self.graph_sequence, kmer)
        except ValueError as exception:
            raise KeyError('Could not retrieve kmer: ' + kmer_string) from exception

        return kmer_comparator.kmer_object

    def __len__(self):
        return max(0, self.n_records)

    def __iter__(self):
        self.graph_handle.seek(self.graph_sequence.body_start)
        return kmer_generator_from_stream_and_header(self.graph_handle, self.header)

    def get_kmer_for_string(self, kmer_string):
        """Will compute the revcomp of kmer string before getting a kmer"""
        return self[lexlo(kmer_string)]


# copied from https://docs.python.org/3.6/library/bisect.html
def index_and_retrieve(sequence, value):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(sequence, value)
    if i != len(sequence):
        val = sequence[i]
        if val == value:
            return val
    raise ValueError("Could not find '{}'".format(value))


@attr.s()
class KmerRecordSequence(Sequence):
    graph_handle = attr.ib()
    header = attr.ib()
    body_start = attr.ib()
    n_records = attr.ib()
    record_size = attr.ib(init=False)
    num_colors = attr.ib(init=False)
    kmer_size = attr.ib(init=False)
    kmer_container_size = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.record_size = self.header.record_size
        self.kmer_size = self.header.kmer_size
        self.num_colors = self.header.num_colors
        self.kmer_container_size = self.header.kmer_container_size

    def __getitem__(self, item):
        if not isinstance(item, int):
            raise TypeError("Index must be of type int")
        if item >= self.n_records or item < 0:
            raise IndexError("Index ({}) is out of range".format(item))
        self.graph_handle.seek(self.body_start + self.record_size * item)
        kmer_bytes = self.graph_handle.read(self.record_size)
        return KmerByStringComparator(
            kmer_object=Kmer(KmerData(
                kmer_bytes,
                kmer_size=self.kmer_size,
                num_colors=self.num_colors,
                kmer_container_size_in_uint64ts=self.kmer_container_size,
            ))
        )

    def __len__(self):
        return max(0, self.n_records)

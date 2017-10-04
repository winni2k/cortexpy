from bisect import bisect_left
from collections import Sequence
from io import SEEK_END

import attr

from pycortex.utils import revcomp
from pycortex.graph.parser.header import Header
from pycortex.kmer import KmerByStringComparator, Kmer


class RandomAccessError(KeyError):
    """Raise this if a random access parser could not find a kmer"""


@attr.s(slots=True)
class RandomAccess(object):
    fh = attr.ib()
    header = attr.ib(init=False)
    graph_sequence = attr.ib(init=False)

    def __attrs_post_init__(self):
        assert self.fh.seekable()
        self.fh.seek(0)
        self.header = Header.from_stream(self.fh)
        body_start_stream_position = self.fh.tell()

        self.fh.seek(0, SEEK_END)
        body_size = self.fh.tell() - body_start_stream_position
        if body_size % self.header.record_size != 0:
            raise RuntimeError(
                "Body size ({}) % Record size ({}) != 0".format(body_size, self.header.record_size))
        n_records = body_size // self.header.record_size
        self.graph_sequence = KmerRecordSequence(fh=self.fh,
                                                 body_start=body_start_stream_position,
                                                 cortex_header=self.header,
                                                 n_records=n_records)

    def get_kmer(self, kmer_string):
        kmer = KmerByStringComparator(kmer=kmer_string)
        try:
            kmer_comparator = index(self.graph_sequence, kmer, retrieve=True)
        except ValueError as e:
            raise RandomAccessError('Could not retrieve kmer: ' + kmer_string) from e

        return kmer_comparator.kmer_object

    def get_kmer_for_string(self, kmer_string):
        """Will compute the revcomp of kmer string before getting a kmer"""
        kmer_string_revcomp = revcomp(kmer_string)
        if kmer_string < kmer_string_revcomp:
            return self.get_kmer(kmer_string)
        else:
            return self.get_kmer(kmer_string_revcomp)


# copied from https://docs.python.org/3.6/library/bisect.html
def index(a, x, retrieve=False):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    if i != len(a):
        val = a[i]
        if val == x:
            if retrieve:
                return val
            else:
                return i
    raise ValueError("Could not find '{}'".format(x))


@attr.s()
class KmerRecordSequence(Sequence):
    fh = attr.ib()
    cortex_header = attr.ib()
    body_start = attr.ib()
    n_records = attr.ib()
    record_size = attr.ib(init=False)
    num_colors = attr.ib(init=False)
    kmer_size = attr.ib(init=False)
    kmer_container_size = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.record_size = self.cortex_header.record_size
        self.kmer_size = self.cortex_header.kmer_size
        self.num_colors = self.cortex_header.num_colors
        self.kmer_container_size = self.cortex_header.kmer_container_size

    def __getitem__(self, item):
        if not isinstance(item, int):
            raise TypeError("Index must be of type int")
        if item >= self.n_records or item < 0:
            raise IndexError("Index ({}) is out of range".format(item))
        self.fh.seek(self.body_start + self.record_size * item)
        kmer_bytes = self.fh.read(self.record_size)
        return KmerByStringComparator(
            kmer_object=Kmer(
                kmer_bytes,
                kmer_size=self.kmer_size,
                num_colors=self.num_colors,
                kmer_container_size_in_uint64ts=self.kmer_container_size,
            )
        )

    def __len__(self):
        return self.n_records

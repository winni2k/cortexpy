from itertools import chain

import attr
from collections import Mapping

from cortexpy.utils import lexlo
from .kmer import Kmer, EmptyKmerBuilder
from .kmer_collection import KmerDataCollection


@attr.s(slots=True)
class RandomAccessCollection(Mapping):
    ra_parsers = attr.ib()
    empty_kmer_builders = attr.ib(init=False)
    num_colors = attr.ib(init=False)
    colors = attr.ib(init=False)

    def __attrs_post_init__(self):
        assert all(ra.kmer_size == self.ra_parsers[0].kmer_size for ra in self.ra_parsers)
        self.empty_kmer_builders = [EmptyKmerBuilder(num_colors=ra.num_colors, default_coverage=0)
                                    for ra in
                                    self.ra_parsers]
        self.num_colors = sum(ra.num_colors for ra in self.ra_parsers)
        self.colors = tuple(range(self.num_colors))

    def __getitem__(self, kmer_string):
        kmers = []
        key_errors = [False for _ in range(len(self.ra_parsers))]
        for parser_idx, parser in enumerate(self.ra_parsers):
            try:
                kmer = parser[kmer_string]
            except KeyError:
                key_errors[parser_idx] = True
                kmer = self.empty_kmer_builders[parser_idx].build_or_get(kmer_string)
            kmers.append(kmer)
        if all(key_errors):
            raise KeyError
        return Kmer.from_kmer_data(KmerDataCollection(kmers))

    def __len__(self):
        return max(0, max((parser.n_records for parser in self.ra_parsers)))

    def __iter__(self):
        raise NotImplementedError

    def get_kmer_for_string(self, string):
        """Will compute the revcomp of string before getting a kmer"""
        return self[lexlo(string)]

    @property
    def sample_names(self):
        return chain.from_iterable(ra.sample_names for ra in self.ra_parsers)

    @property
    def kmer_size(self):
        return self.ra_parsers[0].kmer_size

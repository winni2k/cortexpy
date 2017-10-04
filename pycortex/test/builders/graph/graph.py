from io import BytesIO

import attr
from attr import Factory
from pycortex.test.builders.graph.header import CortexGraphHeaderBuilder
from pycortex.test.builders.graph.body import CortexGraphBodyBuilder


@attr.s
class CortexGraphBuilder(object):
    header = attr.ib(Factory(CortexGraphHeaderBuilder))
    body = attr.ib(Factory(CortexGraphBodyBuilder))

    def with_sorted_kmers(self):
        self.body.sort_kmers = True
        return self

    def with_kmer_size(self, size):
        self.body.kmer_size = size
        return self

    def with_kmer_record(self, record):
        self.body.with_kmer_record(record)
        return self

    def with_num_colors(self, n_colors):
        self.header.with_num_colors(n_colors)
        return self

    def build(self):
        self.header.with_kmer_size(self.body.kmer_size)
        self.header.with_kmer_container_size(self.body.kmer_container_size)

        return BytesIO(self.header.build().getvalue() + self.body.build().getvalue())

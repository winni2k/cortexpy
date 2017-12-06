import attr
import numpy as np


@attr.s(slots=True)
class KmerGraphExpectation(object):
    graph = attr.ib()

    def has_node(self, node):
        assert node in self.graph
        return KmerNodeExpectation(self.graph.node[node])

    def has_n_nodes(self, n):
        assert len(self.graph) == n
        return self

    def has_n_edges(self, n):
        assert len(self.graph.edges) == n
        return self

    def has_nodes(self, *nodes):
        assert set(self.graph.nodes) == set(nodes)
        return self

    def has_edges(self, *edges):
        assert set(self.graph.edges) == set(edges)
        return self

    def has_edge(self, source, target, color):
        assert (source, target, color) in self.graph.edges
        return self


@attr.s
class KmerNodeExpectation(object):
    kmer_node = attr.ib()

    def has_coverages(self, *coverages):
        assert np.all(self.kmer_node['kmer'].coverage == np.array(coverages))
        return self

    def is_missing(self):
        print(self.kmer_node['kmer'].coverage)
        assert all(c == 0 for c in self.kmer_node['kmer'].coverage)
        return self

    def is_not_missing(self):
        assert not self.kmer_node.get('is_missing', False)
        return self

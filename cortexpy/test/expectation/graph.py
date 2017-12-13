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

    def has_nodes(self, *expected_nodes):
        expected_nodes = set(expected_nodes)
        assert set(self.graph.nodes) == expected_nodes
        assert len(self.graph.nodes) == len(expected_nodes)
        return self

    def has_edges(self, *edges):
        expected_edges = []
        for edge in edges:
            if isinstance(edge, str):
                edge = edge.split(' ')
                edge[2] = int(edge[2])
            expected_edges.append(tuple(edge))

        assert set(self.graph.edges) == set(expected_edges)
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

import collections
import pprint

import attr

from cortexpy.test.expectation.graph import KmerGraphExpectation


@attr.s
class CollapsedKmerNodeExpectation(object):
    kmer_node = attr.ib()

    def has_coverages_by_kmer(self, *coverages):
        print(self.kmer_node)
        for expected, actual in zip(coverages, self.kmer_node['coverage']):
            assert expected == list(actual)
        return self


@attr.s(slots=True)
class CollapsedKmerUnitgGraphExpectation(object):
    graph = attr.ib()
    repr_counts = attr.ib(attr.Factory(collections.Counter))
    contig_counts = attr.ib(attr.Factory(collections.Counter))
    nodes_by_repr = attr.ib(attr.Factory(dict))
    graph_expectation = attr.ib(init=False)

    def __attrs_post_init__(self):
        pprint.pprint([self.graph.node[n] for n in self.graph.nodes])
        for node, data in self.graph.nodes.data():
            self.repr_counts[data['repr']] += 1
            self.contig_counts[data['unitig']] += 1
        for node, data in self.graph.nodes.data():
            if self.repr_counts[data['repr']] == 1:
                self.nodes_by_repr[data['repr']] = node
        self.graph_expectation = KmerGraphExpectation(self.graph)

    def is_empty(self):
        self.has_n_nodes(0)
        self.has_n_edges(0)
        return self

    def has_n_nodes(self, n):
        self.graph_expectation.has_n_nodes(n)
        return self

    def has_n_missing_kmers(self, n):
        print([self.graph.node[n] for n in self.graph.nodes])
        missing_nodes = [self.graph.node[n] for n in self.graph.nodes if
                         self.graph.node[n].get('is_missing', False)]
        assert len(missing_nodes) == n
        return self

    def has_contigs(self, *contigs):
        assert set(contigs) == set(self.contig_counts.keys())
        return self

    def has_reprs(self, *kmer_reprs):
        assert set(kmer_reprs) == set(self.repr_counts.keys())
        return self

    def has_kmers(self, *kmer_reprs):
        raise NotImplementedError

    def has_n_kmers_with_repr(self, n, kmer_repr):
        assert kmer_repr in self.repr_counts
        assert self.repr_counts[kmer_repr] == n
        return self

    def has_one_node_with_repr(self, kmer_repr):
        self.has_n_kmers_with_repr(1, kmer_repr)
        return CollapsedKmerNodeExpectation(self.graph.node[self.nodes_by_repr[kmer_repr]])

    def has_n_edges(self, n):
        self.graph_expectation.has_n_edges(n)
        return self

    def has_n_missing_edges(self, n):
        missing_edge_counter = 0
        for _, _, data in self.graph.edges.data():
            if data.get('is_missing', False):
                missing_edge_counter += 1
        assert missing_edge_counter == n
        return self

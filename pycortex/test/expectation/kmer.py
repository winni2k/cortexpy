import pprint

import attr
import collections

import numpy as np


@attr.s
class CollapsedKmerNodeExpectation(object):
    kmer_node = attr.ib()

    def has_coverages(self, *coverages):
        assert np.all(self.kmer_node['coverage'] == np.array(coverages, dtype=np.uint32))
        return self

    def is_missing(self):
        raise NotImplementedError

    def is_not_missing(self):
        raise NotImplementedError


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


@attr.s(slots=True)
class CollapsedKmerUnitgGraphExpectation(object):
    graph = attr.ib()
    repr_counts = attr.ib(attr.Factory(collections.Counter))
    nodes_by_repr = attr.ib(attr.Factory(dict))

    def __attrs_post_init__(self):
        pprint.pprint([self.graph.node[n] for n in self.graph.nodes])
        for node, data in self.graph.nodes.data():
            self.repr_counts[data['repr']] += 1
        for node, data in self.graph.nodes.data():
            if self.repr_counts[data['repr']] == 1:
                self.nodes_by_repr[data['repr']] = node

    def has_n_nodes(self, n):
        assert len(self.graph) == n
        return self

    def has_n_missing_kmers(self, n):
        print([self.graph.node[n] for n in self.graph.nodes])
        missing_nodes = [self.graph.node[n] for n in self.graph.nodes if
                         self.graph.node[n].get('is_missing', False)]
        assert len(missing_nodes) == n
        return self

    def has_kmers(self, *kmer_reprs):
        assert set(self.repr_counts.keys()) == set(kmer_reprs)
        return self

    def has_n_kmers_with_repr(self, n, kmer_repr):
        assert kmer_repr in self.repr_counts
        assert self.repr_counts[kmer_repr] == n
        return self

    def has_one_node_with_repr(self, kmer_repr):
        self.has_n_kmers_with_repr(1, kmer_repr)
        return CollapsedKmerNodeExpectation(self.graph.node[self.nodes_by_repr[kmer_repr]])

    def has_n_edges(self, n):
        assert len(self.graph.edges) == n
        return self

    def has_n_missing_edges(self, n):
        missing_edge_counter = 0
        for _, _, data in self.graph.edges.data():
            if data.get('is_missing', False):
                missing_edge_counter += 1
        assert missing_edge_counter == n
        return self

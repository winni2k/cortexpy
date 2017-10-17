import attr


@attr.s
class KmerNodeExpectation(object):
    kmer_node = attr.ib()

    def has_coverages(self, *coverages):
        assert self.kmer_node['coverage'] == list(coverages)
        return self

    def is_missing(self):
        assert self.kmer_node['is_missing']
        return self

    def is_not_missing(self):
        assert not self.kmer_node['is_missing']
        return self


@attr.s(slots=True)
class CollapsedKmerUnitgGraphExpectation(object):
    graph = attr.ib()
    nodes_by_repr = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.nodes_by_repr = {data['repr']: node for node, data in self.graph.nodes.data()}

    def has_n_kmers(self, n):
        assert len(self.graph) == n
        return self

    def has_kmers(self, *kmer_reprs):
        assert set(self.nodes_by_repr.keys()) == set(kmer_reprs)
        return self

    def has_kmer(self, kmer_repr):
        assert kmer_repr in self.nodes_by_repr
        return KmerNodeExpectation(self.graph.node[self.nodes_by_repr[kmer_repr]])

    def has_n_edges(self, n):
        assert len(self.graph.edges) == n
        return self

    def has_n_missing_edges(self, n):
        missing_edge_counter = 0
        for _, _, data in self.graph.edges.data():
            if data['is_missing']:
                missing_edge_counter += 1
        assert missing_edge_counter == n
        return self

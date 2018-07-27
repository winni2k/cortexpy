import attr
import numpy as np


@attr.s(slots=True)
class UnitigExpectation(object):
    unitig_data = attr.ib()

    def has_contig(self, contig):
        assert contig == self.unitig_data.contig
        return self

    def with_left_node(self, node):
        assert node == self.unitig_data.left_node
        return self

    def with_right_node(self, node):
        assert node == self.unitig_data.right_node
        return self

    def with_coverage(self, coverage_matrix):
        assert np.array_equal(coverage_matrix, self.unitig_data['coverage'])
        return self

    # fixme: all nodes should have is_missing set
    def is_not_missing(self):
        assert not self.unitig_data.get('is_missing')
        return self

    def is_cycle(self):
        assert self.unitig_data.is_cycle
        return self

    def is_not_cycle(self):
        assert not self.unitig_data.is_cycle
        return self


@attr.s(slots=True)
class GraphWithUnitigExpectation(object):
    graph = attr.ib()
    unitigs = attr.ib(attr.Factory(list))
    by_contigs = attr.ib(attr.Factory(dict))

    def __attrs_post_init__(self):
        for source, target, unitig in self.graph.edges(data='unitig'):
            if unitig:
                self.unitigs.append((source, target, unitig))
                self.by_contigs[unitig.contig] = unitig

    def has_n_nodes(self, n):
        try:
            self.has_n_unitigs(n)
        except Exception:
            print(list(self.graph))
            raise
        return self

    def has_n_edges(self, n):
        assert n == len(self.graph.edges)
        return self

    def has_n_unitigs(self, n):
        assert n == len(self.unitigs)
        return self

    def has_one_unitig(self):
        return self.has_n_unitigs(1)

    def is_equivalent_to(self, graph):
        assert set(self.graph.edges) == set(graph.edges)
        return self

    def has_unitig_with_edges(self, *expected_edges):
        expected_edge_set = set(expected_edges)
        actual_edge_sets = [set(u.graph.edges()) for _, _, u in self.unitigs]
        assert expected_edge_set in actual_edge_sets
        return UnitigExpectation(self.unitigs[actual_edge_sets.index(expected_edge_set)][2])

    def has_unitig_with_one_node(self, expected_node):
        unitigs_with_node = list(filter(lambda x: x[2].left_node == expected_node, self.unitigs))
        assert len(unitigs_with_node) == 1
        return UnitigExpectation(unitigs_with_node[0][2])

    def has_unitig_with_contig(self, contig):
        assert contig in self.by_contigs
        return UnitigExpectation(self.by_contigs[contig])

import attr
import numpy as np


@attr.s(slots=True)
class UnitigExpectation(object):
    unitig = attr.ib()
    unitig_data = attr.ib()

    def with_left_node(self, node):
        assert self.unitig_data['left_node'] == node
        return self

    def with_right_node(self, node):
        assert self.unitig_data['right_node'] == node
        return self

    def is_missing(self):
        assert self.unitig_data['is_missing']
        return self

    def with_coverage(self, coverage_matrix):
        assert np.array_equal(self.unitig_data['coverage'], coverage_matrix)
        return self

    # fixme: all nodes should have is_missing set
    def is_not_missing(self):
        assert not self.unitig_data.get('is_missing')
        return self


@attr.s(slots=True)
class GraphWithUnitigExpectation(object):
    graph = attr.ib()
    unitigs = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.unitigs = []
        for node, data in self.graph.nodes.data():
            if data.get('is_unitig', False):
                self.unitigs.append((node, data))

    def has_n_nodes(self, n):
        try:
            self.has_n_unitigs(n)
        except Exception:
            print(list(self.graph))
            raise
        return self

    def has_n_edges(self, n):
        assert len(self.graph.edges) == n
        return self

    def has_n_unitigs(self, n):
        assert len(self.unitigs) == n
        return self

    def has_one_unitig(self):
        return self.has_n_unitigs(1)

    def is_equivalent_to(self, graph):
        assert set(self.graph.edges) == set(graph.edges)
        return self

    def has_unitig_with_edges(self, *expected_edges):
        expected_edge_set = set(expected_edges)
        actual_edge_sets = [set(g.edges()) for g, _ in self.unitigs]
        assert expected_edge_set in actual_edge_sets
        return UnitigExpectation(*self.unitigs[actual_edge_sets.index(expected_edge_set)])

    def has_unitig_with_one_node(self, expected_node):
        unitigs_with_node = list(filter(lambda x: x[1]['left_node'] == expected_node, self.unitigs))
        assert len(unitigs_with_node) == 1
        return UnitigExpectation(unitigs_with_node[0][0], unitigs_with_node[0][1])

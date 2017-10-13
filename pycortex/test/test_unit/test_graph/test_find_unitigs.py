import attr
import networkx as nx

from pycortex.graph.serializer import find_unitigs


@attr.s(slots=True)
class UnitigGraphExpectation(object):
    unitig = attr.ib()
    unitig_data = attr.ib()

    def with_left_node(self, node):
        assert self.unitig_data['left_node'] == node
        return self

    def with_right_node(self, node):
        assert self.unitig_data['right_node'] == node
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
        assert len(self.graph) == n
        return self

    def has_n_unitigs(self, n):
        assert len(self.unitigs) == n
        return self

    def has_one_unitig(self):
        return self.has_n_unitigs(1)

    def is_equivalent_to(self, graph):
        assert set(self.graph.edges) == set(graph.edges)
        return self

    def has_unitig_with_edges(self, *expected_edge_set):
        expected_edge_set = set(expected_edge_set)
        actual_edge_sets = [set(g.edges) for g, _ in self.unitigs]
        assert expected_edge_set in actual_edge_sets
        return UnitigGraphExpectation(*self.unitigs[actual_edge_sets.index(expected_edge_set)])


class TestFindUnitigs(object):
    def test_three_node_path_becomes_a_unitig(self):
        graph = nx.DiGraph()
        graph.add_path(range(3))
        expect = GraphWithUnitigExpectation(find_unitigs(graph))

        expect.has_n_nodes(1)
        expect.has_unitig_with_edges(*graph.edges)

    def test_three_node_path_with_mixed_node_order_becomes_a_unitig(self):
        # given
        graph = nx.DiGraph()
        graph.add_edge(0, 2)
        graph.add_edge(1, 0)

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(graph))

        # then
        expect.has_unitig_with_edges(*graph.edges)

    def test_two_node_cycle_becomes_unitig(self):
        # given
        graph = nx.DiGraph()
        graph.add_cycle(range(2))

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(graph))

        # then
        expect.has_unitig_with_edges((0, 1))  # (1,0) would also be appropriate

    def test_two_node_path_becomes_unitig(self):
        # given
        graph = nx.DiGraph()
        graph.add_path(range(2))

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(graph))

        # then
        expect.has_unitig_with_edges((0, 1))

    def test_three_node_cycle_becomes_three_node_unitig(self):
        # given
        graph = nx.DiGraph()
        graph.add_cycle(range(3))

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(graph))

        # then
        (expect
         .has_n_nodes(1)
         .has_one_unitig()
         .has_unitig_with_edges((0, 1), (1, 2))  # Other unitigs are appropriate as well
         .with_left_node(0)
         .with_right_node(2))

    def test_four_node_cycle_becomes_four_node_unitig(self):
        # given
        graph = nx.DiGraph()
        graph.add_cycle(range(4))

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(graph))

        # then
        (expect
         .has_n_nodes(1)
         .has_one_unitig()
         .has_unitig_with_edges((0, 1), (1, 2), (2, 3))
         .with_left_node(0)
         .with_right_node(3))

    def test_path_and_cycle_remains_unchanged(self):
        # given
        graph = nx.DiGraph()
        graph.add_path(range(4))
        graph.add_cycle([1, 2, 4])

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(graph))

        # then
        expect.is_equivalent_to(graph)

    def test_two_node_path_and_three_node_cycle_becomes_two_unitigs(self):
        # given
        graph = nx.DiGraph()
        graph.add_path(range(3))
        graph.add_cycle(range(2, 5))

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(graph))

        # then
        (expect
         .has_n_unitigs(2)
         .has_n_nodes(2))
        (expect
         .has_unitig_with_edges((0, 1))
         .with_left_node(0)
         .with_right_node(1))
        (expect
         .has_unitig_with_edges((2, 3), (3, 4))
         .with_left_node(2)
         .with_right_node(4))

    def test_cycle_and_six_node_path_results_in_one_unitig(self):
        # given
        graph = nx.DiGraph()
        graph.add_cycle([1, 2, 3])
        graph.add_path([4, 6, 0, 1, 2, 5])

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(graph))

        # then
        (expect
         .has_one_unitig()
         .has_unitig_with_edges((4, 6), (6, 0))
         .with_left_node(4)
         .with_right_node(0))

    def test_two_paths_making_bubble_results_in_two_unitigs(self):
        # given
        graph = nx.DiGraph()
        graph.add_path(range(4))
        graph.add_path([1, 4, 2])

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(graph))

        # then
        expect.has_n_nodes(3)
        expect.has_n_unitigs(2)
        expect.has_unitig_with_edges((0, 1)).with_left_node(0).with_right_node(1)
        expect.has_unitig_with_edges((2, 3)).with_left_node(2).with_right_node(3)

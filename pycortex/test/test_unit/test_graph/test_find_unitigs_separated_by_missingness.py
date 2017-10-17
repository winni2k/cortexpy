import networkx as nx

from pycortex.graph.serializer import find_unitigs
from pycortex.test.expectation.unitig_graph import GraphWithUnitigExpectation


class TestWithMissingNodes(object):
    def test_path_two_and_path_two_with_missing_nodes_results_in_two_unitigs(self):
        # given
        graph = nx.DiGraph()
        graph.add_path(range(2))
        graph.add_path(range(2, 4), is_missing=True)

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(graph))

        # then
        expect.has_n_nodes(2)
        expect.has_n_unitigs(2)
        expect.has_unitig_with_edges((0, 1)).with_left_node(0).with_right_node(1).is_not_missing()
        expect.has_unitig_with_edges((2, 3)).with_left_node(2).with_right_node(3).is_missing()

import networkx as nx

from pycortex.graph.serializer import find_unitigs
from pycortex.test.builder.graph.networkx import NetworkxGraphBuilder
from pycortex.test.expectation.unitig_graph import GraphWithUnitigExpectation


class TestWithMissingEdge(object):
    def test_two_edges_and_two_edges_results_in_two_unitigs(self):
        # given
        builder = NetworkxGraphBuilder()
        builder.graph.add_edge(0, 1)
        builder.graph.add_edge(2, 3)

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(builder.build()))

        # then
        expect.has_n_nodes(2)
        expect.has_n_unitigs(2)
        expect.has_unitig_with_edges((0, 1)).with_left_node(0).with_right_node(1)
        expect.has_unitig_with_edges((2, 3)).with_left_node(2).with_right_node(3)

    def test_two_pairs_of_edges_separated_by_single_edge_results_in_two_unitigs(self):
        # given
        colors = [0, 1]
        builder = NetworkxGraphBuilder()
        builder.with_colors(*colors)
        nx.add_path(builder.graph, range(4), key=0)
        builder.add_edge_with_color(0, 1, 1)
        builder.add_edge_with_color(2, 3, 1)

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(builder.build(), colors=colors))

        # then
        expect.has_n_nodes(2)
        expect.has_n_unitigs(2)
        expect.has_unitig_with_edges((0, 1)).with_left_node(0).with_right_node(1)
        expect.has_unitig_with_edges((2, 3)).with_left_node(2).with_right_node(3)

    def test_two_pairs_of_nodes_separated_by_single_edge_results_in_two_unitigs(self):
        # given
        builder = NetworkxGraphBuilder()
        colors = (0, 1)
        builder.with_colors(*colors)
        builder.add_edge_with_color(0, 1, 0)

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(builder.build(), colors=colors))

        # then
        expect.has_n_unitigs(2)
        expect.has_n_edges(1)
        expect.has_unitig_with_one_node(0)
        expect.has_unitig_with_one_node(1)

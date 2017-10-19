import networkx as nx
import pytest

from pycortex.graph.serializer import find_unitigs
from pycortex.test.builder.graph.networkx import NetworkxGraphBuilder
from pycortex.test.expectation.unitig_graph import GraphWithUnitigExpectation


class TestWithMissingNodes(object):
    def test_two_results_in_one_unitig(self):
        # given
        builder = NetworkxGraphBuilder()
        graph = builder.graph
        graph.add_path(range(2), is_missing=True)
        graph.node[0]['is_missing'] = True
        graph.node[1]['is_missing'] = True
        graph = builder.build()

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(graph))

        # then
        expect.has_n_nodes(1)
        expect.has_n_unitigs(1)
        (expect.has_unitig_with_edges((0, 1))
         .with_left_node(0)
         .with_right_node(1)
         .with_coverage([(1,), (1,)]))

    @pytest.mark.xfail(reason="Requires the use of MultiDiGraph under the hood")
    def test_path_two_and_path_two_results_in_two_unitigs(self):
        # given
        graph = nx.DiGraph()
        graph.add_path(range(2))
        graph.add_path(range(2, 4))
        graph.node[2]['is_missing'] = True
        graph.node[3]['is_missing'] = True

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(graph))

        # then
        expect.has_n_nodes(2)
        expect.has_n_unitigs(2)
        expect.has_unitig_with_edges((0, 1)).with_left_node(0).with_right_node(1).is_not_missing()
        expect.has_unitig_with_edges((2, 3)).with_left_node(2).with_right_node(3).is_missing()

    def test_path_two_and_path_three_results_in_two_unitigs(self):
        # given
        builder = NetworkxGraphBuilder()
        graph = builder.graph
        graph.add_edge(0, 1)
        graph.add_path(range(1, 4), is_missing=True)

        (builder
         .with_node_coverage(0, 1)
         .with_node_coverage(1, 1)
         .with_node_coverage(2, 0)
         .with_node_coverage(3, 0))
        graph = builder.build()

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(graph))

        # then
        expect.has_n_nodes(2)
        expect.has_n_unitigs(2)
        expect.has_unitig_with_edges((0, 1)).with_left_node(0).with_right_node(1).is_not_missing()
        expect.has_unitig_with_edges((2, 3)).with_left_node(2).with_right_node(3).is_missing()


class TestWithMissingEdge(object):
    @pytest.mark.xfail(reason="Requires the use of MultiDiGraph under the hood")
    def test_path_two_and_path_two_results_in_two_unitigs(self):
        # given
        graph = nx.DiGraph()
        graph.add_edge(0, 1)
        graph.add_edge(1, 2, is_missing=True)
        graph.add_edge(2, 3)

        # when
        expect = GraphWithUnitigExpectation(find_unitigs(graph))

        # then
        expect.has_n_nodes(2)
        expect.has_n_unitigs(2)
        expect.has_unitig_with_edges((0, 1)).with_left_node(0).with_right_node(1).is_not_missing()
        expect.has_unitig_with_edges((2, 3)).with_left_node(2).with_right_node(3).is_not_missing()

import networkx as nx

from cortexpy.test.driver.graph.find_unitigs import FindUnitigsTestDriver


class TestWithMissingEdge(object):
    def test_two_edges_and_two_edges_results_in_two_unitigs(self):
        # given
        driver = FindUnitigsTestDriver()
        driver.without_test_coverage()
        builder = driver.builder
        builder.graph.add_edge(0, 1)
        builder.graph.add_edge(2, 3)

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(2)
        expect.has_n_unitigs(2)
        expect.has_unitig_with_edges((0, 1)).with_left_node(0).with_right_node(1)
        expect.has_unitig_with_edges((2, 3)).with_left_node(2).with_right_node(3)

    def test_two_pairs_of_edges_separated_by_single_edge_results_in_two_unitigs(self):
        # given
        driver = FindUnitigsTestDriver()
        driver.without_test_coverage()
        builder = driver.graph_builder
        builder.with_colors(0, 1)
        nx.add_path(builder.graph, range(4), key=0)
        builder.add_edge_with_color(0, 1, 1)
        builder.add_edge_with_color(2, 3, 1)

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(2)
        expect.has_n_unitigs(2)
        expect.has_unitig_with_edges((0, 1)).with_left_node(0).with_right_node(1)
        expect.has_unitig_with_edges((2, 3)).with_left_node(2).with_right_node(3)

    def test_two_pairs_of_nodes_separated_by_single_edge_results_in_two_unitigs(self):
        # given
        driver = FindUnitigsTestDriver().without_test_coverage()
        builder = driver.graph_builder
        builder.with_colors(0, 1)
        builder.add_edge_with_color(0, 1, 0)

        # when
        expect = driver.run()

        # then
        expect.has_n_unitigs(2)
        expect.has_n_edges(1)
        expect.has_unitig_with_one_node(0)
        expect.has_unitig_with_one_node(1)

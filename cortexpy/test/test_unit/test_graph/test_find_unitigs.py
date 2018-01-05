from unittest.mock import Mock

import networkx as nx
from hypothesis import given, settings
from hypothesis import strategies as s

from cortexpy.constants import EdgeTraversalOrientation
from cortexpy.test.builder.graph.networkx import add_kmers_to_graph
from cortexpy.test.driver.graph.find_unitigs import FindUnitigsTestDriver


class Test(object):
    def test_three_node_path_becomes_a_unitig(self):
        # given
        driver = FindUnitigsTestDriver()
        driver.graph.add_path(range(3))

        # when
        expect = driver.run()

        expect.has_n_nodes(1)
        expect.has_unitig_with_edges((0, 1), (1, 2))

    def test_three_node_path_with_mixed_node_order_becomes_a_unitig(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_edge(0, 2)
        graph.add_edge(1, 0)

        # when
        expect = driver.run()

        # then
        expect.has_unitig_with_edges((0, 2), (1, 0))

    def test_two_node_cycle_becomes_unitig(self):
        # given
        driver = FindUnitigsTestDriver()
        driver.graph.add_cycle(range(2))

        # when
        expect = driver.run()

        # then
        expect.has_n_unitigs(1)
        expect.has_unitig_with_edges((0, 1))  # (1,0) would also be appropriate

    def test_two_node_path_becomes_unitig(self):
        # given
        driver = FindUnitigsTestDriver()
        driver.graph.add_path(range(2))

        # when
        expect = driver.run()

        # then
        expect.has_unitig_with_edges((0, 1))

    def test_three_node_cycle_becomes_three_node_unitig(self):
        # given
        driver = FindUnitigsTestDriver()
        driver.graph.add_cycle(range(3))

        # when
        expect = driver.run()

        # then
        (expect
         .has_n_nodes(1)
         .has_one_unitig()
         .has_unitig_with_edges((0, 1), (1, 2))  # Other unitigs are appropriate as well
         .with_left_node(0)
         .with_right_node(2))

    def test_four_node_cycle_becomes_four_node_unitig(self):
        # given
        driver = FindUnitigsTestDriver()
        driver.graph.add_cycle(range(4))

        # when
        expect = driver.run()

        # then
        (expect
         .has_n_nodes(1)
         .has_one_unitig()
         .has_unitig_with_edges((0, 1), (1, 2), (2, 3))
         .with_left_node(0)
         .with_right_node(3))

    def test_path_and_cycle_becomes_four_unitigs(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_path(range(4))
        graph.add_edge(2, 4)
        graph.add_edge(4, 1)

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(4)
        expect.has_n_unitigs(4)
        expect.has_unitig_with_edges((1, 2))
        expect.has_unitig_with_one_node(0)
        expect.has_unitig_with_one_node(3)
        expect.has_unitig_with_one_node(4)

    def test_two_node_path_and_three_node_cycle_becomes_two_unitigs(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_path(range(3))
        graph.add_cycle(range(2, 5))

        # when
        expect = driver.run()

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

    def test_cycle_and_six_node_path_results_in_four_unitigs(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_edge(2, 3)
        graph.add_edge(3, 1)
        graph.add_path([4, 6, 0, 1, 2, 5])
        graph = add_kmers_to_graph(graph)

        # when
        expect = driver.run()

        # then
        (expect
         .has_n_unitigs(4)
         .has_unitig_with_edges((4, 6), (6, 0))
         .with_left_node(4)
         .with_right_node(0))
        expect.has_unitig_with_edges((1, 2))
        expect.has_unitig_with_one_node(3)
        expect.has_unitig_with_one_node(5)

    def test_two_paths_making_bubble_results_in_four_unitigs(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_path(range(4))
        graph.add_path([1, 4, 2])

        # when
        expect = driver.run()

        # then
        expect.has_n_unitigs(3)
        expect.has_unitig_with_edges((0, 1)).with_left_node(0).with_right_node(1)
        expect.has_unitig_with_edges((2, 3)).with_left_node(2).with_right_node(3)
        expect.has_unitig_with_one_node(4)


class TestIsUnitigEnd(object):
    def test_single_node_is_end_from_both_sides(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_node(0)
        driver.build()

        # when/then
        for orientation in EdgeTraversalOrientation:
            assert driver.finder.is_unitig_end(0, orientation)

    @given(s.integers(min_value=1, max_value=3))
    @settings(max_examples=3)
    def test_each_end_of_path_is_end(self, num_edges):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        for color in range(num_edges):
            graph.add_edge(0, 1, key=color)
        finder = driver.build()

        # when/then
        assert finder.is_unitig_end(0, EdgeTraversalOrientation.reverse)
        assert not finder.is_unitig_end(0, EdgeTraversalOrientation.original)
        assert not finder.is_unitig_end(1, EdgeTraversalOrientation.reverse)
        assert finder.is_unitig_end(1, EdgeTraversalOrientation.original)

    def test_two_edges_into_one_node(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_edge(0, 1)
        graph.add_edge(2, 1)
        finder = driver.build()

        # when/then
        for node in range(3):
            for orientation in EdgeTraversalOrientation:
                assert finder.is_unitig_end(node, orientation)

    def test_two_edges_out_of_one_node(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_edge(0, 1)
        graph.add_edge(0, 2)
        finder = driver.build()

        # when/then
        for node in range(3):
            for orientation in EdgeTraversalOrientation:
                assert finder.is_unitig_end(node, orientation)

    def test_middle_unconnected_node_in_two_color_graph(self):
        # given
        driver = FindUnitigsTestDriver().with_colors(0, 1)
        graph = driver.graph
        nx.add_path(graph, range(3), key=0)
        finder = driver.build()

        # when/then
        for node in range(3):
            for orientation in EdgeTraversalOrientation:
                assert finder.is_unitig_end(node, orientation)


class TestFindUnitigFromTwoColorGraph(object):
    def test_three_node_path_becomes_a_unitig_and_attributes_are_copied_across(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_path(range(3))
        for node in graph:
            kmer_mock = Mock()
            kmer_mock.coverage = (node + 1, 1)
            graph.node[node]['kmer'] = kmer_mock
            graph.node[node]['bla'] = node * 3
        finder = driver.builder.build()

        # when
        for start_node in range(3):
            unitig = finder.find_unitig_from(start_node)

            # then
            assert unitig.left_node == 0
            assert unitig.right_node == 2
            assert len(unitig.graph) == 3
            assert set(unitig.graph.edges(keys=False)) == {(0, 1), (1, 2)}
            for node in range(3):
                assert unitig.graph.node[node]['kmer'].coverage == (node + 1, 1)
                assert unitig.graph.node[node]['bla'] == node * 3


class TestUnitigGraphCoverage(object):
    @given(s.booleans(), s.booleans(),
           s.lists(s.tuples(s.integers(min_value=0, max_value=1),
                            s.integers(min_value=0, max_value=1)),
                   min_size=2,
                   max_size=2))
    def test_two_node_path_with_missing_links_are_not_joined(self, link_color_0_exists,
                                                             link_color_1_exists,
                                                             kmer_coverages):
        # given
        driver = FindUnitigsTestDriver()
        driver.graph.add_nodes_from([0, 1])
        colors = [0, 1]
        driver.with_colors(*colors)

        for color, link_exists in zip(colors, [link_color_0_exists, link_color_1_exists]):
            if link_exists:
                driver.graph.add_edge(0, 1, key=color)
        driver.build()

        kmers = [driver.graph.node[r]['kmer'] for r in range(2)]
        kmers[0].coverage = kmer_coverages[0]
        kmers[1].coverage = kmer_coverages[1]

        for start_node in range(2):
            # when
            unitig = driver.finder.find_unitig_from(start_node)

            # then
            if ((link_color_0_exists and link_color_1_exists) or
                    (link_color_1_exists and kmer_coverages[0][0] == kmer_coverages[1][0] == 0) or
                    (link_color_0_exists and kmer_coverages[0][1] == kmer_coverages[1][1] == 0)):
                assert len(unitig.coverage) == 2
                assert unitig.coverage[0] == kmer_coverages[0]
                assert unitig.coverage[1] == kmer_coverages[1]
            else:
                assert len(unitig.coverage) == 1
                assert unitig.coverage[0] == kmer_coverages[start_node]

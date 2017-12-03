from unittest.mock import Mock

import networkx as nx
from hypothesis import given, assume, settings
from hypothesis import strategies as s

from pycortex.graph.serializer import find_unitig_from, is_unitig_end, \
    EdgeTraversalOrientation
from pycortex.test.builder.graph.networkx import add_kmers_to_graph
from pycortex.test.driver.graph.find_unitigs import FindUnitigsTestDriver


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
        graph = FindUnitigsTestDriver().graph
        graph.add_node(0)

        # when/then
        for orientation in EdgeTraversalOrientation:
            assert is_unitig_end(0, graph, orientation)

    @given(s.integers(min_value=1, max_value=3))
    @settings(max_examples=3)
    def test_each_end_of_path_is_end(self, num_edges):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        for color in range(num_edges):
            graph.add_edge(0, 1, key=color)
        driver.run()

        # when/then
        assert is_unitig_end(0, graph, EdgeTraversalOrientation.reverse)
        assert not is_unitig_end(0, graph, EdgeTraversalOrientation.original)
        assert not is_unitig_end(1, graph, EdgeTraversalOrientation.reverse)
        assert is_unitig_end(1, graph, EdgeTraversalOrientation.original)

    def test_two_edges_into_one_node(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_edge(0, 1)
        graph.add_edge(2, 1)
        driver.run()

        # when/then
        for node in range(3):
            for orientation in EdgeTraversalOrientation:
                assert is_unitig_end(node, graph, orientation)

    def test_two_edges_out_of_one_node(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_edge(0, 1)
        graph.add_edge(0, 2)
        driver.run()

        # when/then
        for node in range(3):
            for orientation in EdgeTraversalOrientation:
                assert is_unitig_end(node, graph, orientation)

    def test_middle_unconnected_node_in_two_color_graph(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        nx.add_path(graph, range(3), key=0)
        driver.run()

        # when/then
        for node in range(3):
            for orientation in EdgeTraversalOrientation:
                assert is_unitig_end(node, graph, orientation, colors=[0, 1])


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

        # when
        for start_node in range(3):
            unitig = find_unitig_from(start_node, graph)

            # then
            assert unitig.left_node == 0
            assert unitig.right_node == 2
            assert len(unitig.graph) == 3
            assert set(unitig.graph.edges(keys=False)) == {(0, 1), (1, 2)}
            for node in range(3):
                assert unitig.graph.node[node]['kmer'].coverage == (node + 1, 1)
                assert unitig.graph.node[node]['bla'] == node * 3

    @given(s.sampled_from(((0, 1), (1, 0), (1, 1), (0, 0), (0,), (1,), tuple())),
           s.sampled_from(((0, 1), (1, 0), (1, 1), (0, 0), (0,), (1,), tuple())))
    def test_two_node_path_with_differing_missing_kmers_are_joined(self, coverage0, coverage1):
        # given
        assume(len(coverage0) == len(coverage1))
        driver = FindUnitigsTestDriver()
        driver.graph.add_edge(0, 1)

        # when
        driver.run()

        for node, coverage in zip(range(2), [coverage0, coverage1]):
            driver.graph.node[node]['kmer'].coverage = coverage

        # then
        for start_node in range(2):
            unitig = find_unitig_from(start_node, driver.graph)
            assert unitig.left_node == 0
            assert unitig.right_node == 1
            assert len(unitig.graph) == 2
            assert len(unitig.graph.edges) == 1


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
        for color, link_exists in zip(colors, [link_color_0_exists, link_color_1_exists]):
            if link_exists:
                driver.graph.add_edge(0, 1, key=color)
        driver.run()

        kmers = [driver.graph.node[r]['kmer'] for r in range(2)]
        kmers[0].coverage = kmer_coverages[0]
        kmers[1].coverage = kmer_coverages[1]

        for start_node in range(2):
            # when
            unitig = find_unitig_from(start_node, driver.graph, colors=colors, test_coverage=True)

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

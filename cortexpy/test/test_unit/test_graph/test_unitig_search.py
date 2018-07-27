from hypothesis import given, strategies as s

from cortexpy.constants import EdgeTraversalOrientation
from cortexpy.graph.serializer.unitig import UnitigSearch
from cortexpy.test.builder.graph.cortex import CortexGraphBuilder


class TestIsUnitigEnd(object):
    @given(s.integers(min_value=1, max_value=2))
    def test_two_connected_nodes_are_a_unitig(self, num_colors):
        # given
        colors = list(range(num_colors))
        coverage = [1] + [0 for _ in range(1, num_colors)]

        b = CortexGraphBuilder()
        b.with_colors(*colors)
        b.with_node_coverage('AAA', *coverage)
        b.with_node_coverage('AAT', *coverage)
        b.add_edge_with_color('AAA', 'AAT', 0)

        # when
        graph = b.build()

        # then
        search = UnitigSearch.from_node_and_graph('AAA', graph)
        assert not search.is_unitig_end('AAA', EdgeTraversalOrientation.original)
        assert search.is_unitig_end('AAA', EdgeTraversalOrientation.reverse)

        search = UnitigSearch.from_node_and_graph('AAT', graph)
        assert search.is_unitig_end('AAT', EdgeTraversalOrientation.original)
        assert not search.is_unitig_end('AAT', EdgeTraversalOrientation.reverse)

    @given(s.sampled_from(('AAA', 'AAT')),
           s.sampled_from(EdgeTraversalOrientation))
    def test_two_colors_to_one_color_is_node_end(self, start_node, orientation):
        # given
        b = CortexGraphBuilder()
        b.with_colors(0, 1)
        b.with_node_coverage('AAA', 1, 1)
        b.with_node_coverage('AAT', 1, 0)
        b.add_edge_with_color('AAA', 'AAT', 0)

        # when
        graph = b.build()

        # then
        search = UnitigSearch.from_node_and_graph(start_node, graph)
        assert search.is_unitig_end(start_node, orientation)

    @given(s.sampled_from(('AAA', 'AAT', 'AAC')),
           s.sampled_from(EdgeTraversalOrientation))
    def test_one_color_to_two_nodes_is_node_end(self, start_node, orientation):
        # given
        b = CortexGraphBuilder()
        b.with_colors(0)
        b.with_node_coverage('AAA', 1)
        b.with_node_coverage('AAT', 1)
        b.with_node_coverage('AAC', 1)
        b.add_edge_with_color('AAA', 'AAT', 0)
        b.add_edge_with_color('AAA', 'AAC', 0)

        # when
        graph = b.build()

        # then
        search = UnitigSearch.from_node_and_graph(start_node, graph)
        assert search.is_unitig_end(start_node, orientation)

    def test_two_colors_to_two_nodes_in_both_colors_is_node_end(self):
        # given
        b = CortexGraphBuilder()
        b.with_colors(0, 1)
        b.with_node_coverage('AAA', 1, 1)
        b.with_node_coverage('AAT', 1, 1)
        b.with_node_coverage('AAC', 1, 1)
        for color in (0, 1):
            b.add_edge_with_color('AAA', 'AAT', color)
            b.add_edge_with_color('AAA', 'AAC', color)

        # when
        graph = b.build()

        # then
        for orientation in EdgeTraversalOrientation:
            for node in ['AAA', 'AAT', 'AAC']:
                search = UnitigSearch.from_node_and_graph(node, graph)
                assert search.is_unitig_end(node, orientation)

    def test_two_colors_to_two_nodes_in_first_color_is_node_end(self):
        # given
        b = CortexGraphBuilder()
        b.with_colors(0, 1)
        b.with_node_coverage('AAA', 1, 1)
        b.with_node_coverage('AAT', 1, 0)
        b.with_node_coverage('AAC', 1, 0)
        b.add_edge_with_color('AAA', 'AAT', 0)
        b.add_edge_with_color('AAA', 'AAC', 0)

        # when
        graph = b.build()

        # then
        for orientation in EdgeTraversalOrientation:
            for node in ['AAA', 'AAT', 'AAC']:
                search = UnitigSearch.from_node_and_graph(node, graph)
                assert search.is_unitig_end(node, orientation)

    def test_two_colors_to_two_nodes_in_each_color_is_node_end(self):
        # given
        b = CortexGraphBuilder()
        b.with_colors(0, 1)
        b.with_node_coverage('AAA', 1, 1)
        b.with_node_coverage('AAT', 1, 0)
        b.with_node_coverage('AAC', 0, 1)
        b.add_edge_with_color('AAA', 'AAT', 0)
        b.add_edge_with_color('AAA', 'AAC', 1)

        # when
        graph = b.build()

        # then
        for orientation in EdgeTraversalOrientation:
            for node in ['AAA', 'AAT', 'AAC']:
                search = UnitigSearch.from_node_and_graph(node, graph)
                assert search.is_unitig_end(node, orientation)

    @given(s.sampled_from(('AAA', 'AAT', 'ATA')),
           s.sampled_from(EdgeTraversalOrientation))
    def test_middle_unconnected_node_for_one_color_in_two_color_graph_results_in_three_unitigs(
        self, start_node, orientation
    ):
        # given
        b = CortexGraphBuilder()
        b.with_colors(0, 1)
        b.with_node_coverage('AAA', 1, 1)
        b.with_node_coverage('AAT', 1, 1)
        b.with_node_coverage('ATA', 1, 1)
        b.add_edge_with_color('AAA', 'AAT', 0)
        b.add_edge_with_color('AAT', 'ATA', 0)

        # when
        graph = b.build()
        search = UnitigSearch.from_node_and_graph(start_node, graph)
        assert search.is_unitig_end(start_node, orientation)

    @given(s.data(),
           s.integers(min_value=1, max_value=3),
           s.sampled_from(EdgeTraversalOrientation))
    def test_multi_color_y_graph_always_has_three_unitigs(self, data, num_colors, orientation):
        # given
        colors = list(range(num_colors))
        edge_colors = [data.draw(s.sampled_from(colors)) for _ in range(3)]

        b = CortexGraphBuilder()
        b.with_colors(*colors)
        for edge_color in edge_colors:
            b.add_edge_with_color_and_coverage('AAA', 'AAT', edge_color, color_coverage=1)
        b.add_edge_with_color_and_coverage('AAT', 'ATA', edge_colors[1], color_coverage=1)
        b.add_edge_with_color_and_coverage('AAT', 'ATC', edge_colors[2], color_coverage=1)

        # when
        graph = b.build()

        # then
        for node in ['ATA', 'ATC']:
            search = UnitigSearch.from_node_and_graph(node, graph)
            search.is_unitig_end(node, orientation)
        search = UnitigSearch.from_node_and_graph('AAA', graph)
        assert not search.is_unitig_end('AAA', EdgeTraversalOrientation.original)
        assert search.is_unitig_end('AAA', EdgeTraversalOrientation.reverse)
        assert search.is_unitig_end('AAT', EdgeTraversalOrientation.original)
        assert not search.is_unitig_end('AAT', EdgeTraversalOrientation.reverse)


class TestNextUnitigNodeIter(object):
    @given(s.sampled_from(('AAA', 'AAT', 'ATA')))
    def test_returns_two_nodes_in_order_regardless_of_starting_node(self, start_node):
        # given
        b = CortexGraphBuilder()
        b.with_colors(0)
        b.add_path('AAA', 'AAT', 'ATA')

        # when
        graph = b.build()

        # then
        search = UnitigSearch.from_node_and_graph(start_node, graph)
        assert ['AAT', 'ATA'] == list(
            search.next_unitig_node_iter(EdgeTraversalOrientation.original,
                                         start_node='AAA'))
        assert ['AAT', 'AAA'] == list(search.next_unitig_node_iter(EdgeTraversalOrientation.reverse,
                                                                   start_node='ATA'))

    def test_returns_two_nodes_in_order(self):
        # given
        b = CortexGraphBuilder()
        b.with_colors(0)
        b.add_path('AAA', 'AAT', 'ATA')

        # when
        graph = b.build()

        # then
        search = UnitigSearch.from_node_and_graph('AAA', graph)
        assert ['AAT', 'ATA'] == list(
            search.next_unitig_node_iter(EdgeTraversalOrientation.original))

        search = UnitigSearch.from_node_and_graph('ATA', graph)
        assert ['AAT', 'AAA'] == list(
            search.next_unitig_node_iter(EdgeTraversalOrientation.reverse))

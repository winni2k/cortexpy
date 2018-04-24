from hypothesis import given, settings
from hypothesis import strategies as s

from cortexpy.constants import EdgeTraversalOrientation
from cortexpy.test.driver.graph.find_unitigs import FindUnitigsTestDriver


class Test(object):
    def test_three_node_path_becomes_a_unitig(self):
        # given
        driver = FindUnitigsTestDriver()
        nodes = ['AAA', 'AAC', 'ACC']
        driver.graph.add_path(nodes)

        # when
        expect = driver.run()

        expect.has_n_nodes(1)
        expect.has_unitig_with_edges(('AAA', 'AAC'), ('AAC', 'ACC')).is_not_cycle()

    def test_three_node_path_with_mixed_node_order_becomes_a_unitig(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_edge('AAA', 'AAG')
        graph.add_edge('CAA', 'AAA')

        # when
        expect = driver.run()

        # then
        expect.has_unitig_with_edges(('AAA', 'AAG'), ('CAA', 'AAA')).is_not_cycle()

    def test_three_node_cycle_becomes_unitig(self):
        # given
        driver = FindUnitigsTestDriver()
        g = driver.graph
        edges = [('ACG', 'CGA'), ('CGA', 'GAC'), ('GAC', 'ACG')]
        for e in edges:
            g.add_edge(*e)

        # when
        expect = driver.run()

        # then
        expect.has_n_unitigs(1)
        expect.has_unitig_with_edges(('ACG', 'CGA'), ('CGA', 'GAC')).is_cycle()

    def test_two_node_path_becomes_unitig(self):
        # given
        driver = FindUnitigsTestDriver()
        driver.graph.add_path(['AAA', 'AAC'])

        # when
        expect = driver.run()

        # then
        expect.has_unitig_with_edges(('AAA', 'AAC'))

    def test_four_node_cycle_becomes_four_node_unitig(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_edge('AAA', 'AAC')
        graph.add_edge('AAC', 'ACA')
        graph.add_edge('ACA', 'CAA')
        graph.add_edge('CAA', 'AAA')

        # when
        expect = driver.run()

        # then
        (expect
         .has_n_nodes(1)
         .has_one_unitig()
         .has_unitig_with_edges(('AAA', 'AAC'), ('AAC', 'ACA'), ('ACA', 'CAA'))
         .with_left_node('AAA')
         .with_right_node('CAA')
         .is_cycle())

    def test_path_and_cycle_becomes_four_unitigs(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_edge('AAA', 'AAC')
        graph.add_edge('AAC', 'ACA')
        graph.add_edge('ACA', 'CAA')
        graph.add_edge('CAA', 'AAA')
        graph.add_edge('TAA', 'AAA')
        graph.add_edge('AAC', 'ACC')

        # when
        expect = driver.run()

        # then
        expect.has_n_unitigs(4)
        expect.has_unitig_with_edges(('AAA', 'AAC')).is_not_cycle()
        expect.has_unitig_with_one_node('TAA').is_not_cycle()
        expect.has_unitig_with_one_node('ACC').is_not_cycle()
        expect.has_unitig_with_edges(('ACA', 'CAA')).is_not_cycle()

    def test_two_node_path_and_four_node_cycle_becomes_two_unitigs(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_edge('AAA', 'AAC')
        graph.add_edge('AAC', 'ACA')
        graph.add_edge('ACA', 'CAA')
        graph.add_edge('CAA', 'AAA')
        graph.add_edge('TAA', 'AAA')
        graph.add_edge('ATA', 'TAA')

        # when
        expect = driver.run()

        # then
        expect.has_n_unitigs(2)
        (expect
         .has_unitig_with_edges(('ATA', 'TAA'))
         .with_left_node('ATA')
         .with_right_node('TAA')
         .is_not_cycle())
        (expect
         .has_unitig_with_edges(('AAA', 'AAC'), ('AAC', 'ACA'), ('ACA', 'CAA'))
         .with_left_node('AAA')
         .with_right_node('CAA')
         .is_cycle())

    def test_two_paths_making_bubble_results_in_four_unitigs(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_path(['CAA', 'AAA', 'AAT', 'ATC', 'TCC', 'CCC'])
        graph.add_path(['CAA', 'AAA', 'AAG', 'AGC', 'GCC', 'CCC'])

        # when
        expect = driver.run()

        # then
        expect \
            .has_unitig_with_edges(('CAA', 'AAA')) \
            .with_left_node('CAA') \
            .with_right_node('AAA') \
            .is_not_cycle()
        expect.has_unitig_with_edges(('AAT', 'ATC'), ('ATC', 'TCC')).with_left_node(
            'AAT').with_right_node('TCC').is_not_cycle()
        expect.has_unitig_with_edges(('AAG', 'AGC'), ('AGC', 'GCC')).with_left_node(
            'AAG').with_right_node('GCC').is_not_cycle()
        expect.has_unitig_with_one_node('CCC').is_not_cycle()
        expect.has_n_unitigs(4)


class TestIsUnitigEnd(object):
    def test_single_node_is_end_from_both_sides(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_node('AAA')
        driver.build()

        # when/then
        for orientation in EdgeTraversalOrientation:
            assert driver.finder.is_unitig_end('AAA', orientation)

    @given(s.integers(min_value=1, max_value=3))
    @settings(max_examples=3)
    def test_each_end_of_path_is_end(self, num_edges):
        # given
        colors = list(range(num_edges))
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.with_colors(*colors)
        for color in colors:
            graph.add_edge('AAA', 'AAC', key=color)
        finder = driver.build()

        # when/then
        assert finder.is_unitig_end('AAA', EdgeTraversalOrientation.reverse)
        assert not finder.is_unitig_end('AAA', EdgeTraversalOrientation.original)
        assert not finder.is_unitig_end('AAC', EdgeTraversalOrientation.reverse)
        assert finder.is_unitig_end('AAC', EdgeTraversalOrientation.original)

    def test_two_edges_into_one_node(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_edge('AAA', 'AAC')
        graph.add_edge('GAA', 'AAC')
        finder = driver.build()

        # when/then
        for node in ['AAA', 'AAC', 'GAA']:
            for orientation in EdgeTraversalOrientation:
                assert finder.is_unitig_end(node, orientation)

    def test_two_edges_out_of_one_node(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_edge('AAA', 'AAC')
        graph.add_edge('AAA', 'AAG')
        finder = driver.build()

        # when/then
        for node in ['AAA', 'AAC', 'AAG']:
            for orientation in EdgeTraversalOrientation:
                assert finder.is_unitig_end(node, orientation)

    def test_middle_unconnected_node_in_two_color_graph(self):
        # given
        driver = FindUnitigsTestDriver().with_colors(0, 1)
        graph = driver.graph
        nodes = ['AAA', 'AAC', 'ACC']
        graph.add_path(nodes, color=0, coverage=1)
        finder = driver.build()

        # when/then
        for node in nodes:
            for orientation in EdgeTraversalOrientation:
                assert finder.is_unitig_end(node, orientation)


class TestFindUnitigFromTwoColorGraph(object):
    def test_three_node_path_becomes_a_unitig_and_attributes_are_copied_across(self):
        # given
        nodes = ['AAA', 'AAC', 'ACC']
        colors = [0, 1]

        driver = FindUnitigsTestDriver()
        g_builder = driver.graph
        g_builder.with_colors(*colors)
        for color in colors:
            g_builder.add_path(nodes, color=color)
        for idx, node in enumerate(nodes):
            g_builder.with_node_coverage(node, (idx + 1, 1))
        finder = driver.builder.build()

        # when
        for start_node in nodes:
            unitig = finder.find_unitig_from(start_node)

            # then
            assert 'AAA' == unitig.left_node
            assert 'ACC' == unitig.right_node
            assert 3 == len(unitig.graph)
            assert {('AAA', 'AAC'), ('AAC', 'ACC')} == set(unitig.graph.edges(keys=False))
            for idx, node in enumerate(nodes):
                assert unitig.graph.node[node]['kmer'].coverage == (idx + 1, 1)


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
        colors = [0, 1]
        driver.with_colors(*colors)

        nodes = ['AAA', 'AAC']
        for node in nodes:
            driver.graph.add_node(node)

        for color, link_exists in zip(colors, [link_color_0_exists, link_color_1_exists]):
            if link_exists:
                driver.graph.add_edge(*nodes, key=color)
        driver.build()

        driver.graph.with_node_coverage(nodes[0], kmer_coverages[0])
        driver.graph.with_node_coverage(nodes[1], kmer_coverages[1])

        for idx, start_node in enumerate(nodes):
            # when
            unitig = driver.finder.find_unitig_from(start_node)

            # then
            if (
                (link_color_0_exists and link_color_1_exists) or
                (link_color_1_exists and kmer_coverages[0][0] == kmer_coverages[1][0] == 0) or
                (link_color_0_exists and kmer_coverages[0][1] == kmer_coverages[1][1] == 0)
            ):
                assert 2 == len(unitig.coverage)
                assert unitig.coverage[0] == kmer_coverages[0]
                assert unitig.coverage[1] == kmer_coverages[1]
            else:
                assert len(unitig.coverage) == 1
                assert unitig.coverage[0] == kmer_coverages[idx]

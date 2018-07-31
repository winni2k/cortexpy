import attr
from hypothesis import given, assume
from hypothesis import strategies as s

from cortexpy.graph.serializer.unitig import UnitigFinder
from cortexpy.test.builder.graph.cortex import CortexGraphBuilder
from cortexpy.test.driver.graph.find_unitigs import FindUnitigsTestDriver
from cortexpy.test.expectation.unitig_graph import UnitigExpectation


@attr.s(slots=True)
class FindUnitigFromTestDriver(object):
    start_node = attr.ib()
    builder = attr.ib(attr.Factory(CortexGraphBuilder))

    def __getattr__(self, item):
        return getattr(self.builder, item)

    def run(self):
        self.builder.make_consistent(self.start_node)
        finder = UnitigFinder.from_graph(self.builder.build())
        return UnitigExpectation(finder.find_unitig_from(self.start_node))


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
        expect.has_unitig_with_contig('AAACC')

    def test_three_node_path_with_mixed_node_order_becomes_a_unitig(self):
        # given
        driver = FindUnitigsTestDriver()
        graph = driver.graph
        graph.add_edge_with_coverage('AAA', 'AAG', 1)
        graph.add_edge_with_coverage('CAA', 'AAA', 1)

        # when
        expect = driver.run()

        # then
        expect.has_unitig_with_edges(('AAA', 'AAG'), ('CAA', 'AAA')).is_not_cycle()
        expect.has_unitig_with_contig('CAAAG')

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
        expect.has_unitig_with_contig('ACGAC')

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
        driver.with_seeds('AAA')

        # when
        expect = driver.run()

        # then
        expect \
            .has_n_nodes(1) \
            .has_one_unitig() \
            .has_unitig_with_edges(('AAA', 'AAC'), ('AAC', 'ACA'), ('ACA', 'CAA')) \
            .with_left_node('AAA') \
            .with_right_node('CAA') \
            .is_cycle()

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

    @given(s.data(),
           s.integers(min_value=1, max_value=23),
           s.sampled_from(('AAA', 'AAT', 'ATA', 'ATC')))
    def test_y_graph_becomes_three_unitigs(self, data, num_colors, start_node):
        # given
        colors = list(range(num_colors))
        edge_colors = [data.draw(s.sampled_from(colors)) for _ in range(3)]

        d = FindUnitigsTestDriver()
        d.with_colors(*colors)
        for color in colors:
            d.add_edge('AAA', 'AAT', color=color)
        d.add_edge('AAT', 'ATA', color=edge_colors[1])
        d.add_edge('AAT', 'ATC', color=edge_colors[1])

        # when
        expect = d.run()

        expect.has_unitig_with_contig('AAAT')
        expect.has_unitig_with_edges(('AAA', 'AAT'))
        expect.has_unitig_with_contig('ATA')
        expect.has_unitig_with_contig('ATC')
        expect.has_n_nodes(3)


class TestFindUnitigFrom(object):
    @given(s.data(),
           s.integers(min_value=1, max_value=3),
           s.sampled_from(('AAA', 'AAT', 'ATA', 'ATC')))
    def test_in_three_color_y_graph_always_finds_one_of_three_unitigs(self, data, num_colors,
                                                                      start_node):
        # given
        colors = list(range(num_colors))
        edge_colors = [data.draw(s.sampled_from(colors)) for _ in range(3)]

        contig_by_start_node = {'AAA': 'AAAT', 'AAT': 'AAAT', 'ATA': 'ATA', 'ATC': 'ATC'}
        d = FindUnitigFromTestDriver(start_node)
        d.with_colors(*colors)
        for edge_color in edge_colors:
            d.add_edge_with_color_and_coverage('AAA', 'AAT', edge_color, color_coverage=1)
        d.add_edge_with_color_and_coverage('AAT', 'ATA', edge_colors[1], color_coverage=1)
        d.add_edge_with_color_and_coverage('AAT', 'ATC', edge_colors[2], color_coverage=1)

        # when
        expect = d.run()

        expect.has_contig(contig_by_start_node[start_node])

    def test_three_node_path_becomes_a_unitig_and_attributes_are_copied_across_two_colors(self):
        # given
        nodes = ['AAA', 'AAC', 'ACC']
        colors = [0, 1]

        driver = FindUnitigsTestDriver()
        g_builder = driver.graph
        g_builder.with_colors(*colors)
        for color in colors:
            g_builder.add_path(nodes, color=color)
        for idx, node in enumerate(nodes):
            g_builder.with_node_coverage(node, idx + 1, 1)
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
    @given(s.data(), s.booleans(), s.booleans())
    def test_two_node_path_with_missing_links_are_not_joined(self, data, link_color_0_exists,
                                                             link_color_1_exists):
        kmer_coverages = [[0, 0], [0, 0]]
        if link_color_0_exists:
            kmer_coverages[0][0] = kmer_coverages[1][0] = 1
        else:
            kmer_coverages[0][0] = int(data.draw(s.booleans()))
            kmer_coverages[1][0] = int(data.draw(s.booleans()))
        if link_color_1_exists:
            kmer_coverages[0][1] = kmer_coverages[1][1] = 1
        else:
            kmer_coverages[0][1] = int(data.draw(s.booleans()))
            kmer_coverages[1][1] = int(data.draw(s.booleans()))

        kmer_coverages = [tuple(c) for c in kmer_coverages]
        assume(all(c != (0, 0) for c in kmer_coverages))

        # given
        driver = FindUnitigsTestDriver()
        colors = [0, 1]
        driver.with_colors(*colors)

        nodes = ['AAA', 'AAC']
        for idx, node in enumerate(nodes):
            driver.with_node_coverage(node, *kmer_coverages[idx])
        for color, link_exists in zip(colors, [link_color_0_exists, link_color_1_exists]):
            if link_exists:
                driver.add_edge_with_color_and_coverage(*nodes, color, 1)

        driver.build()

        for idx, start_node in enumerate(nodes):
            # when
            unitig = driver.finder.find_unitig_from(start_node)

            # then
            if (
                (link_color_0_exists and link_color_1_exists) or
                (link_color_1_exists and kmer_coverages[0][0] == kmer_coverages[1][0] == 0) or
                (link_color_0_exists and kmer_coverages[0][1] == kmer_coverages[1][1] == 0)
            ):
                assert unitig.coverage[0] == kmer_coverages[0]
                assert unitig.coverage[1] == kmer_coverages[1]
                assert 2 == len(unitig.coverage)
            else:
                assert len(unitig.coverage) == 1
                assert unitig.coverage[0] == kmer_coverages[idx]

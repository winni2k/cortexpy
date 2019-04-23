import pytest

import cortexpy.graph
import cortexpy.graph.parser
from cortexpy.graph.parser.random_access import SlurpedRandomAccess
from cortexpy.test.driver.graph.traversal import EngineTestDriver
from cortexpy.test.expectation import KmerGraphExpectation


@pytest.fixture(params=('slurped', 'not_slurped'))
def driver(request):
    if request.param == 'slurped':
        return EngineTestDriver(ra_constructor=SlurpedRandomAccess.from_handle)
    return EngineTestDriver()


class Test:
    def test_raises_on_empty(self, driver):
        # given
        driver \
            .with_kmer_size(3) \
            .with_start_kmer_string('AAA')

        # when/then
        with pytest.raises(KeyError):
            driver.run()

        assert len(driver.traverser.graph) == 0

    def test_three_connected_kmers_returns_graph_with_three_kmers(self, driver):
        # given
        (driver
         .with_kmer_size(3)
         .with_kmer('AAA 1 .......T')
         .with_kmer('AAT 1 a....C..')
         .with_kmer('ATC 1 a.......')
         .with_start_kmer_string('AAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('AAA', 'AAT', 'ATC')
         .has_edges('AAA AAT 0', 'AAT ATC 0'))

    @pytest.mark.parametrize('traversal_color', (0, 1))
    def test_three_connected_kmers_returns_graph_with_three_kmers_for_two_colors(self,
                                                                                 traversal_color,
                                                                                 driver):
        # given
        driver \
            .with_kmer_size(3) \
            .with_num_colors(2) \
            .with_traversal_colors(traversal_color) \
            .with_start_kmer_string('AAA')
        driver \
            .with_kmer('AAA 1 1 .......T .......T') \
            .with_kmer('AAT 1 1 a....C.. a....C..') \
            .with_kmer('ATC 1 1 a....... a.......')

        # when
        expect = driver.run()

        # then
        for node in ['AAA', 'AAT', 'ATC']:
            expect.has_node(node).has_coverages(1, 1)
        expect.has_n_nodes(3)
        expect.has_edges('AAA AAT 0',
                         'AAA AAT 1',
                         'AAT ATC 0',
                         'AAT ATC 1')

    def test_four_connected_kmers_in_star_returns_graph_with_four_kmers(self, driver):
        # given
        (driver
         .with_kmer_size(3)
         .with_kmer('AAA 1 .......T')
         .with_kmer('AAT 1 a....CG.')
         .with_kmer('ATC 1 a.......')
         .with_kmer('ATG 1 a.......')
         .with_start_kmer_string('AAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('AAA', 'AAT', 'ATC', 'ATG')
         .has_edges('AAA AAT 0', 'AAT ATC 0', 'AAT ATG 0'))

    def test_cycle_is_traversed_once(self, driver):
        # given
        (driver
         .with_kmer_size(3)
         .with_kmer('CAA 1 ....A...')
         .with_kmer('AAA 1 .c.t...T')
         .with_kmer('AAT 1 a...A...')
         .with_kmer('ATA 1 a...A...')
         .with_kmer('TAA 1 a...A...')
         .with_start_kmer_string('CAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('CAA', 'AAA', 'AAT', 'ATA', 'TAA')
         .has_n_edges(5))

    def test_cycle_and_branch_is_traversed_once(self, driver):
        # given
        (driver
         .with_kmer_size(3)
         .with_kmer('CAA 1 ....A...')
         .with_kmer('AAA 1 .c.t.C.T')
         .with_kmer('AAC 1 a.......')
         .with_kmer('AAT 1 a...A...')
         .with_kmer('ATA 1 a...A...')
         .with_kmer('TAA 1 a...A...')
         .with_start_kmer_string('CAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('CAA', 'AAA', 'AAT', 'ATA', 'TAA', 'AAC')
         .has_n_edges(6))

    def test_two_cycles_are_traversed_once(self, driver):
        # given
        (driver
         .with_kmer_size(3)
         .with_kmer('CAA 1 a...A...')
         .with_kmer('AAA 1 .c.t.C.T')
         .with_kmer('AAC 1 a...A...')
         .with_kmer('ACA 1 a...A...')
         .with_kmer('AAT 1 a...A...')
         .with_kmer('ATA 1 a...A...')
         .with_kmer('TAA 1 a...A...')
         .with_start_kmer_string('CAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('CAA', 'AAA', 'AAT', 'ATA', 'TAA', 'AAC', 'ACA')
         .has_n_edges(8))

    def test_two_cycles_are_traversed_once_in_revcomp(self, driver):
        # given
        (driver
         .with_kmer_size(3)
         .with_kmer('CAA 1 a...A...')
         .with_kmer('AAA 1 .c.t.C.T')
         .with_kmer('AAC 1 a...A...')
         .with_kmer('ACA 1 a...A...')
         .with_kmer('AAT 1 a...A...')
         .with_kmer('ATA 1 a...A...')
         .with_kmer('TAA 1 a...A...')
         .with_start_kmer_string('TTG'))

        # when
        expect = driver.run()

        # then
        expect \
            .has_nodes('TAA', 'ATA', 'AAT', 'ACA', 'AAC', 'AAA', 'CAA') \
            .has_n_edges(8)

    def test_exploration_of_two_colors_returns_all_kmers(self, driver):
        # given
        driver.with_kmer_size(3) \
            .with_num_colors(2) \
            .with_traversal_colors(0, 1) \
            .with_start_kmer_string('AAA')
        driver \
            .with_kmer('AAA 0 1 ........ .......T') \
            .with_kmer('AAT 1 1 .....C.. a.......') \
            .with_kmer('ATC 1 0 a....... ........')

        # when
        expect = driver.run()

        # then
        expect.has_node_coverages('AAA 0 1',
                                  'AAT 1 1',
                                  'ATC 1 0')
        expect.has_edges('AAA AAT 1',
                         'AAT ATC 0')

    def test_exploration_of_two_colors_on_initial_branch_point_returns_all_kmers(self, driver):
        # given
        driver.with_kmer_size(3) \
            .with_num_colors(2) \
            .with_traversal_colors(0, 1) \
            .with_start_kmer_string('AAA')
        driver \
            .with_kmer('AAA 1 1 .......T .....C..') \
            .with_kmer('AAT 1 0 a....... ........') \
            .with_kmer('AAC 0 1 ........ a.......')

        # when
        expect = driver.run()

        # then
        expect.has_node_coverages('AAA 1 1',
                                  'AAT 1 0',
                                  'AAC 0 1')
        expect.has_edges('AAA AAT 0',
                         'AAA AAC 1', )


class TestTraversalOrientationBoth:
    def test_with_three_linked_kmers_returns_graph_of_three_kmers(self, driver):
        # given
        (driver
         .with_kmer_size(3)
         .with_kmer('AAA 1 .......T')
         .with_kmer('AAT 1 a....C..')
         .with_kmer('ATC 1 a.......')
         .with_start_kmer_string('ATC')
         .with_traversal_orientation('both'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('AAA', 'AAT', 'ATC')
         .has_n_edges(2))


class TestReverseOrientation:
    def test_three_connected_kmers_returns_graph_with_three_kmers(self, driver):
        # given
        (driver
         .with_kmer_size(3)
         .with_kmer('AAA 1 .......T')
         .with_kmer('AAT 1 a....C..')
         .with_kmer('ATC 1 a.......')
         .with_start_kmer_string('ATC')
         .with_traversal_orientation('reverse'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('AAA', 'AAT', 'ATC')
         .has_edges(('AAA', 'AAT', 0), ('AAT', 'ATC', 0)))

    @pytest.mark.parametrize('start_kmer_string', ('AAA', 'AAT', 'ATA', 'TAA'))
    def test_cycle_is_traversed_once(self, driver, start_kmer_string):
        # given
        (driver
         .with_kmer_size(3)
         .with_kmer('CAA 1 ....A...')
         .with_kmer('AAA 1 .c.t...T')
         .with_kmer('AAT 1 a...A...')
         .with_kmer('ATA 1 a...A...')
         .with_kmer('TAA 1 a...A...')
         .with_traversal_orientation('reverse')
         .with_start_kmer_string(start_kmer_string))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('CAA', 'AAA', 'AAT', 'ATA', 'TAA')
         .has_n_edges(5))


class TestBothOrientation(object):
    @pytest.mark.parametrize('start_kmer_string', ('CCA', 'CAA', 'AAA', 'AAT', 'ATA', 'TAA'))
    def test_cycle_and_branch_are_traversed_once(self, driver, start_kmer_string):
        # given
        driver \
            .with_kmer_size(3) \
            .with_kmer('CCA 1 ....A...') \
            .with_kmer('CAA 1 .c..A...') \
            .with_kmer('AAA 1 .c.t...T') \
            .with_kmer('AAT 1 a...A...') \
            .with_kmer('ATA 1 a...A...') \
            .with_kmer('TAA 1 a...A...') \
            .with_traversal_orientation('both') \
            .with_start_kmer_string(start_kmer_string)

        # when
        expect = driver.run()

        # then
        expect \
            .has_nodes('CCA', 'CAA', 'AAA', 'AAT', 'ATA', 'TAA') \
            .has_edges('CAA CCA 0',
                       'AAA CAA 0',
                       'AAA AAT 0',
                       'AAT ATA 0',
                       'ATA TAA 0',
                       'AAA TAA 0')

    @pytest.mark.parametrize('start_kmer_index', range(5))
    def test_star_with_two_colors(self, driver, start_kmer_index):
        kmers = ('CAA', 'GAA', 'AAA', 'AAT', 'AAC')
        start_kmer_string = kmers[start_kmer_index]

        # given
        (driver
         .with_kmer_size(3)
         .with_num_colors(2)
         .with_kmer('CAA 1 1 ....A... ....A...')
         .with_kmer('GAA 1 1 ....A... ....A...')
         .with_kmer('AAA 1 1 .cg..C.T .cg.....')
         .with_kmer('AAT 1 0 a....... ........')
         .with_kmer('AAC 1 0 a....... ........')
         .with_traversal_orientation('both')
         .with_start_kmer_string(start_kmer_string))

        # when
        expect = driver.run()

        # then
        for node in ['CAA', 'GAA', 'AAA']:
            expect.has_node(node).has_coverages(1, 1)
        for node in ['AAT', 'AAC']:
            expect.has_node(node).has_coverages(1, 0)
        expect \
            .has_n_nodes(len(kmers)) \
            .has_edges('AAA CAA 0',
                       'AAA CAA 1',
                       'AAA GAA 0',
                       'AAA GAA 1',
                       'AAA AAT 0',
                       'AAA AAC 0')

    @pytest.mark.parametrize('start_kmer_index', range(4))
    def test_branch_with_two_colors(self, driver, start_kmer_index):
        kmers = ('CAA', 'GAA', 'AAA', 'AAC')
        start_kmer_string = kmers[start_kmer_index]

        # given
        driver \
            .with_kmer_size(3) \
            .with_num_colors(2) \
            .with_kmer('CAA 0 1 ........ ....A...') \
            .with_kmer('GAA 1 0 ....A... ........') \
            .with_kmer('AAA 1 1 ..g..C.. .c......') \
            .with_kmer('AAC 1 0 a....... ........') \
            .with_traversal_orientation('both') \
            .with_start_kmer_string(start_kmer_string) \
            .with_traversal_colors(0, 1)

        # when
        expect = driver.run()

        # then
        expect.has_node('CAA').has_coverages('0 1')
        expect.has_node('GAA').has_coverages('1 0')
        expect.has_node('AAA').has_coverages('1 1')
        expect.has_node('AAC').has_coverages('1 0')
        expect \
            .has_n_nodes(len(kmers)) \
            .has_edges('AAA CAA 1',
                       'AAA GAA 0',
                       'AAA AAC 0')


class TestTraverseFromEachKmerIn:
    def test_does_not_raise_on_empty(self, driver):
        # given
        (driver
         .with_kmer_size(3)
         .with_start_string('AAA'))

        # when/then
        driver.run()

        assert len(driver.traverser.graph) == 0

    def test_cycle_and_branch_are_traversed_once(self, driver):
        # given
        start_string = 'CCAAATAA'
        (driver
         .with_kmer_size(3)
         .with_kmer('CCA 1 ....A...')
         .with_kmer('CAA 1 .c..A...')
         .with_kmer('AAA 1 .c.t...T')
         .with_kmer('AAT 1 a...A...')
         .with_kmer('ATA 1 a...A...')
         .with_kmer('TAA 1 a...A...')
         .with_traversal_orientation('both')
         .with_start_string(start_string))

        # when
        expect = driver.run()

        # then
        expect \
            .has_nodes('CCA', 'CAA', 'AAA', 'AAT', 'ATA', 'TAA') \
            .has_edges('CAA CCA 0',
                       'AAA CAA 0',
                       'AAA AAT 0',
                       'AAT ATA 0',
                       'ATA TAA 0',
                       'AAA TAA 0')


class TestMaxNodes(object):
    def test_of_two_returns_with_two_nodes_plus_edges(self, driver):
        # given
        (driver
         .with_kmer_size(3)
         .with_max_nodes(2)
         .with_kmer('AAA 1 .......T')
         .with_kmer('AAT 1 a....CG.')
         .with_kmer('ATC 1 a......T')
         .with_kmer('ATG 1 a.......')
         .with_kmer('AGA 1 .......T')
         .with_start_kmer_string('AAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('AAA', 'AAT', 'ATC', 'ATG')
         .has_n_edges(3))


class TestStartStringSize:
    def test_raises_when_string_wrong_size(self, driver):
        for start_string in ['AAA', 'AAAAAAA']:
            # given
            driver.with_kmer_size(5).with_start_kmer_string(start_string)

            with pytest.raises(AssertionError):
                driver.run()


class TestEdgeAnnotation:
    def test_with_single_kmer_and_link_annotates_etra_links(self, driver):
        driver.with_kmer_size(3) \
            .with_kmer('AAA 1 1 ........ .c...C..') \
            .with_start_kmer_string('AAA') \
            .with_traversal_colors(0)

        # when
        graph_expectation = driver.run()
        annotated_graph = cortexpy.graph.traversal.engine \
            .annotate_kmer_graph_edges(graph_expectation.graph)
        expect = KmerGraphExpectation(annotated_graph)

        # then
        expect.has_node('CAA').has_coverages(1, 1)
        expect.has_node('AAA').has_coverages(1, 1)
        expect.has_node('AAC').has_coverages(1, 1)
        expect.has_n_nodes(3)
        expect.has_edges('AAA AAC 1', 'AAA CAA 1')

    def test_for_revcomp_with_single_kmer_and_link_annotates_etra_links(self, driver):
        driver \
            .with_kmer_size(3) \
            .with_kmer('AAA 1 1 ........ ..g..C..') \
            .with_traversal_colors(0) \
            .with_start_kmer_string('TTT')

        # when
        graph_expectation = driver.run()
        annotated_graph = cortexpy.graph.traversal.engine \
            .annotate_kmer_graph_edges(graph_expectation.graph)
        expect = KmerGraphExpectation(annotated_graph)

        # then
        expect.has_node('TTC').has_coverages(1, 1)
        expect.has_node('TTT').has_coverages(1, 1)
        expect.has_node('GTT').has_coverages(1, 1)
        expect.has_n_nodes(3)
        expect.has_edges('AAA GAA 1', 'AAA AAC 1')


class TestFixture:
    def test_multi_color_traversal_bug(self, driver):
        start_kmer_string = 'ATCTGAT'
        driver \
            .with_kmer_size(7) \
            .with_kmer('ATCAGAT 1 1 .....C.. .....C..') \
            .with_kmer('GATCTGA 1 1 a......T ..g....T') \
            .with_kmer('AGATCTG 1 0 ....A... ........') \
            .with_kmer('CAGATCC 0 1 ........ ...t....') \
            .with_start_kmer_string(start_kmer_string) \
            .with_traversal_colors(0, 1) \
            .with_traversal_orientation('both')

        # when
        graph_expectation = driver.run()
        annotated_graph = cortexpy.graph.traversal.engine \
            .annotate_kmer_graph_edges(graph_expectation.graph)
        expect = KmerGraphExpectation(annotated_graph)

        # then
        # expect.has_node('CAA').has_coverages(0, 0)
        # expect.has_node('AAA').has_coverages(1, 1)
        expect.has_node('GGATCTG').has_coverages('0 1')
        expect.has_n_nodes(4)
        # expect.has_edges('AAA AAC 1', 'CAA AAA 1')

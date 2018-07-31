from unittest import mock

import attr
import pytest
from hypothesis import given, strategies as s

from cortexpy.graph.traversal import branch
import cortexpy.graph.parser
import cortexpy.test.builder as builder
import cortexpy.test.expectation as expectation
from cortexpy.constants import EdgeTraversalOrientation


@attr.s(slots=True)
class BranchExpectation(object):
    traversed_branch = attr.ib()
    start_kmer_string = attr.ib()
    mocked_rap_getitem = attr.ib(None)
    graph_expectation = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.graph_expectation = expectation.graph.KmerGraphExpectation(self.traversed_branch.graph)

    def __getattr__(self, item):
        return getattr(self.graph_expectation, item)

    def has_first_kmer_string(self, kmer_string):
        assert self.traversed_branch.first_kmer_string == kmer_string
        return self

    def has_last_kmer_string(self, kmer_string):
        assert kmer_string in self.traversed_branch.graph
        assert self.traversed_branch.last_kmer_string == kmer_string
        self.has_first_kmer_string(self.start_kmer_string)
        return self

    def has_neighbor_kmer_strings(self, *neighbor_kmer_strings):
        assert set(self.traversed_branch.neighbor_kmer_strings) == set(neighbor_kmer_strings)
        return self

    def has_reverse_neighbor_kmer_strings(self, *neighbor_kmer_strings):
        assert set(self.traversed_branch.reverse_neighbor_kmer_strings) == set(
            neighbor_kmer_strings)
        return self

    def called_stream_read_n_times(self, n):
        if self.mocked_rap_getitem is None:
            raise Exception('RandomAccess parser was not mocked')
        assert n == self.mocked_rap_getitem.call_count
        return self


@attr.s(slots=True)
class BranchTestDriver(object):
    graph_builder = attr.ib(attr.Factory(builder.Graph))
    start_kmer_string = attr.ib(None)
    traversal_color = attr.ib(0)
    other_stopping_colors = attr.ib(attr.Factory(set))
    traversal_orientation = attr.ib(EdgeTraversalOrientation.original)
    parent_graph = attr.ib(attr.Factory(set))
    expected_start_kmer_string = attr.ib(None)
    warmup_ra_parser = attr.ib(False)

    def with_kmer(self, *args, **kwargs):
        self.graph_builder.with_kmer(*args, **kwargs)
        return self

    def with_kmer_size(self, n):
        self.graph_builder.with_kmer_size(n)
        return self

    def with_num_colors(self, n):
        self.graph_builder.with_num_colors(n)
        return self

    def with_start_kmer_string(self, start_kmer_string):
        self.start_kmer_string = start_kmer_string
        return self

    def with_traversal_color(self, color):
        self.traversal_color = color
        return self

    def with_other_stopping_colors(self, *colors):
        self.other_stopping_colors = set(colors)
        return self

    def with_parent_graph_nodes(self, *graph_nodes):
        self.parent_graph = set(graph_nodes)
        return self

    def with_reverse_traversal_orientation(self):
        self.traversal_orientation = EdgeTraversalOrientation.reverse
        return self

    def with_ra_parser_warmup(self):
        self.warmup_ra_parser = True
        return self

    def run(self):
        if self.expected_start_kmer_string is None:
            self.expected_start_kmer_string = self.start_kmer_string
        stream = self.graph_builder.build()
        ra_parser = cortexpy.graph.parser.random_access.RandomAccess(stream)
        traverser = branch.Traverser(
            ra_parser,
            traversal_color=self.traversal_color,
            other_stopping_colors=self.other_stopping_colors)
        if self.warmup_ra_parser:
            for k_string in list(ra_parser):
                ra_parser[k_string]
        with mock.patch.object(stream, 'read', wraps=stream.read) as mocked_read:
            traversed_branch = traverser \
                .traverse_from(self.start_kmer_string,
                               parent_graph=self.parent_graph,
                               orientation=self.traversal_orientation)
            return BranchExpectation(traversed_branch,
                                     start_kmer_string=self.expected_start_kmer_string,
                                     mocked_rap_getitem=mocked_read)


class Test(object):
    def test_raises_on_empty_graph_returns_empty_graph(self):
        # given
        driver = BranchTestDriver().with_kmer_size(3).with_start_kmer_string('AAA')

        # when
        with pytest.raises(KeyError):
            driver.run()

    def test_two_unconnected_kmers_returns_graph_with_one_kmer(self):
        # given
        driver = BranchTestDriver()
        driver \
            .with_kmer_size(3) \
            .with_kmer('AAA') \
            .with_kmer('AAT') \
            .with_start_kmer_string('AAA')

        # when
        expect = driver.run()

        # then
        expect \
            .has_nodes('AAA') \
            .has_n_edges(0)

    def test_two_connected_kmers_returns_graph_with_two_kmers(self):
        # given
        driver = BranchTestDriver()
        (driver
         .with_kmer_size(3)
         .with_kmer('AAA 1 .......T')
         .with_kmer('AAT 1 a.......')
         .with_start_kmer_string('AAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('AAA', 'AAT')
         .has_edges(('AAA', 'AAT', 0)))

    def test_two_connected_kmers_returns_graph_with_two_kmers_and_two_colors(self):
        # given
        driver = BranchTestDriver()
        (driver
         .with_kmer_size(3)
         .with_kmer('AAA 1 1 .......T .......T')
         .with_kmer('AAT 1 1 a....... a.......')
         .with_start_kmer_string('AAA'))

        # when
        expect = driver.run()

        # then
        expect.has_node('AAA').has_coverages(1, 1)
        expect.has_node('AAT').has_coverages(1, 1)
        expect.has_n_nodes(2)
        expect.has_edges('AAA AAT 0', 'AAA AAT 1')

    def test_two_connected_kmers_with_other_edges_returns_graph_with_one_kmer(self):
        # given
        driver = BranchTestDriver()
        (driver
         .with_kmer_size(3)
         .with_kmer('AAA 1 ac..A..T')
         .with_kmer('AAT 1 a....C.T')
         .with_start_kmer_string('AAA'))

        # when
        expect = driver.run()

        # then
        expect \
            .has_neighbor_kmer_strings('AAT', 'AAA') \
            .has_reverse_neighbor_kmer_strings('AAA', 'CAA') \
            .has_nodes('AAA') \
            .has_edges('AAA AAA 0')

    def test_raises_when_two_connected_kmers_with_missing_exiting_kmer(self):
        # given
        driver = BranchTestDriver()
        (driver
         .with_kmer_size(3)
         .with_kmer('AAA 1 .......T')
         .with_kmer('AAT 1 a....C..')
         .with_start_kmer_string('AAA'))

        # when/then
        with pytest.raises(KeyError):
            driver.run()

    def test_three_connected_kmers_returns_graph_with_three_kmers(self):
        # given
        driver = BranchTestDriver()
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
         .has_n_edges(2))

    def test_three_connected_kmers_returns_graph_with_three_kmers_as_revcomp(self):
        # given
        driver = BranchTestDriver()
        (driver
         .with_kmer_size(3)
         .with_kmer('AAA 1 .......T')
         .with_kmer('AAT 1 a....C..')
         .with_kmer('ATC 1 a.......')
         .with_start_kmer_string('GAT'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('GAT', 'ATT', 'TTT')
         .has_n_edges(2))

    def test_four_kmers_one_revcomp_returns_graph_with_four_kmers_as_revcomp(self):
        # given
        driver = BranchTestDriver()
        (driver
         .with_kmer_size(3)
         .with_kmer('AAA 1 .......T')
         .with_kmer('AAT 1 a....C..')
         .with_kmer('ATC 1 a.....G.')
         .with_kmer('CGA 1 .......T')
         .with_start_kmer_string('CGA'))

        # when
        expect = driver.run()

        # then
        expect \
            .has_last_kmer_string('TTT') \
            .has_nodes('CGA', 'GAT', 'ATT', 'TTT') \
            .has_n_edges(3)

    @given(s.sampled_from(('AAA', 'AAT', 'ATA', 'TAA')), s.booleans(),
           s.integers(min_value=1, max_value=3))
    def test_cycle_is_traversed_once(self, start_kmer_string, reverse_orientation, num_colors):
        # given
        d = BranchTestDriver()
        d \
            .with_kmer_size(3) \
            .with_kmer('AAA', 1, '...t...T', repeat_color_edges_n_times=num_colors) \
            .with_kmer('AAT', 1, 'a...A...', repeat_color_edges_n_times=num_colors) \
            .with_kmer('ATA', 1, 'a...A...', repeat_color_edges_n_times=num_colors) \
            .with_kmer('TAA', 1, 'a...A...', repeat_color_edges_n_times=num_colors) \
            .with_start_kmer_string(start_kmer_string) \
            .with_ra_parser_warmup()

        if reverse_orientation:
            d.with_reverse_traversal_orientation()

        # when
        expect = d.run()

        # then
        expect \
            .has_nodes('AAA', 'AAT', 'ATA', 'TAA') \
            .has_n_edges(4 * num_colors)
        expect.called_stream_read_n_times(8)

    def test_one_seen_kmer_returns_empty_graph(self):
        # given
        driver = BranchTestDriver()
        (driver
         .with_kmer_size(3)
         .with_kmer('AAA')
         .with_parent_graph_nodes('AAA')
         .with_start_kmer_string('AAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_first_kmer_string(None)
         .has_nodes()
         .has_n_edges(0))


class TestReverseOrientation(object):
    def test_two_connected_kmers_returns_graph_with_two_kmers(self):
        # given
        driver = BranchTestDriver()
        (driver
         .with_kmer_size(3)
         .with_kmer('AAA 1 .......T')
         .with_kmer('AAT 1 a.......')
         .with_start_kmer_string('AAT')
         .with_reverse_traversal_orientation())

        # when
        expect = driver.run()

        # then
        (expect
         .has_last_kmer_string('AAA')
         .has_nodes('AAA', 'AAT')
         .has_edges(('AAA', 'AAT', 0))
         .has_n_edges(1))


class TestJunctions(object):
    def test_stops_at_junction_with_two_in_one_out(self):
        # given
        driver = BranchTestDriver()
        (driver
         .with_kmer_size(3)
         .with_kmer('CAA 1 ....A...')
         .with_kmer('TAA 1 ....A...')
         .with_kmer('AAA 1 .c.t...T')
         .with_kmer('AAT 1 a.......')
         .with_start_kmer_string('CAA'))
        expected_nodes = ['AAA', 'CAA']
        expected_neighbor_kmer_strings = ['AAT']
        expected_reverse_neighbor_kmer_strings = ['TAA']

        # when
        expect = driver.run()

        # then
        expect.has_nodes(*expected_nodes)
        expect \
            .has_last_kmer_string('AAA') \
            .has_neighbor_kmer_strings(*expected_neighbor_kmer_strings) \
            .has_reverse_neighbor_kmer_strings(*expected_reverse_neighbor_kmer_strings) \
            .has_edges('CAA AAA 0')

    def test_when_reverse_traverse_stops_at_junction_with_two_out_one_in(self):
        # given
        driver = BranchTestDriver()
        (driver
         .with_kmer_size(3)
         .with_kmer('TAA 1 ....A...')
         .with_kmer('AAA 1 ...t.C.T')
         .with_kmer('AAT 1 a.......')
         .with_kmer('AAC 1 a.......')
         .with_reverse_traversal_orientation()
         .with_start_kmer_string('AAC'))
        expected_nodes = ['AAA', 'AAC']
        expected_neighbor_kmer_strings = ['TAA']
        expected_reverse_neighbor_kmer_strings = ['AAT']

        # when
        expect = driver.run()

        # then
        expect.has_nodes(*expected_nodes)
        (expect
         .has_last_kmer_string('AAA')
         .has_neighbor_kmer_strings(*expected_neighbor_kmer_strings)
         .has_reverse_neighbor_kmer_strings(*expected_reverse_neighbor_kmer_strings)
         .has_n_edges(len(expected_nodes) - 1))

        if len(expected_nodes) == 2:
            expect.has_edges(('AAA', 'AAC', 0))

    def test_stops_at_junction_with_two_in(self):
        # given
        driver = BranchTestDriver()
        driver \
            .with_kmer_size(3) \
            .with_kmer('CAA 1 ....A...') \
            .with_kmer('TAA 1 ....A...') \
            .with_kmer('AAA 1 .c.t....') \
            .with_start_kmer_string('CAA')
        expected_reverse_neighbor_kmer_strings = ['TAA']

        # when
        expect = driver.run()

        # then
        expect.has_nodes('CAA', 'AAA')
        expect \
            .has_last_kmer_string('AAA') \
            .has_neighbor_kmer_strings() \
            .has_reverse_neighbor_kmer_strings(*expected_reverse_neighbor_kmer_strings) \
            .has_n_edges(1)

    def test_with_continuing_kmer_stops_at_junction_with_second_color(self):
        # given
        driver = BranchTestDriver()
        driver \
            .with_kmer_size(3) \
            .with_kmer('CAA 0 1 ........ ....A...') \
            .with_kmer('GAA 1 0 ....A... ........') \
            .with_kmer('AAA 1 1 ..g..C.. .c......') \
            .with_kmer('AAC 1 0 a....... ........') \
            .with_start_kmer_string('AAC') \
            .with_other_stopping_colors(1) \
            .with_reverse_traversal_orientation()

        # when
        expect = driver.run()

        # then
        expect.has_nodes('AAA', 'AAC')
        expect \
            .has_last_kmer_string('AAA') \
            .has_neighbor_kmer_strings('GAA') \
            .has_reverse_neighbor_kmer_strings() \
            .has_n_edges(1)

    def test_without_continuing_kmer_stops_at_junction_with_second_color(self):
        # given
        driver = BranchTestDriver()
        driver \
            .with_kmer_size(3) \
            .with_kmer('CAA 0 1 ........ ....A...') \
            .with_kmer('AAA 1 1 .....C.. .c......') \
            .with_kmer('AAC 1 0 a....... ........') \
            .with_start_kmer_string('AAC') \
            .with_other_stopping_colors(1) \
            .with_reverse_traversal_orientation()

        # when
        expect = driver.run()

        # then
        expect.has_nodes('AAA', 'AAC')
        expect \
            .has_last_kmer_string('AAA') \
            .has_neighbor_kmer_strings() \
            .has_reverse_neighbor_kmer_strings() \
            .has_n_edges(1)

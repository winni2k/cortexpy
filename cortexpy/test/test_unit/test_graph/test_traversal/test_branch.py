import attr
import pytest

import cortexpy.graph.parser
import cortexpy.test.builder as builder
import cortexpy.test.expectation as expectation
from cortexpy.graph.serializer import EdgeTraversalOrientation


@attr.s(slots=True)
class BranchExpectation(object):
    traversed_branch = attr.ib()
    start_kmer_string = attr.ib()
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


@attr.s(slots=True)
class BranchTestDriver(object):
    graph_builder = attr.ib(attr.Factory(builder.Graph))
    start_kmer_string = attr.ib(None)
    traversal_color = attr.ib(0)
    traversal_orientation = attr.ib(EdgeTraversalOrientation.original)
    parent_graph = attr.ib(attr.Factory(set))

    def with_kmer(self, *args):
        self.graph_builder.with_kmer(*args)
        return self

    def with_kmer_size(self, n):
        self.graph_builder.with_kmer_size(n)
        return self

    def with_start_kmer_string(self, start_kmer_string):
        self.start_kmer_string = start_kmer_string
        return self

    def with_traversal_color(self, color):
        self.traversal_color = color
        return self

    def with_parent_graph_nodes(self, *graph_nodes):
        self.parent_graph = set(graph_nodes)
        return self

    def with_reverse_traversal_orientation(self):
        self.traversal_orientation = EdgeTraversalOrientation.reverse
        return self

    def run(self):
        assert self.start_kmer_string is not None
        random_access_parser = cortexpy.graph.parser.RandomAccess(self.graph_builder.build())
        traversed_branch = (cortexpy.graph.traversal
                            .Branch(random_access_parser,
                                    self.traversal_color)
                            .traverse_from(self.start_kmer_string,
                                           parent_graph=self.parent_graph,
                                           orientation=self.traversal_orientation))
        return BranchExpectation(traversed_branch, start_kmer_string=self.start_kmer_string)


class Test(object):
    def test_raises_on_empty_graph_returns_empty_graph(self):
        # given
        driver = BranchTestDriver().with_kmer_size(3).with_start_kmer_string('AAA')

        # when
        with pytest.raises(KeyError):
            driver.run()

    def test_two_unconnected_kmers_returns_graph_with_one_kmer(self):
        # given
        driver = (BranchTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA')
                  .with_kmer('AAT')
                  .with_start_kmer_string('AAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('AAA')
         .has_n_edges(0))

    def test_two_connected_kmers_returns_graph_with_two_kmers(self):
        # given
        driver = (BranchTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 0, '.......T')
                  .with_kmer('AAT', 0, 'a.......')
                  .with_start_kmer_string('AAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('AAA', 'AAT')
         .has_edges(('AAA', 'AAT', 0)))

    def test_two_connected_kmers_with_other_edges_returns_graph_with_two_kmers(self):
        # given
        driver = (BranchTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 0, 'ac.....T')
                  .with_kmer('AAT', 0, 'a....C.T')
                  .with_start_kmer_string('AAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_neighbor_kmer_strings('ATC', 'ATT')
         .has_nodes('AAA', 'AAT')
         .has_n_edges(1))

    def test_raises_when_two_connected_kmers_with_missing_exiting_kmer(self):
        # given
        driver = (BranchTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 0, '.......T')
                  .with_kmer('AAT', 0, 'a....C..')
                  .with_start_kmer_string('AAA'))

        # when/then
        with pytest.raises(KeyError):
            driver.run()

    def test_three_connected_kmers_returns_graph_with_three_kmers(self):
        # given
        driver = (BranchTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 0, '.......T')
                  .with_kmer('AAT', 0, 'a....C..')
                  .with_kmer('ATC', 0, 'a.......')
                  .with_start_kmer_string('AAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('AAA', 'AAT', 'ATC')
         .has_n_edges(2))

    def test_three_connected_kmers_returns_graph_with_three_kmers_as_revcomp(self):
        # given
        driver = (BranchTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 0, '.......T')
                  .with_kmer('AAT', 0, 'a....C..')
                  .with_kmer('ATC', 0, 'a.......')
                  .with_start_kmer_string('GAT'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('GAT', 'ATT', 'TTT')
         .has_n_edges(2))

    def test_four_kmers_one_revcomp_returns_graph_with_four_kmers_as_revcomp(self):
        # given
        driver = (BranchTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 0, '.......T')
                  .with_kmer('AAT', 0, 'a....C..')
                  .with_kmer('ATC', 0, 'a.....G.')
                  .with_kmer('CGA', 0, '.......T')
                  .with_start_kmer_string('CGA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_last_kmer_string('TTT')
         .has_nodes('CGA', 'GAT', 'ATT', 'TTT')
         .has_n_edges(3))

    def test_cycle_is_traversed_once(self):
        # given
        driver = (BranchTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('CAA', 0, '....A...')
                  .with_kmer('AAA', 0, '.c.t...T')
                  .with_kmer('AAT', 0, 'a...A...')
                  .with_kmer('ATA', 0, 'a...A...')
                  .with_kmer('TAA', 0, 'a...A...')
                  .with_start_kmer_string('CAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_last_kmer_string('TAA')
         .has_nodes('CAA', 'AAA', 'AAT', 'ATA', 'TAA')
         .has_n_edges(4))

    def test_one_seen_kmer_returns_empty_graph(self):
        # given
        driver = (BranchTestDriver()
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
        driver = (BranchTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 0, '.......T')
                  .with_kmer('AAT', 0, 'a.......')
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

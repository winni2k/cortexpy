import attr
import pytest

import cortexpy.graph.parser
import cortexpy.test.builder as builder
import cortexpy.test.expectation as expectation


@attr.s(slots=True)
class EngineTestDriver(object):
    graph_builder = attr.ib(attr.Factory(builder.Graph))
    start_kmer_string = attr.ib(None)
    max_nodes = attr.ib(1000)

    def with_kmer(self, *args):
        self.graph_builder.with_kmer(*args)
        return self

    def with_kmer_size(self, n):
        self.graph_builder.with_kmer_size(n)
        return self

    def with_start_kmer_string(self, start_kmer_string):
        self.start_kmer_string = start_kmer_string
        return self

    def with_max_nodes(self, max_nodes):
        self.max_nodes = max_nodes
        return self

    def run(self):
        random_access_parser = cortexpy.graph.parser.RandomAccess(self.graph_builder.build())
        graph = (cortexpy.graph.traversal
                 .Engine(random_access_parser, max_nodes=self.max_nodes)
                 .traverse_from(self.start_kmer_string))
        return expectation.graph.KmerGraphExpectation(graph)


class Test(object):
    def test_raises_on_empty(self):
        # given
        driver = (EngineTestDriver()
                  .with_kmer_size(3)
                  .with_start_kmer_string('AAA'))

        # when/then
        with pytest.raises(KeyError):
            driver.run()

    def test_three_connected_kmers_returns_graph_with_three_kmers(self):
        # given
        driver = (EngineTestDriver()
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

    def test_four_connected_kmers_in_star_returns_graph_with_four_kmers(self):
        # given
        driver = (EngineTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 0, '.......T')
                  .with_kmer('AAT', 0, 'a....CG.')
                  .with_kmer('ATC', 0, 'a.......')
                  .with_kmer('ATG', 0, 'a.......')
                  .with_start_kmer_string('AAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('AAA', 'AAT', 'ATC', 'ATG')
         .has_n_edges(3))

    def test_cycle_is_traversed_once(self):
        # given
        driver = (EngineTestDriver()
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
         .has_nodes('CAA', 'AAA', 'AAT', 'ATA', 'TAA')
         .has_n_edges(5))

    def test_cycle_and_branch_is_traversed_once(self):
        # given
        driver = (EngineTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('CAA', 0, '....A...')
                  .with_kmer('AAA', 0, '.c.t.C.T')
                  .with_kmer('AAC', 0, 'a.......')
                  .with_kmer('AAT', 0, 'a...A...')
                  .with_kmer('ATA', 0, 'a...A...')
                  .with_kmer('TAA', 0, 'a...A...')
                  .with_start_kmer_string('CAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('CAA', 'AAA', 'AAT', 'ATA', 'TAA', 'AAC')
         .has_n_edges(6))

    def test_two_cycles_are_traversed_once(self):
        # given
        driver = (EngineTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('CAA', 0, 'a...A...')
                  .with_kmer('AAA', 0, '.c.t.C.T')
                  .with_kmer('AAC', 0, 'a...A...')
                  .with_kmer('ACA', 0, 'a...A...')
                  .with_kmer('AAT', 0, 'a...A...')
                  .with_kmer('ATA', 0, 'a...A...')
                  .with_kmer('TAA', 0, 'a...A...')
                  .with_start_kmer_string('CAA'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('CAA', 'AAA', 'AAT', 'ATA', 'TAA', 'AAC', 'ACA')
         .has_n_edges(8))

    def test_two_cycles_are_traversed_once_in_revcomp(self):
        # given
        driver = (EngineTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('CAA', 0, 'a...A...')
                  .with_kmer('AAA', 0, '.c.t.C.T')
                  .with_kmer('AAC', 0, 'a...A...')
                  .with_kmer('ACA', 0, 'a...A...')
                  .with_kmer('AAT', 0, 'a...A...')
                  .with_kmer('ATA', 0, 'a...A...')
                  .with_kmer('TAA', 0, 'a...A...')
                  .with_start_kmer_string('TTG'))

        # when
        expect = driver.run()

        # then
        (expect
         .has_nodes('TTG', 'TTT', 'ATT', 'TAT', 'TTA', 'GTT', 'TGT')
         .has_n_edges(8))


class TestMaxNodes(object):
    def test_of_two_raises_with_three_connected_kmers(self):
        # given
        driver = (EngineTestDriver()
                  .with_kmer_size(3)
                  .with_max_nodes(2)
                  .with_kmer('AAA', 0, '.......T')
                  .with_kmer('AAT', 0, 'a....C..')
                  .with_kmer('ATC', 0, 'a.......')
                  .with_start_kmer_string('AAA'))

        # when/then
        with pytest.raises(OverflowError):
            driver.run()

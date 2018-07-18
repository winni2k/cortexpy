import pytest

from cortexpy.test.driver.graph.serializer import CollapseKmerUnitigsTestDriver


class TestCreatesSingleUnitig(object):
    def test_with_missing_kmer(self):
        # given
        driver = CollapseKmerUnitigsTestDriver() \
            .with_kmer_size(3) \
            .retrieve_contig('GTT')

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(1)
        expect.has_one_node_with_repr('GTT').has_coverages([0, 1])

    def test_with_one_kmer(self):
        # given
        driver = CollapseKmerUnitigsTestDriver() \
            .with_kmer_size(3) \
            .with_kmer('AAC', 1) \
            .retrieve_contig('GTT')

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(1)
        expect.has_one_node_with_repr('GTT').has_coverages([1, 1])

    def test_with_two_linked_kmers(self):
        driver = CollapseKmerUnitigsTestDriver() \
            .with_kmer_size(3) \
            .with_kmer('AAA 0 .....C..') \
            .with_kmer('AAC 0 a.......') \
            .retrieve_contig('AAAC')

        # when
        expect = driver.run()

        # then
        assert expect.has_kmers('AAAC')
        expect.has_n_nodes(1).has_one_node_with_repr('AAAC').has_coverages([0, 1], [0, 1])

    def test_with_two_linked_kmers_retrieving_revcomp(self):
        driver = CollapseKmerUnitigsTestDriver() \
            .with_kmer_size(3) \
            .with_kmer('AAA 0 .....C..') \
            .with_kmer('AAC 0 a.......') \
            .retrieve_contig('GTTT')

        # when
        expect = driver.run()

        # then
        assert expect.has_kmers('GTTT')
        expect.has_n_nodes(1).has_one_node_with_repr('GTTT').has_coverages([0, 1], [0, 1])

    def test_with_three_linked_kmers(self):
        driver = CollapseKmerUnitigsTestDriver() \
            .with_kmer_size(3) \
            .with_kmer('AAA', 0, '.....C..') \
            .with_kmer('AAC', 0, 'a....C..') \
            .with_kmer('ACC', 0, 'a.......') \
            .retrieve_contig('AAACC')

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(1)
        expect.has_one_node_with_repr('AAACC')

    def test_with_two_unlinked_missing_kmers(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('CAA')
                  .with_kmer('AAC')
                  .retrieve_contig('GTTG'))

        # when
        expect = driver.run()

        # then
        expect.has_one_node_with_repr('GTTG').has_coverages([0, 1], [0, 1])
        expect.has_n_nodes(1)
        expect.has_n_edges(0)


class TestTraverse(object):
    def test_with_missing_kmer(self):
        # given
        driver = CollapseKmerUnitigsTestDriver() \
            .with_kmer_size(3) \
            .traverse_with_start_kmer_and_colors('GTT', 0)

        # when/then
        with pytest.raises(KeyError):
            driver.run()

    def test_with_one_kmer(self):
        # given
        driver = CollapseKmerUnitigsTestDriver() \
            .with_kmer_size(3) \
            .with_kmer('AAC', 1) \
            .traverse_with_start_kmer_and_colors('GTT', 0)

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(1)
        expect.has_one_node_with_repr('GTT').has_coverages([1])

    def test_with_two_linked_kmers(self):
        driver = CollapseKmerUnitigsTestDriver() \
            .with_kmer_size(3) \
            .with_kmer('AAA 0 .....C..') \
            .with_kmer('AAC 0 a.......') \
            .traverse_with_start_kmer_and_colors('AAA', 0)

        # when
        expect = driver.run()

        # then
        assert expect.has_kmers('AAAC')
        expect.has_n_nodes(1).has_one_node_with_repr('AAAC').has_coverages([0], [0])

    def test_with_two_linked_kmers_retrieving_revcomp(self):
        driver = CollapseKmerUnitigsTestDriver() \
            .with_kmer_size(3) \
            .with_kmer('AAA 0 .....C..') \
            .with_kmer('AAC 0 a.......') \
            .traverse_with_start_kmer_and_colors('GTT', 0)

        # when
        expect = driver.run()

        # then
        assert expect.has_kmers('GTTT')
        expect.has_n_nodes(1).has_one_node_with_repr('GTTT').has_coverages([0], [0])

    def test_with_three_linked_kmers(self):
        driver = CollapseKmerUnitigsTestDriver() \
            .with_kmer_size(3) \
            .with_kmer('AAA', 0, '.....C..') \
            .with_kmer('AAC', 0, 'a....C..') \
            .with_kmer('ACC', 0, 'a.......') \
            .traverse_with_start_kmer_and_colors('AAA', 0)

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(1)
        expect.has_one_node_with_repr('AAACC')

    def test_with_two_unlinked_missing_kmers(self):
        # given
        driver = CollapseKmerUnitigsTestDriver() \
            .with_kmer_size(3) \
            .with_kmer('CAA') \
            .with_kmer('AAC') \
            .traverse_with_start_kmer_and_colors('CAA', 0)

        # when
        expect = driver.run()

        # then
        expect.has_one_node_with_repr('CAA').has_coverages([0], [0])
        expect.has_n_nodes(1)
        expect.has_n_edges(0)


class TestCreatesMultipleUnitigs(object):
    def test_with_two_unlinked_kmers_creates_two_unitigs(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('CAA', 1)
                  .with_kmer('AAC', 1)
                  .retrieve_contig('GTTG'))

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(2)
        expect.has_one_node_with_repr('GTT').has_coverages([1, 1])
        expect.has_one_node_with_repr('G').has_coverages([1, 1])
        expect.has_n_edges(1)

    def test_with_two_node_path_and_three_node_cycle_results_in_two_unitigs(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 1, '.....C..')
                  .with_kmer('AAC', 1, 'a.....G.')
                  .with_kmer('ACG', 1, 'a.g.A...')
                  .with_kmer('CGA', 1, 'a....C..')
                  .with_kmer('GAC', 1, '.c....G.')
                  .retrieve_contig('AAACGAC'))

        # when
        expect = driver.run()

        # then
        expect.has_kmers('AAAC', 'GAC')
        expect.has_one_node_with_repr('AAAC').has_coverages([1, 1], [1, 1])
        expect.has_one_node_with_repr('GAC').has_coverages([1, 1], [1, 1], [1, 1])
        expect.has_n_edges(3)
        expect.has_n_missing_edges(0)

    def test_four_node_path_with_one_node_bubble_in_three_nodes_to_three_unitigs(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAC', 1, '.....C..')
                  .with_kmer('ACC', 1, 'a....CG.')
                  .with_kmer('CCC', 1, 'a.....G.')
                  .with_kmer('CCG', 1, 'ac..A...')
                  .with_kmer('CGA', 1, '.c......')
                  .retrieve_contig('AACCCGA'))

        # when
        expect = driver.run()

        # then
        expect.has_one_node_with_repr('AACC').has_coverages([1, 1], [1, 1])
        expect.has_one_node_with_repr('C').has_coverages([1, 1])
        expect.has_one_node_with_repr('GA').has_coverages([1, 1])
        expect.has_n_nodes(3)

    def test_two_kmers_one_kmer_apart_do_not_collapse(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('GAA', 1, '........')
                  .with_kmer('ACC', 1, '........')
                  .retrieve_contig('GGTTC'))

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(3)
        expect.has_n_edges(2)
        expect.has_one_node_with_repr('GGT').has_coverages([1, 1])
        expect.has_one_node_with_repr('T').has_coverages([0, 1])
        expect.has_one_node_with_repr('C').has_coverages([1, 1])

    def test_two_linked_kmers_with_incoming_edge_and_missing_kmer_returns_three_unitigs(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 1, '.c...C..')
                  .with_kmer('AAC', 1, 'a.......')
                  .retrieve_contig('GTTTA'))

        # when
        expect = driver.run()

        # then
        expect.has_one_node_with_repr('A').has_coverages([0, 1])
        expect.has_one_node_with_repr('G').has_coverages([0, 0])
        expect.has_one_node_with_repr('GTTT').has_coverages([1, 1], [1, 1])

        expect.has_n_nodes(3)
        expect.has_n_edges(2)

    def test_unlinked_kmers_followed_by_two_linked_kmers_collapse_to_two_unitigs(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 1, '.....C..')
                  .with_kmer('AAC', 1, 'a.......')
                  .retrieve_contig('GTTTAA'))

        # when
        expect = driver.run()

        # then
        expect.has_one_node_with_repr('GTTT').has_coverages([1, 1], [1, 1])
        expect.has_one_node_with_repr('AA').has_coverages([0, 2], [0, 2])
        expect.has_n_nodes(2)
        expect.has_n_edges(1)

    def test_linked_kmers_with_outgoing_edge_surrounded_by_missing_kmers_returns_four_unitigs(self):
        # given
        driver = (CollapseKmerUnitigsTestDriver()
                  .with_kmer_size(3)
                  .with_kmer('AAA', 1, '.....C..')
                  .with_kmer('AAC', 1, 'a....C..')
                  .retrieve_contig('CAAGTTTAGG'))

        # when
        expect = driver.run()

        # then
        expect.has_one_node_with_repr('TT').has_coverages([1, 1], [1, 1])
        expect.has_one_node_with_repr('CAAGT').has_coverages([0, 1], [0, 1], [0, 1])
        expect.has_one_node_with_repr('GGT').has_coverages([0, 0])
        expect.has_one_node_with_repr('AGG').has_coverages([0, 1])
        expect.has_n_nodes(4)
        expect.has_n_edges(3)

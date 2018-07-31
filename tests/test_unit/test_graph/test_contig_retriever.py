from cortexpy.graph.contig_retriever import ContigRetriever
import cortexpy.test.builder as builder
from cortexpy.test.expectation.graph import KmerGraphExpectation


class TestGetKmers(object):
    def test_two_nodes_linking_to_self(self):
        # given
        graph_builder = builder.Graph().with_kmer_size(3)

        # when
        kmer_list = ContigRetriever(graph_builder.build()).get_kmers('TTAA')

        # then
        assert len(kmer_list) == 2
        kmer = kmer_list[0][0]
        assert kmer_list[0][0].kmer == 'TAA'
        assert kmer_list[0][1] == 'TTA'
        assert kmer_list[1][0].kmer == 'TAA'
        assert kmer_list[1][1] == 'TAA'
        assert kmer_list[1][0] is kmer

        assert kmer.edges[1].is_edge('t')
        for letter in 'acgACGT':
            assert not kmer.edges[1].is_edge(letter)


class TestGetKmerGraph(object):
    def test_with_no_kmer_returns_missing_kmer(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        retriever = ContigRetriever(graph_builder.build())

        # when
        expect = KmerGraphExpectation(retriever.get_kmer_graph('AAA'))

        # then
        expect.has_n_nodes(1) \
            .has_n_edges(0) \
            .has_node('AAA') \
            .has_coverages(0, 1)

    def test_with_one_kmer_returns_one_kmer(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '........')
        retriever = ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('AAA')

        # then
        assert len(kmer_graph.edges) == 0
        assert list(kmer_graph) == ['AAA']

    def test_with_one_kmer_asking_for_longer_contig_returns_one_kmer_with_coverage_2(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '........')
        retriever = ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('AAAA')

        # then
        assert 1 == len(kmer_graph.edges)
        assert list(kmer_graph) == ['AAA']
        assert [1, 2] == list(kmer_graph.nodes['AAA']['kmer'].coverage)

    def test_with_two_linked_kmers_returns_two_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '.....C..')
        graph_builder.with_kmer('AAC', 1, 'a.......')
        retriever = ContigRetriever(graph_builder.build())

        # when
        expect = KmerGraphExpectation(retriever.get_kmer_graph('AAA'))

        # then
        expect.has_nodes('AAA', 'AAC').has_n_edges(1)
        expect.has_edge('AAA', 'AAC', 0)

    def test_with_three_linked_kmers_and_two_colors_returns_three_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3)
                         .with_num_colors(2))
        graph_builder.with_kmer('AAA', [1, 1], ['.....C..', '.......T'])
        graph_builder.with_kmer('AAC', [1, 0], ['a.......', '........'])
        graph_builder.with_kmer('AAT', [0, 1], ['........', 'a.......'])
        retriever = ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('AAAC')

        # then
        assert set(kmer_graph) == {'AAA', 'AAC', 'AAT'}
        assert set(kmer_graph.edges) == {('AAA', 'AAC', 0), ('AAA', 'AAT', 1), ('AAA', 'AAC', 2)}

    def test_with_two_unlinked_kmers_request_for_three_returns_three_linked_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '........')
        graph_builder.with_kmer('ACC', 1, '........')
        retriever = ContigRetriever(graph_builder.build())

        # when
        expect = KmerGraphExpectation(retriever.get_kmer_graph('AAACC'))

        # then
        expect.has_nodes('AAA', 'AAC', 'ACC')
        expect.has_n_edges(2)
        expect.has_edge('AAA', 'AAC', 1)
        expect.has_edge('AAC', 'ACC', 1)

    def test_with_two_neighboring_unlinked_kmers_returns_two_kmers_linked_by_missing_edge(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '........')
        graph_builder.with_kmer('AAC', 1, '........')
        retriever = ContigRetriever(graph_builder.build())

        # when
        expect = KmerGraphExpectation(retriever.get_kmer_graph('AAAC'))

        # then
        expect.has_nodes('AAA', 'AAC')
        expect.has_n_edges(1)
        expect.has_edge('AAA', 'AAC', 1)


class TestGetKmerGraphComplex(object):
    def test_two_node_path_and_three_node_cycle(self):
        # given
        colors = [0, 1]
        graph_builder = (builder
                         .Graph()
                         .with_kmer_size(3)
                         .with_kmer('AAA', 1, '.....C..')
                         .with_kmer('AAC', 1, 'a.....G.')
                         .with_kmer('ACG', 1, 'a.g.A...')
                         .with_kmer('CGA', 1, 'a....C..')
                         .with_kmer('GAC', 1, '.c....G.')
                         )

        retriever = ContigRetriever(graph_builder.build())

        # when
        expect = KmerGraphExpectation(retriever.get_kmer_graph('AAACGAC'))

        # then
        for color in colors:
            expect.has_edge('AAA', 'AAC', color)
            expect.has_edge('AAC', 'ACG', color)
            expect.has_edge('ACG', 'CGA', color)
            expect.has_edge('CGA', 'GAC', color)
        expect.has_edge('GAC', 'ACG', 0)
        expect.has_n_edges(9)

    def test_two_nodes_linking_to_self(self):
        # given
        graph_builder = builder.Graph().with_kmer_size(3)

        # when
        expect = KmerGraphExpectation(
            ContigRetriever(graph_builder.build()).get_kmer_graph('TTAA'))

        # then
        expect.has_edge('TTA', 'TAA', 1)
        expect.has_n_edges(1)


class TestGetKmerGraphRevcomp(object):
    def test_with_one_kmer_returns_kmer(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '........')
        retriever = ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('TTT')

        # then
        assert len(kmer_graph.edges) == 0
        assert list(kmer_graph.nodes) == ['TTT']

    def test_with_two_linked_kmers_returns_two_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '.....C..')
        graph_builder.with_kmer('AAC', 1, 'a.......')
        retriever = ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('GTTT')

        # then
        assert set(kmer_graph.nodes) == {'GTT', 'TTT'}
        assert set(kmer_graph.edges) == {('GTT', 'TTT', 0), ('GTT', 'TTT', 1)}

    def test_with_three_linked_kmers_and_two_colors_returns_three_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3)
                         .with_num_colors(2))
        graph_builder.with_kmer('AAA', [1, 1], ['.....C..', '.......T'])
        graph_builder.with_kmer('AAC', [1, 0], ['a.......', '........'])
        graph_builder.with_kmer('AAT', [0, 1], ['........', 'a.......'])
        retriever = ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('GTTT')

        # then
        assert set(kmer_graph.nodes) == {'TTT', 'GTT', 'ATT'}
        assert set(kmer_graph.edges) == {('GTT', 'TTT', 0),
                                         ('GTT', 'TTT', 2),
                                         ('ATT', 'TTT', 1)}

        assert list(kmer_graph.nodes['TTT']['kmer'].coverage) == [1, 1, 1]
        assert list(kmer_graph.nodes['GTT']['kmer'].coverage) == [1, 0, 1]

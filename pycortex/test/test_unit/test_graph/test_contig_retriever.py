import pycortex.test.builder as builder
import pycortex.graph as graph


class TestGetKmerGraph(object):
    def test_with_no_kmer_returns_missing_kmer(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        retriever = graph.ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('AAA')

        # then
        assert len(kmer_graph.edges) == 0
        assert list(kmer_graph.nodes) == ['AAA']

    def test_with_one_kmer_returns_one_kmer(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '........')
        retriever = graph.ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('AAA')

        # then
        assert len(kmer_graph.edges) == 0
        assert list(kmer_graph.nodes) == ['AAA']

    def test_with_two_linked_kmers_returns_two_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '.....C..')
        graph_builder.with_kmer('AAC', 1, 'a.......')
        retriever = graph.ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('AAAC')

        # then
        assert set(kmer_graph.nodes) == {'AAA', 'AAC'}
        assert list(kmer_graph.edges) == [('AAA', 'AAC')]

    def test_with_three_linked_kmers_and_two_colors_returns_three_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3)
                         .with_num_colors(2))
        graph_builder.with_kmer('AAA', [1, 1], ['.....C..', '.......T'])
        graph_builder.with_kmer('AAC', [1, 0], ['a.......', '........'])
        graph_builder.with_kmer('AAT', [0, 1], ['........', 'a.......'])
        retriever = graph.ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('AAAC')

        # then
        assert set(kmer_graph.nodes) == {'AAA', 'AAC', 'AAT'}
        assert set(kmer_graph.edges) == {('AAA', 'AAC'), ('AAA', 'AAT')}


class TestGetKmerGraphRevcomp(object):
    def test_with_one_kmer_returns_kmer(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '........')
        retriever = graph.ContigRetriever(graph_builder.build())

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
        retriever = graph.ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('GTTT')

        # then
        assert set(kmer_graph.nodes) == {'GTT', 'TTT'}
        assert list(kmer_graph.edges) == [('GTT', 'TTT')]

    def test_with_three_linked_kmers_and_two_colors_returns_three_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3)
                         .with_num_colors(2))
        graph_builder.with_kmer('AAA', [1, 1], ['.....C..', '.......T'])
        graph_builder.with_kmer('AAC', [1, 0], ['a.......', '........'])
        graph_builder.with_kmer('AAT', [0, 1], ['........', 'a.......'])
        retriever = graph.ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('GTTT')

        # then
        assert set(kmer_graph.nodes) == {'TTT', 'GTT', 'ATT'}
        assert set(kmer_graph.edges) == {('GTT', 'TTT'), ('ATT', 'TTT')}

        assert list(kmer_graph.nodes['TTT']['kmer'].coverage) == [1, 1]
        assert list(kmer_graph.nodes['GTT']['kmer'].coverage) == [1, 0]

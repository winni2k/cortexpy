import pycortex.test.builders as builders

import pycortex.graph as graph


class TestGetKmerGraph(object):
    def test_with_one_kmer_returns_one_kmer(self):
        # given
        builder = (builders.graph.Graph()
                   .with_kmer_size(3))
        builder.with_kmer('AAA', 1, '........')
        retriever = graph.ContigRetriever(builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('AAA')

        # then
        assert len(kmer_graph.edges) == 0
        assert list(kmer_graph.nodes) == ['AAA']

    def test_with_two_linked_kmers_returns_two_kmers(self):
        # given
        builder = (builders.graph.Graph()
                   .with_kmer_size(3))
        builder.with_kmer('AAA', 1, '.....C..')
        builder.with_kmer('AAC', 1, 'a.......')
        retriever = graph.ContigRetriever(builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('AAAC')

        # then
        assert set(kmer_graph.nodes) == {'AAA', 'AAC'}
        assert list(kmer_graph.edges) == [('AAA', 'AAC')]

    def test_with_three_linked_kmers_and_two_colors_returns_three_kmers(self):
        # given
        builder = (builders.graph.Graph()
                   .with_kmer_size(3)
                   .with_num_colors(2))
        builder.with_kmer('AAA', [1, 1], ['.....C..', '.......T'])
        builder.with_kmer('AAC', [1, 0], ['a.......', '........'])
        builder.with_kmer('AAT', [0, 1], ['........', 'a.......'])
        retriever = graph.ContigRetriever(builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('AAAC')

        # then
        assert set(kmer_graph.nodes) == {'AAA', 'AAC', 'AAT'}
        assert set(kmer_graph.edges) == {('AAA', 'AAC'), ('AAA', 'AAT')}

    def test_with_one_kmer_returns_revcomp_kmer_if_asked(self):
        # given
        builder = (builders.graph.Graph()
                   .with_kmer_size(3))
        builder.with_kmer('AAA', 1, '........')
        retriever = graph.ContigRetriever(builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('TTT')

        # then
        assert len(kmer_graph.edges) == 0
        assert list(kmer_graph.nodes) == ['TTT']

    def test_with_two_linked_kmers_returns_two_kmers_revcomp(self):
        # given
        builder = (builders.graph.Graph()
                   .with_kmer_size(3))
        builder.with_kmer('AAA', 1, '.....C..')
        builder.with_kmer('AAC', 1, 'a.......')
        retriever = graph.ContigRetriever(builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('GTTT')

        # then
        assert set(kmer_graph.nodes) == {'GTT', 'TTT'}
        assert list(kmer_graph.edges) == [('GTT', 'TTT')]

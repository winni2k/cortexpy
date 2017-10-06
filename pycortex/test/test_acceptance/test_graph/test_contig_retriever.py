from pycortex import graph
from pycortex.test import builder


class TestGetKmerGraph(object):
    def test_with_three_linked_kmers_and_two_colors_returns_three_kmers(self, tmpdir):
        # given
        kmer_size = 3
        output_graph = (builder.Mccortex(kmer_size)
                        .with_dna_sequence('AAAC')
                        .with_dna_sequence('AAAT')
                        .build(tmpdir))

        retriever = graph.ContigRetriever(open(output_graph, 'rb'))

        # when
        kmer_graph = retriever.get_kmer_graph('GTTT')

        # then
        assert set(kmer_graph.nodes) == {'TTT', 'GTT', 'ATT'}
        assert set(kmer_graph.edges) == {('GTT', 'TTT'), ('ATT', 'TTT')}

import pytest

from cortexpy.graph.contig_retriever import ContigRetriever
from cortexpy.test import builder


class TestGetKmerGraph:
    @pytest.mark.parametrize('builder_class,contig_retriever',
                             [(builder.Mccortex, ContigRetriever.from_cortex), ])
    # (builder.Bifrost, ContigRetriever.from_bifrost)])
    def test_with_three_linked_kmers_and_two_colors_returns_three_kmers(self, builder_class,
                                                                        contig_retriever, tmpdir):
        # given
        kmer_size = 3
        output_graph = (builder_class(kmer_size)
                        .with_dna_sequence('AAAC')
                        .with_dna_sequence('AAAT')
                        .build(tmpdir))

        retriever = contig_retriever(open(output_graph, 'rb'))

        # when
        kmer_graph = retriever.get_kmer_graph('GTTT')

        # then
        assert set(kmer_graph.nodes) == {'TTT', 'GTT', 'ATT'}
        assert set(kmer_graph.edges) == {
            ('GTT', 'TTT', 0),
            ('GTT', 'TTT', 1),
            ('ATT', 'TTT', 0)
        }

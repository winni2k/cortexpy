import pycortex.test.builders as builders

import pycortex.graph as graph


class TestContigRetriever(object):
    def test_retrieves_one_kmer_from_graph_with_one_kmer(self):
        # given
        builder = (builders.graph.Graph()
                   .with_kmer_size(3))
        builder.with_kmer('AAA', 1, '........')
        retriever = graph.ContigRetriever(builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('AAA')

        # then
        assert len(kmer_graph.edges) == 0
        assert len(kmer_graph.nodes) == 1
        assert 'AAA' in kmer_graph.nodes

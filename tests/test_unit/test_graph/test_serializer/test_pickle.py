import io

import networkx as nx

from cortexpy.graph.contig_retriever import ContigRetriever
from cortexpy.test import builder as builder


class Test(object):
    def test_two_linked_kmers_pickle_ok(self):
        # given
        color_names = 'samp1', 'samp2'
        graph_builder = builder.Graph() \
            .with_kmer_size(3) \
            .with_num_colors(2) \
            .with_color_names(*color_names) \
            .with_kmer('AAA', [1, 1], ['.....C..', '.......T']) \
            .with_kmer('AAC', [1, 0], ['a.......', '........'])
        retriever = ContigRetriever(graph_builder.build())
        kmer_graph = retriever.get_kmer_graph('GTTT')

        # when
        buffer = io.BytesIO()
        nx.write_gpickle(kmer_graph, buffer)
        buffer.seek(0)
        unpickled_kmer_graph = nx.read_gpickle(buffer)

        # then
        assert len(unpickled_kmer_graph) == len(kmer_graph)
        unpickle_node_data = unpickled_kmer_graph.nodes(data=True)
        for node, data in kmer_graph.nodes(data=True):
            assert unpickle_node_data[node] == data

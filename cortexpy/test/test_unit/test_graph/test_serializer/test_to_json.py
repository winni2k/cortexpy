import json

import cortexpy.graph.serializer.serializer
from cortexpy.graph.interactor import Interactor
from cortexpy.graph.contig_retriever import ContigRetriever
from cortexpy.graph.parser.streaming import load_cortex_graph
from cortexpy.test import builder as builder, expectation


class TestFromKmerGraph(object):
    def test_two_linked_kmers_are_jsonifiable(self):
        # given
        color_names = ['samp1', 'samp2']
        graph_builder = builder.Graph() \
            .with_kmer_size(3) \
            .with_num_colors(2) \
            .with_color_names(*color_names) \
            .with_kmer('AAA 1 1 .....C.. ........') \
            .with_kmer('AAC 1 0 a....... ........')

        retriever = ContigRetriever(graph_builder.build())
        graph = retriever.get_kmer_graph('GTTT')

        # when
        kmer_json = cortexpy.graph.serializer.serializer.Serializer(graph).to_json()
        expect = expectation.JsonGraph.from_string(kmer_json)

        # then
        kmer_data = json.loads(kmer_json)  # does not raise
        assert kmer_data['graph']['colors'] == [0, 1, 2]
        assert kmer_data['graph']['sample_names'] == color_names + ['retrieved_contig']

        expect.has_n_nodes(2)
        expect.has_n_edges(2)


class TestFromCortexGraph(object):
    def test_two_linked_kmers_are_jsonifiable(self):
        # given
        colors = (0, 1)
        color_names = ['samp1', 'samp2']
        graph_builder = builder.Graph() \
            .with_kmer_size(3) \
            .with_num_colors(2) \
            .with_color_names(*color_names) \
            .with_kmer('AAA 1 1 .....C.. ........') \
            .with_kmer('AAC 1 0 a....... ........')

        graph = load_cortex_graph(graph_builder.build())
        graph = Interactor(graph) \
            .make_graph_nodes_consistent(seed_kmer_strings=['GTT']) \
            .graph
        kmer_json = cortexpy.graph.serializer.serializer.Serializer(graph).to_json()

        # when
        expect = expectation.JsonGraph.from_string(kmer_json)

        # then
        kmer_data = json.loads(kmer_json)  # does not raise
        assert kmer_data['graph']['colors'] == list(colors)
        assert kmer_data['graph']['sample_names'] == color_names

        expect.has_n_nodes(2)
        expect.has_n_edges(1)

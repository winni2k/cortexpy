import pytest

from cortexpy.graph.interactor import Interactor
from cortexpy.test.builder.graph.cortex import get_cortex_builder


class Test(object):
    @pytest.mark.parametrize('seed', ['AAA', 'TTT'])
    def test_single_kmer_revcomp_seed(self, seed):
        # given
        b = get_cortex_builder()
        b.with_kmer('AAA 1 ......G.')
        b.with_kmer('AAG 1 a.......')
        cdb = b.build()

        # when
        graph = Interactor(cdb).make_graph_nodes_consistent([seed]).graph

        # then
        if seed == 'AAA':
            assert [] == list(graph.in_edges(seed))
            assert [('AAA', 'AAG')] == list(graph.out_edges(seed))
        else:
            assert [('CTT', 'TTT')] == list(graph.in_edges(seed))
            assert [] == list(graph.out_edges(seed))

    def test_gets_correct_neighbors_of_kmer(self):
        # given
        b = get_cortex_builder()
        b.with_kmer('AAC 1 .......T')
        b.with_kmer('ACT 1 a.....G.')
        b.with_kmer('CAG 1 .......T')
        cdb = b.build()
        seed = 'AAC'

        # when
        graph = Interactor(cdb).make_graph_nodes_consistent([seed]).graph

        # then
        assert ['CTG'] == list(graph['ACT'])
        assert ['CTG'] == list(graph.succ['ACT'])
        assert ['AAC'] == list(graph.pred['ACT'])
        assert [('ACT', 'CTG')] == list(graph.out_edges('ACT'))
        assert [('AAC', 'ACT')] == list(graph.in_edges('ACT'))

        assert [] == list(graph['CTG'])
        assert [] == list(graph.succ['CTG'])
        assert ['ACT'] == list(graph.pred['CTG'])
        assert [] == list(graph.out_edges('CTG'))
        assert [('ACT', 'CTG')] == list(graph.in_edges('CTG'))

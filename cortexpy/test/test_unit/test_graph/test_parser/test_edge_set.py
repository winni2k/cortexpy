import numpy as np
import pytest

from cortexpy.edge_set import EdgeSet


class TestIsEdge(object):
    def test_with_all_true(self):
        es = EdgeSet(np.ones(8))
        for letter in 'acgtACGT':
            assert es.is_edge(letter)

    def test_with_none_true(self):
        es = EdgeSet(np.zeros(8))
        for letter in 'acgtACGT':
            assert not es.is_edge(letter)


class TestAddEdge(object):
    def test_adds_each_edge(self):
        es = EdgeSet(np.zeros(8))
        for letter in 'acgtACGT':
            assert not es.is_edge(letter)
            es.add_edge(letter)
            assert es.is_edge(letter)


class TestRemoveEdge(object):
    def test_removes_each_edge(self):
        es = EdgeSet(np.ones(8))
        for letter in 'acgtACGT':
            assert es.is_edge(letter)
            es.remove_edge(letter)
            assert not es.is_edge(letter)


class TestGetitem(object):
    def test_works(self):
        es = EdgeSet(np.ones(8))
        for edge_idx in range(8):
            assert es[edge_idx]


class TestIncomingOutgoingEdges(object):
    def test_with_all_incoming_and_no_outgiong(self):
        es = EdgeSet(np.concatenate([np.ones(4), np.zeros(4)]))
        assert list(es.incoming) == [1, 1, 1, 1]
        assert list(es.outgoing) == [0, 0, 0, 0]


class TestIncomingOutgoingKmers(object):
    def test_no_incoming_or_outgoing(self):
        es = EdgeSet(np.zeros(8))
        assert 0 == len(es.get_incoming_kmers('AAA'))
        assert 0 == len(es.get_outgoing_kmers('AAA'))

    def test_all_incoming_no_outgoing(self):
        es = EdgeSet(np.concatenate([np.ones(4), np.zeros(4)]))
        assert es.get_incoming_kmers('AAA') == ['AAA', 'CAA', 'GAA', 'TAA']
        assert 0 == len(es.get_outgoing_kmers('AAA'))

    def test_no_incoming_all_outgoing(self):
        es = EdgeSet(np.concatenate([np.zeros(4), np.ones(4)]))
        assert 0 == len(es.get_incoming_kmers('AAA'))
        assert {'AAA', 'AAC', 'AAG', 'AAT'} == set(es.get_outgoing_kmers('AAA'))

    def test_incoming_returns_lexicographically_lowest_kmers(self):
        es = EdgeSet(np.zeros(8))
        es.add_edge('t')
        assert ['TAA'] == es.get_incoming_kmers('TAA')

    def test_incoming_strings_does_not_return_lexicographically_lowest_kmers(self):
        es = EdgeSet(np.zeros(8))
        es.add_edge('t')
        assert ['TTA'] == list(es.get_incoming_kmer_strings('TAA'))

    def test_outgoing_returns_lexicographically_lowest_kmers(self):
        es = EdgeSet(np.zeros(8))
        es.add_edge('G')
        print(es)
        assert ['CCG'] == es.get_outgoing_kmers('ACG')

    def test_outgoing_strings_does_not_return_lexicographically_lowest_kmer(self):
        es = EdgeSet(np.zeros(8))
        es.add_edge('G')
        assert ['CGG'] == list(es.get_outgoing_kmer_strings('ACG'))

    def test_raises_on_non_lexlo_kmer(self):
        es = EdgeSet(np.zeros(8))
        with pytest.raises(AssertionError):
            es.get_outgoing_kmers('TTT')
        with pytest.raises(AssertionError):
            es.get_incoming_kmers('TTT')


class TestStr(object):
    def test_empty_kmer(self):
        es = EdgeSet(np.zeros(8))
        for as_revcomp in [True, False]:
            assert es.to_str(as_revcomp=as_revcomp) == '........'

    def test_with_a_and_c(self):
        es = EdgeSet(np.zeros(8))
        es.add_edge('A')
        es.add_edge('c')
        assert '.c..A...' == es.to_str()
        assert '...T..g.' == es.to_str(as_revcomp=True)

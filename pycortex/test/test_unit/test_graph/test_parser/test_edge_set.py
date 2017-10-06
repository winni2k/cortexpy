import numpy as np

from pycortex.edge_set import EdgeSet


class TestIsEdge(object):
    def test_with_all_true(self):
        es = EdgeSet(np.ones(8))
        for letter in 'acgtACGT':
            assert es.is_edge(letter)

    def test_with_none_true(self):
        es = EdgeSet(np.zeros(8))
        for letter in 'acgtACGT':
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
        assert len(es.get_incoming_kmers('AAA')) == 0
        assert len(es.get_outgoing_kmers('AAA')) == 0

    def test_all_incoming_no_outgoing(self):
        es = EdgeSet(np.concatenate([np.ones(4), np.zeros(4)]))
        assert es.get_incoming_kmers('AAA') == ['AAA', 'CAA', 'GAA', 'TAA']
        assert len(es.get_outgoing_kmers('AAA')) == 0

    def test_no_incoming_all_outgoing(self):
        es = EdgeSet(np.concatenate([np.zeros(4), np.ones(4)]))
        assert len(es.get_incoming_kmers('AAA')) == 0
        assert es.get_outgoing_kmers('AAA') == ['AAA', 'AAC', 'AAG', 'AAT']

import pytest

from cortexpy.test.driver import command


class Test:
    def test_traverses_two_subgraphs_into_three_records(self, tmpdir):
        # given
        d = command.Traverse(tmpdir)
        d.with_records('CCCGC', 'CCCGA', 'AAAT')
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_record_ids('g0_p0', 'g0_p1', 'g0_p2')
        expect.has_n_groups(1)

    def test_traverses_into_two_records_with_custom_graph_idx(self, tmpdir):
        # given
        d = command.Traverse(tmpdir)
        d.with_records('CCCGC', 'CCCGA')
        d.with_graph_index(7)
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_record_ids('g7_p0', 'g7_p1')
        expect.has_n_groups(1)


class TestCycle:
    @pytest.mark.xfail(reason="Due to a bug in how kmers are collapsed to unitigs")
    def test_single_kmer_without_seed_provides_no_output(self, tmpdir):
        # given
        d = command.Traverse(tmpdir)
        d.with_records('CCCC')
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_no_records()

    def test_multiple_kmers_without_seed_provides_no_output(self, tmpdir):
        # given
        d = command.Traverse(tmpdir)
        d.with_records('CCCACCC')
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_no_records()

import pytest

from cortexpy.test.driver import command


class Test:
    def test_traverses_two_subgraphs_into_three_records(self, tmpdir):
        # given
        d = command.TraverseDriver(tmpdir)
        d.with_records('CCCGC', 'CCCGA', 'AAAT')
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_record_ids('g0_p0', 'g0_p1', 'g0_p2')
        expect.has_n_groups(1)

    def test_traverses_into_two_records_with_custom_graph_idx(self, tmpdir):
        # given
        d = command.TraverseDriver(tmpdir)
        d.with_records('CCCGC', 'CCCGA')
        d.with_graph_index(7)
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_record_ids('g7_p0', 'g7_p1')
        expect.has_n_groups(1)


class TestLinks:
    def test_traverses_tangle_into_two_paths(self, tmpdir):
        # given
        d = command.TraverseDriver(tmpdir)
        d.with_kmer_size(5)
        d.with_records('CAAAACCCCC', 'TAAAACCCCT')
        d.with_link_records('CAAAACCCCT', 'TAAAACCCCC')

        # when
        expect = d.run()

        # then
        expect.has_records('CAAAACCCCT', 'TAAAACCCCC')

    def test_traverses_two_bubbles_into_two_paths(self, tmpdir):
        # given
        d = command.TraverseDriver(tmpdir)
        d.with_kmer_size(5)
        d.with_records('AACCACAAAACCCCCGAACA', 'AACCATAAAACCCCTGAACA')
        d.with_link_records('AACCACAAAACCCCTGAACA', 'AACCATAAAACCCCCGAACA')

        # when
        expect = d.run()

        # then
        expect.has_records('AACCACAAAACCCCTGAACA', 'AACCATAAAACCCCCGAACA')

    def test_traverses_tangle_with_non_lexlo_kmer_into_two_paths(self, tmpdir):
        # given
        d = command.TraverseDriver(tmpdir)
        d.with_kmer_size(5)
        d.with_records('CAAATGGGGG', 'TAAATGGGGT')
        d.with_link_records('CAAATGGGGT', 'TAAATGGGGG')

        # when
        expect = d.run()

        # then
        expect.has_records('CAAATGGGGT', 'TAAATGGGGG')


class TestCycle:
    @pytest.mark.xfail(reason="Due to a bug in how kmers are collapsed to unitigs")
    def test_single_kmer_without_seed_provides_no_output(self, tmpdir):
        # given
        d = command.TraverseDriver(tmpdir)
        d.with_records('CCCC')
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_no_records()

    def test_multiple_kmers_without_seed_provides_no_output(self, tmpdir):
        # given
        d = command.TraverseDriver(tmpdir)
        d.with_records('CCCACCC')
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_no_records()

    def test_multiple_kmers_with_extra_start_kmer_reports_pseudo_linear_contig(self, tmpdir):
        # given
        d = command.TraverseDriver(tmpdir)
        d.with_records('CCCACCC')
        d.with_extra_start_kmer('CAC')
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_record('CACCCA')
        expect.has_n_records(1)

    def test_contig_with_extra_start_kmer_reports_shorter_contig(self, tmpdir):
        # given
        d = command.TraverseDriver(tmpdir)
        d.with_records('ACCCAAA')
        d.with_extra_start_kmer('CCC')
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_record('CCCAAA')
        expect.has_record('ACC')
        expect.has_n_records(2)

    def test_contig_with_extra_start_kmer_that_is_start_reports_same_contig(self, tmpdir):
        # given
        d = command.TraverseDriver(tmpdir)
        d.with_records('ACCCAAA')
        d.with_extra_start_kmer('ACC')
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_record('ACCCAAA')
        expect.has_n_records(1)

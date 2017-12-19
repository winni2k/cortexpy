import os

from Bio.Seq import reverse_complement

from cortexpy.test import builder, runner, expectation

if os.environ.get('CI'):
    SPAWN_PROCESS = True
else:
    SPAWN_PROCESS = False


class Test(object):
    def test_outputs_fasta(self, tmpdir):
        # given
        record = 'ATTCC'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        # when
        completed_process = (runner
                             .Cortexpy(SPAWN_PROCESS)
                             .view_traversal(output_format='fasta', graph=output_graph,
                                             contig=record))
        stdout = completed_process.stdout

        # then
        assert completed_process.returncode == 0, completed_process
        expect = expectation.Fasta(stdout)
        (expect
         .has_record('ATT')
         .has_record('TTC')
         .has_record('TCC')
         .has_n_records(3))


class TestContigs(object):
    def test_outputs_multiple_combinations(self, tmpdir):
        # given
        records = [
            'CAACC',
            'AAACA',
            'AAACT',
        ]
        kmer_size = 3
        maker = builder.Mccortex().with_kmer_size(kmer_size)
        for rec in records:
            maker.with_dna_sequence(rec)

        output_graph = maker.build(tmpdir)

        # when
        completed_process = (runner
                             .Cortexpy(SPAWN_PROCESS)
                             .view_traversal(output_format='fasta', graph=output_graph,
                                             contig='AAA', output_type='contigs'))
        stdout = completed_process.stdout

        # then
        assert completed_process.returncode == 0, completed_process
        expect = expectation.Fasta(stdout)
        for record in records:
            expect.has_record(record)
        expect.has_record('CAACT')
        expect.has_record('CAACA')
        expect.has_record('AAACC')
        expect.has_n_records(6)

    def test_outputs_contigs_for_bubble(self, tmpdir):
        # given
        records = [
            'AAACCC',
            'AAAGCC',
        ]
        kmer_size = 3
        maker = builder.Mccortex().with_kmer_size(kmer_size)
        for rec in records:
            maker.with_dna_sequence(rec)

        output_graph = maker.build(tmpdir)

        # when
        completed_process = (runner
                             .Cortexpy(SPAWN_PROCESS)
                             .view_traversal(output_format='fasta', graph=output_graph,
                                             contig='AAA', output_type='contigs'))
        stdout = completed_process.stdout

        # then
        assert completed_process.returncode == 0, completed_process
        expect = expectation.Fasta(stdout)
        for record in records:
            expect.has_record(record)
        expect.has_n_records(2)

    def test_dal19_data(self, tmpdir):
        # given
        reverse_complement
        records = [
            'CCCCGAGGGAAGCTCTATGAATTCGCCAATCCCAGTATGCAAAAAATGTTGGAGAGGTATCAAAAGTATTCACAAGAAAGT',
            'GTATGCAAAAAATGTTGGAGAGGTATCAAAAGTATTCACAAGAAAGTGACATA',
            'TAAAAAATGTTGGAGAGGTATCAAAAGTATTCACAAGAAAGTGACATAGATAACACTACCAAAGAGCAAGACTATCAG'
        ]
        expected_records = [records[0] + records[1][47:] + records[2][48:], records[2]]
        kmer_size = 47
        maker = builder.Mccortex().with_kmer_size(kmer_size)
        for rec in records:
            maker.with_dna_sequence(rec)

        output_graph = maker.build(tmpdir)

        # when
        completed_process = \
            (runner.Cortexpy(SPAWN_PROCESS)
             .view_traversal(output_format='fasta', graph=output_graph,
                             contig='CAAAAAATGTTGGAGAGGTATCAAAAGTATTCACAAGAAAGTGACAT',
                             output_type='contigs')
             )
        stdout = completed_process.stdout

        # then
        assert completed_process.returncode == 0, completed_process
        expect = expectation.Fasta(stdout)
        for record in expected_records:
            expect.has_record(record)
        expect.has_n_records(2)

    def test_outputs_warning_on_max_nodes_succeeded(self, tmpdir):
        # given
        records = [
            'CAACC',
        ]
        kmer_size = 3
        maker = builder.Mccortex().with_kmer_size(kmer_size)
        for rec in records:
            maker.with_dna_sequence(rec)

        output_graph = maker.build(tmpdir)

        # when
        stderr = (runner
                  .Cortexpy(True)
                  .view_traversal(graph=output_graph, contig='CAA', max_nodes=1)).stderr

        # then
        assert 'Max nodes (1) exceeded: 3 nodes found' in stderr

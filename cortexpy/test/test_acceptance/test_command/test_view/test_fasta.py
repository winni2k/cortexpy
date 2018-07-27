import os

from cortexpy.test import builder, runner, expectation
from cortexpy.test.driver import command

if os.environ.get('CI'):
    SPAWN_PROCESS = True
else:
    SPAWN_PROCESS = False


class Test(object):
    def test_outputs_fasta_of_kmers(self, tmpdir):
        # given
        record = 'ATTCC'
        kmer_size = 3
        output_graph = builder.Mccortex() \
            .with_dna_sequence(record) \
            .with_kmer_size(kmer_size) \
            .build(tmpdir)

        # when
        completed_process = runner.Cortexpy(SPAWN_PROCESS) \
            .view_traversal(kmers=True, to_json=False, graph=output_graph, contig=record)

        stdout = completed_process.stdout

        # then
        assert completed_process.returncode == 0, completed_process
        expect = expectation.Fasta(stdout)
        expect.has_record('AAT')
        expect.has_record('GAA')
        expect.has_record('GGA')
        expect.has_n_records(3)


class TestContigs(object):
    def test_outputs_multiple_combinations(self, tmpdir):
        # given
        records = [
            'CAACC',
            'AAACA',
            'AAACT',
        ]
        n_paths = 6
        kmer_size = 3
        maker = builder.Mccortex().with_kmer_size(kmer_size)
        for rec in records:
            maker.with_dna_sequence(rec)

        output_graph = maker.build(tmpdir)

        # when
        completed_process = (
            runner.Cortexpy(SPAWN_PROCESS).view_traversal(to_json=False,
                                                          kmers=False,
                                                          graph=output_graph,
                                                          contig='AAA')
        )
        stdout = completed_process.stdout

        # then
        assert completed_process.returncode == 0, completed_process
        expect = expectation.Fasta(stdout)
        ids = ['gx_p{}'.format(i) for i in range(n_paths)]
        for record in records:
            expect.has_record(record).has_id_in(*ids)
        expect.has_record('CAACT').has_id_in(*ids)
        expect.has_record('CAACA').has_id_in(*ids)
        expect.has_record('AAACC').has_id_in(*ids)
        expect.has_n_records(6)

    def test_outputs_contigs_for_bubble(self, tmpdir):
        # given
        records = [
            'AAACCC',
            'AAAGCCC',
        ]
        kmer_size = 3
        maker = builder.Mccortex().with_kmer_size(kmer_size)
        for rec in records:
            maker.with_dna_sequence(rec)

        output_graph = maker.build(tmpdir)

        # when
        completed_process = (
            runner.Cortexpy(SPAWN_PROCESS).view_traversal(output_format='fasta', graph=output_graph,
                                                          contig='AAA', output_type='contigs')
        )
        stdout = completed_process.stdout

        # then
        assert completed_process.returncode == 0, completed_process
        expect = expectation.Fasta(stdout)
        for record in records:
            expect.has_record(record)
        expect.has_n_records(2)

    def test_dal19_data(self, tmpdir):
        # given
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
        completed_process = runner.Cortexpy(SPAWN_PROCESS) \
            .view_traversal(to_json=False,
                            kmers=False,
                            graph=output_graph,
                            contig='CAAAAAATGTTGGAGAGGTATCAAAAGTATTCACAAGAAAGTGACAT')
        stdout = completed_process.stdout

        # then
        assert completed_process.returncode == 0, completed_process
        expect = expectation.Fasta(stdout)
        for record in expected_records:
            expect.has_record(record)
        expect.has_n_records(2)

    def test_raises_on_max_nodes_exceeded(self, tmpdir):
        # given
        query = 'CAA'
        records = ['CAACC']
        kmer_size = 3
        maker = builder.Mccortex().with_kmer_size(kmer_size)
        for rec in records:
            maker.with_dna_sequence(rec)

        output_graph = maker.build(tmpdir)

        # when
        completed_process = runner.Cortexpy(True).view_traversal(graph=output_graph, contig=query,
                                                                 max_nodes=1)

        # then
        assert 0 != completed_process.returncode
        assert 'Max nodes (1) exceeded: 3 nodes found' in completed_process.stderr

    def test_raises_with_kmer(self, tmpdir):
        # given
        query = 'CAACC'
        records = [query]
        kmer_size = 3
        maker = builder.Mccortex().with_kmer_size(kmer_size)
        for rec in records:
            maker.with_dna_sequence(rec)

        output_graph = maker.build(tmpdir)

        # when
        completed_process = runner.Cortexpy(True).traverse(graphs=[output_graph],
                                                           contig=query,
                                                           out='-',
                                                           max_nodes=1)

        # then
        assert 0 != completed_process.returncode
        assert ('Max nodes (1) exceeded: 3 nodes found') in completed_process.stderr

    def test_does_not_raise_without_max_nodes(self, tmpdir):
        # given
        query = 'CAACC'
        records = [query]
        kmer_size = 3
        maker = builder.Mccortex().with_kmer_size(kmer_size)
        for rec in records:
            maker.with_dna_sequence(rec)

        output_graph = maker.build(tmpdir)

        # when
        completed_process = runner.Cortexpy(True).traverse(graphs=[output_graph],
                                                           contig=query,
                                                           out=tmpdir / 'discarded.pickle',
                                                           max_nodes=None)

        # then
        assert 0 == completed_process.returncode


class TestTraversal(object):
    def test_traverses_two_subgraphs_as_two_subgraphs(self, tmpdir):
        # given
        d = command.ViewTraversal(tmpdir)
        d.with_records('CCCGC', 'CCCGA', 'AAAT')
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_n_records(3)
        expect.has_n_groups(1)
        expect.has_record_ids('gx_p0', 'gx_p1', 'gx_p2')

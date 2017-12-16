import os

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
                             .view_traversal(output_type='fasta', graph=output_graph,
                                             contig=record))
        stdout = completed_process.stdout

        # then
        assert completed_process.returncode == 0, completed_process
        expect = expectation.Fasta(stdout)
        (expect.has_n_records(3)
         .has_record('ATT')
         .has_record('TTC')
         .has_record('TCC'))

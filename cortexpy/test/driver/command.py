import attr
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from cortexpy.test import builder, expectation
from cortexpy.test import runner

if os.environ.get('CI'):
    SPAWN_PROCESS = True
else:
    SPAWN_PROCESS = False


@attr.s(slots=True)
class Assemble(object):
    tmpdir = attr.ib()
    records = attr.ib(attr.Factory(list))
    initial_seqs = attr.ib(attr.Factory(list))
    mccortex_builder = attr.ib(attr.Factory(builder.Mccortex))

    def with_records(self, *records):
        self.records = records
        return self

    def with_initial_sequences(self, *seqs):
        self.initial_seqs = seqs
        return self

    def with_kmer_size(self, size):
        self.mccortex_builder.with_kmer_size(size)
        return self

    def run(self):
        inital_seq_fasta = self.tmpdir / 'initial_seqs.fasta'
        initial_seqs = [SeqRecord(Seq(rec), id=str(idx)) for idx, rec in enumerate(self.initial_seqs)]
        SeqIO.write(initial_seqs, str(inital_seq_fasta), "fasta")
        for rec in self.records:
            self.mccortex_builder.with_dna_sequence(rec)
        output_graph = self.mccortex_builder.build(self.tmpdir)

        completed_process = (
            runner.Cortexpy(SPAWN_PROCESS).assemble(graph=output_graph,
                                                    initial_seqs=inital_seq_fasta)
        )
        assert completed_process.returncode == 0, completed_process
        return expectation.Fasta(completed_process.stdout)

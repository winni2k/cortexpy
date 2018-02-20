from pathlib import Path

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
    initial_seqs = attr.ib(attr.Factory(list))
    mccortex_builder = attr.ib(attr.Factory(builder.Mccortex))

    def with_records(self, *records):
        for rec in records:
            self.mccortex_builder.with_dna_sequence(rec)
        return self

    def with_dna_sequence(self, sequence, **kwargs):
        self.mccortex_builder.with_dna_sequence(sequence, **kwargs)
        return self

    def with_initial_sequences(self, *seqs):
        self.initial_seqs = seqs
        return self

    def with_kmer_size(self, size):
        self.mccortex_builder.with_kmer_size(size)
        return self

    def run(self):
        inital_seq_fasta = self.tmpdir / 'initial_seqs.fasta'
        initial_seqs = [SeqRecord(Seq(rec), id=str(idx)) for idx, rec in
                        enumerate(self.initial_seqs)]
        SeqIO.write(initial_seqs, str(inital_seq_fasta), "fasta")
        output_graph = self.mccortex_builder.build(self.tmpdir)

        completed_process = (
            runner.Cortexpy(SPAWN_PROCESS).assemble(graph=output_graph,
                                                    initial_seqs=inital_seq_fasta)
        )
        assert completed_process.returncode == 0, completed_process
        return expectation.Fasta(completed_process.stdout)


@attr.s(slots=True)
class Prune(object):
    tmpdir = attr.ib()
    mccortex_builder = attr.ib(attr.Factory(builder.Mccortex))
    min_tip_length = attr.ib(None)
    last_record = attr.ib(None)
    kmer_size = attr.ib(None)

    def with_records(self, *records):
        for rec in records:
            self.mccortex_builder.with_dna_sequence(rec)
            self.last_record = rec
        return self

    def prune_tips_less_than(self, n):
        self.min_tip_length = n
        return self

    def with_kmer_size(self, size):
        self.mccortex_builder.with_kmer_size(size)
        self.kmer_size = size
        return self

    def run(self):
        import networkx as nx
        mccortex_graph = self.mccortex_builder.build(self.tmpdir)

        cortexpy_graph = self.tmpdir / 'cortexpy_graph.pickle'
        initial_contig = self.last_record[0:self.kmer_size]
        ctp_runner = runner.Cortexpy(SPAWN_PROCESS)
        ctp_runner.traverse(graph=mccortex_graph, out=cortexpy_graph, contig=initial_contig)

        pruned_graph = Path(cortexpy_graph).with_suffix('.pruned.pickle')
        completed_process = ctp_runner.prune(graph=cortexpy_graph, out=pruned_graph,
                                             remove_tips=self.min_tip_length)

        assert completed_process.returncode == 0, completed_process

        return expectation.KmerGraphExpectation(nx.read_gpickle(str(pruned_graph)))

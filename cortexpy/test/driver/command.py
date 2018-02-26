from pathlib import Path

import Bio
import attr
import os

import networkx as nx
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
    records = attr.ib(attr.Factory(list))
    kmer_size = attr.ib(None)

    def with_records(self, *records):
        for rec in records:
            self.mccortex_builder.with_dna_sequence(rec)
            self.records.append(SeqRecord(Seq(rec)))
        return self

    def prune_tips_less_than(self, n):
        self.min_tip_length = n
        return self

    def with_kmer_size(self, size):
        self.mccortex_builder.with_kmer_size(size)
        self.kmer_size = size
        return self

    def run(self):
        mccortex_graph = self.mccortex_builder.build(self.tmpdir)

        cortexpy_graph = self.tmpdir / 'cortexpy_graph.pickle'
        contig_fasta = self.tmpdir / 'initial_contigs.fa'
        with open(str(contig_fasta), 'w') as fh:
            SeqIO.write(self.records, fh, 'fasta')
        ctp_runner = runner.Cortexpy(SPAWN_PROCESS)
        ctp_runner.traverse(graph=mccortex_graph, out=cortexpy_graph, contig=contig_fasta,
                            contig_fasta=True, subgraphs=True)

        pruned_graph = Path(cortexpy_graph).with_suffix('.pruned.pickle')
        completed_process = ctp_runner.prune(graph=cortexpy_graph, out=pruned_graph,
                                             remove_tips=self.min_tip_length)

        assert completed_process.returncode == 0, completed_process

        subgraphs = list(load_graph_stream(str(pruned_graph)))
        return expectation.graph.KmerGraphsExpectation(subgraphs)


@attr.s(slots=True)
class Traverse(object):
    """Runner for traverse acceptance tests"""
    tmpdir = attr.ib()
    mccortex_builder = attr.ib(attr.Factory(builder.Mccortex))
    traversal_contigs = attr.ib(None)
    added_records = attr.ib(attr.Factory(list))
    output_subgraphs = attr.ib(False)
    cortexpy_graph = attr.ib(init=False)

    def with_records(self, *records):
        for rec in records:
            self.mccortex_builder.with_dna_sequence(rec)
            self.added_records.append(rec)
        return self

    def with_initial_contigs(self, *contigs):
        self.traversal_contigs = contigs

    def with_subgraph_output(self):
        self.output_subgraphs = True
        return self

    def with_kmer_size(self, size):
        self.mccortex_builder.with_kmer_size(size)
        return self

    def run(self):
        mccortex_graph = self.mccortex_builder.build(self.tmpdir)
        contig_fasta = self.tmpdir / 'cortexpy_initial_contigs.fa'
        if self.traversal_contigs:
            initial_contigs = self.traversal_contigs
        else:
            initial_contigs = self.added_records
        with open(str(contig_fasta), 'w') as fh:
            Bio.SeqIO.write([SeqRecord(Seq(s)) for s in initial_contigs], fh, 'fasta')

        self.cortexpy_graph = self.tmpdir / 'cortexpy_graph.pickle'
        ctp_runner = runner.Cortexpy(SPAWN_PROCESS)
        completed_process = ctp_runner.traverse(graph=mccortex_graph,
                                                out=self.cortexpy_graph,
                                                contig=contig_fasta,
                                                contig_fasta=True,
                                                subgraphs=self.output_subgraphs)

        subgraphs = list(load_graph_stream(self.cortexpy_graph))
        assert completed_process.returncode == 0, completed_process

        return expectation.graph.KmerGraphsExpectation(subgraphs)


@attr.s(slots=True)
class ViewTraversal(object):
    """Runner for view of traversal acceptance tests"""
    tmpdir = attr.ib()
    traverse_driver = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.traverse_driver = Traverse(self.tmpdir)

    def with_records(self, *records):
        self.traverse_driver.with_records(*records)
        return self

    def with_subgraph_output(self):
        self.traverse_driver.with_subgraph_output()
        return self

    def with_kmer_size(self, size):
        self.traverse_driver.with_kmer_size(size)
        return self

    def run(self):
        self.traverse_driver.run()
        ret = runner.Cortexpy(SPAWN_PROCESS).view(
            cortexpy_graph=self.traverse_driver.cortexpy_graph)
        assert ret.returncode == 0, ret
        return expectation.Fasta(ret.stdout)


def load_graph_stream(path):
    with open(str(path), 'rb') as fh:
        while True:
            try:
                yield nx.read_gpickle(fh)
            except EOFError:
                break

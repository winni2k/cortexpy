import logging
from collections import OrderedDict
from pathlib import Path

import attr
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from attr import Factory

import cortexpy.test.runner as runner

logger = logging.getLogger(__name__)


@attr.s(slots=True)
class Mccortex:
    kmer_size = attr.ib(3)
    sequences = attr.ib(Factory(OrderedDict))
    mccortex_bin = attr.ib('mccortex')

    def with_kmer_size(self, kmer_size):
        self.kmer_size = kmer_size
        return self

    def with_dna_sequence(self, sequence, *, name='sample_0'):
        if name not in self.sequences:
            self.sequences[name] = []
        self.sequences[name].append(sequence)
        return self

    def build(self, tmpdir):
        mccortex_args = ['build', '--force', '--sort', '--kmer', str(self.kmer_size)]
        counter = 0
        for name, dna_sequences in self.sequences.items():
            input_fasta = str(tmpdir.join('input.{}.fasta'.format(name)))
            with open(input_fasta, 'w') as fh:
                for sequence in dna_sequences:
                    fh.write(
                        SeqRecord(Seq(sequence), id=str(counter), description='').format('fasta'))
                    counter += 1
            mccortex_args.extend(['--sample', name, '-1', input_fasta])

        output_graph = str(tmpdir.join('output.ctx'))
        mccortex_args.append(output_graph)

        ret = runner.Mccortex(self.kmer_size, mccortex_bin=self.mccortex_bin).run(mccortex_args)
        logger.debug('\n' + ret.stdout.decode())
        logger.debug('\n' + ret.stderr.decode())

        ret = runner.Mccortex(self.kmer_size, mccortex_bin=self.mccortex_bin).view(output_graph)
        logger.debug('\n' + ret.stdout.decode())
        logger.debug('\n' + ret.stderr.decode())

        return output_graph


@attr.s(slots=True)
class MccortexLinks:
    kmer_size = attr.ib(3)
    link_sequences = attr.ib(Factory(list))
    mccortex_bin = attr.ib('mccortex')

    def with_link_dna_sequence(self, sequence):
        self.link_sequences.append(sequence)
        return self

    def with_link_dna_sequences(self, *seqs):
        for seq in seqs:
            self.with_link_dna_sequence(seq)
        return self

    def build(self, tmpdir, mccortex_graph):
        mccortex_graph = Path(mccortex_graph)
        input_fasta = str(tmpdir / f'link_input.fasta')
        SeqIO.write([SeqRecord(Seq(s), id=str(i), description='') for i, s in
                     enumerate(self.link_sequences)], input_fasta, 'fasta')
        out_links = mccortex_graph.with_suffix('.ctp.gz')

        args = f'thread -1 {input_fasta} -o {out_links} {mccortex_graph}:0'.split()
        ret = runner.Mccortex(self.kmer_size, mccortex_bin=self.mccortex_bin).run(args)
        print(ret.stdout.decode())
        print(ret.stderr.decode(),file=sys.stderr)
        assert 0 == ret.returncode
        return out_links


@attr.s(slots=True)
class MccortexGraphLinks:
    graph_builder = attr.ib(attr.Factory(Mccortex))
    link_builder = attr.ib(attr.Factory(MccortexLinks))

    def __getattr__(self, name):
        if name in ['with_link_dna_sequence', 'with_link_dna_sequences']:
            return getattr(self.link_builder, name)
        if name in ['with_kmer_size', 'with_dna_sequence']:
            return getattr(self.graph_builder, name)
        raise NotImplementedError(f'Could not find attribute {name}')

    def build(self, tmpdir):
        graph = self.graph_builder.build(tmpdir)
        return graph, self.link_builder.build(tmpdir, graph)

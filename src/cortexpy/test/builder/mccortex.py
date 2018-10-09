import logging
from collections import OrderedDict
from pathlib import Path

import attr
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from attr import Factory

import cortexpy.test.runner as runner

logger = logging.getLogger(__name__)


@attr.s(slots=True)
class Mccortex(object):
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
                    fh.write(SeqRecord(Seq(sequence), id=str(counter), description='').format('fasta'))
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
class MccortexLinks(object):
    mccortex_graph = attr.ib()
    kmer_size = attr.ib(3)
    link_sequences = attr.ib(Factory(list))
    mccortex_bin = attr.ib('mccortex')

    def __attrs_post_init__(self):
        self.mccortex_graph = Path(self.mccortex_graph)

    def with_link_dna_sequence(self, sequence):
        self.link_sequences.append(sequence)
        return self

    def build(self, tmpdir):
        input_fasta = str(tmpdir / f'link_input.fasta')
        SeqIO.write([SeqRecord(Seq(s), id=str(i), description='') for i,s in enumerate(self.link_sequences)], input_fasta, 'fasta')
        out_links = self.mccortex_graph.with_suffix('.ctp.gz')

        args = f'thread -1 {input_fasta} -o {out_links} {self.mccortex_graph}'.split()
        ret = runner.Mccortex(self.kmer_size, mccortex_bin=self.mccortex_bin).run(args)
        logger.info('\n' + ret.stdout.decode())
        logger.warning('\n' + ret.stderr.decode())
        return out_links

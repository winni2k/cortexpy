import cortexpy.test.runner as runner
from collections import OrderedDict

import attr
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from attr import Factory
import logging

logger = logging.getLogger(__name__)


@attr.s
class Mccortex(object):
    kmer_size = attr.ib(3)
    sequences = attr.ib(Factory(OrderedDict))
    mccortex_bin = attr.ib(None)

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
        for name, dna_sequences in self.sequences.items():
            input_fasta = str(tmpdir.join('input.{}.fasta'.format(name)))
            with open(input_fasta, 'w') as fh:
                for sequence in dna_sequences:
                    fh.write(SeqRecord(Seq(sequence)).format('fasta'))
            mccortex_args.extend(['--sample', name, '-1', input_fasta])

        output_graph = str(tmpdir.join('output.ctx'))
        mccortex_args.append(output_graph)

        ret = runner.Mccortex(self.kmer_size, mccortex_bin=self.mccortex_bin).run(mccortex_args)
        logger.debug('\n' + ret.stderr.decode())

        ret = runner.Mccortex(self.kmer_size, mccortex_bin=self.mccortex_bin).view(output_graph)
        logger.debug('\n' + ret.stdout.decode())

        return output_graph

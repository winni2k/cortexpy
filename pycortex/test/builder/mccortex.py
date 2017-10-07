import pycortex.test.runner as runner
from collections import defaultdict

import attr
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from attr import Factory


@attr.s
class Mccortex(object):
    kmer_size = attr.ib(3)
    sequences = attr.ib(Factory(lambda: defaultdict(list)))

    def with_kmer_size(self, kmer_size):
        self.kmer_size = kmer_size
        return self

    def with_dna_sequence(self, sequence, *, name='sample_0'):
        self.sequences[name].append(sequence)
        return self

    def build(self, tmpdir):
        mccortex_args = ['build', '--sort', '--kmer', str(self.kmer_size)]
        for name, dna_sequences in self.sequences.items():
            input_fasta = str(tmpdir.join('input.{}.fasta'.format(name)))
            with open(input_fasta, 'w') as fh:
                for sequence in dna_sequences:
                    fh.write(SeqRecord(Seq(sequence)).format('fasta'))
            mccortex_args.extend(['--sample', name, '-1', input_fasta])

        output_graph = str(tmpdir.join('output.ctx'))
        mccortex_args.append(output_graph)

        runner.Mccortex(self.kmer_size).run(mccortex_args)

        return output_graph

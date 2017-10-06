import pycortex.test.runner as runner

import attr
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from attr import Factory


@attr.s
class Mccortex(object):
    kmer_size = attr.ib(3)
    sequences = attr.ib(Factory(list))

    def with_kmer_size(self, kmer_size):
        self.kmer_size = kmer_size
        return self

    def with_dna_sequence(self, name, sequence):
        self.sequences.append((name, sequence))
        return self

    def build(self, tmpdir):
        mccortex_args = ['build', '--sort', '--kmer', str(self.kmer_size)]
        input_fasta = str(tmpdir.join('input.fasta'))
        with open(input_fasta, 'w') as fh:
            for name, dna_sequence in self.sequences:
                fh.write(SeqRecord(Seq(dna_sequence)).format('fasta'))
                mccortex_args.extend(['--sample', name, '-1', input_fasta])

        output_graph = str(tmpdir.join('output.ctx'))
        mccortex_args.append(output_graph)

        runner.Mccortex(self.kmer_size).run(mccortex_args)

        return output_graph

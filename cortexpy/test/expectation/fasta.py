import io

import sys
from Bio import SeqIO
import attr


@attr.s(slots=True)
class Fasta(object):
    fasta_string = attr.ib()
    fasta_records = attr.ib(init=False)

    def __attrs_post_init__(self):
        print(self.fasta_string, file=sys.stderr)
        self.fasta_records = [rec for rec in SeqIO.parse(io.StringIO(self.fasta_string), "fasta")]
        print(self.fasta_records, file=sys.stderr)

    def has_n_records(self, n):
        assert len(self.fasta_records) == n
        return self

    def has_record(self, record):
        assert next((rec for rec in self.fasta_records if rec.seq == record), None).seq == record
        return self

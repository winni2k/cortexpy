import io
from collections import defaultdict

from Bio import SeqIO
import attr
import logging

logger = logging.getLogger(__name__)


@attr.s(slots=True)
class BioSeqRecord(object):
    record = attr.ib()

    def has_id_in(self, *ids):
        ids = {str(val) for val in ids}
        assert self.record.id in ids
        return self


@attr.s(slots=True)
class Fasta(object):
    fasta_string = attr.ib()
    fasta_records = attr.ib(init=False)
    fasta_record_dict = attr.ib(init=False)
    fasta_groups = attr.ib(attr.Factory(lambda: defaultdict(lambda: 0)))

    def __attrs_post_init__(self):
        self.fasta_records = [rec for rec in SeqIO.parse(io.StringIO(self.fasta_string), "fasta")]
        self.fasta_record_dict = {str(rec.seq): rec for rec in self.fasta_records}
        for rec in self.fasta_records:
            print(rec.id)
            group_id = rec.id.split('_')[0]
            print(group_id)
            self.fasta_groups[group_id] += 1
            print(self.fasta_groups)

    def has_n_records(self, n):
        assert n == len(self.fasta_records)
        return self

    def has_n_groups(self, n):
        assert n == len(self.fasta_groups.keys())
        return self

    def has_record_ids(self, *ids):
        assert set(ids) == {rec.id for rec in self.fasta_records}
        return self

    def has_record(self, expected_record_seq):
        assert expected_record_seq in self.fasta_record_dict.keys()
        return BioSeqRecord(self.fasta_record_dict[expected_record_seq])

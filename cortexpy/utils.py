from datetime import datetime

import attr
from Bio.Seq import reverse_complement
from functools import lru_cache


@lru_cache()
def revcomp(dna_string):
    return reverse_complement(dna_string)


@lru_cache()
def lexlo(kmer_string):
    """Returns the lexicographically lowest version of the kmer string"""
    alt_kmer_string = revcomp(kmer_string)
    if alt_kmer_string < kmer_string:
        return alt_kmer_string
    return kmer_string


@attr.s(slots=True)
class IntervalLogger(object):
    logger = attr.ib()
    min_log_interval_seconds = attr.ib(0)
    last_log_time = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.last_log_time = datetime.now()

    def _ok_to_log(self):
        return (datetime.now() - self.last_log_time).total_seconds() > self.min_log_interval_seconds

    def info(self, *args, **kwargs):
        if self._ok_to_log():
            self.last_log_time = datetime.now()
            return self.logger.info(*args, **kwargs)

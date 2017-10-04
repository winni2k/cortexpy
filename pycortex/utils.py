from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


def revcomp(kmer_string):
    return str(Seq(kmer_string, IUPAC.unambiguous_dna).reverse_complement())

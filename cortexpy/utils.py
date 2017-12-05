from Bio.Seq import reverse_complement


def revcomp(dna_string):
    return reverse_complement(dna_string)


def lexlo(kmer_string):
    """Returns the lexicographically lowest version of the kmer string"""
    alt_kmer_string = revcomp(kmer_string)
    if alt_kmer_string < kmer_string:
        return alt_kmer_string
    return kmer_string

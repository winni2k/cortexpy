import attr


@attr.s(slots=True)
class Header(object):
    kmer_size = attr.ib()
    kmer_container_size = attr.ib()
    num_colors = attr.ib()

import attr


@attr.s(slots=True)
class LinksExpectation:
    links = attr.ib()

    def has_n_link_groups(self, n):
        assert n == len(self.links.body.values())
        return self

    def has_link_group_for_kmer(self, kmer):
        assert kmer in self.links.body
        return LinkGroupExpectation(self.links.body[kmer])


@attr.s(slots=True)
class LinkGroupExpectation:
    group = attr.ib()

    def has_links(self, *expected_links):
        links = [
            f'{l.orientation.name} {l.num_juncs} {",".join([str(c) for c in l.counts])} {l.juncs}'
            for l in self.group.link_lines
        ]
        assert sorted(list(expected_links)) == sorted(links)
        return self

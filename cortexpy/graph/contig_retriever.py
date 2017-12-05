from collections import defaultdict

import attr
import networkx as nx

from cortexpy.graph import parser as parser
from cortexpy.graph.parser.kmer import EmptyKmerBuilder, flip_kmer_string_to_match


@attr.s(slots=True)
class ContigRetriever(object):
    graph_handle = attr.ib()
    graph_parser = attr.ib(init=False)
    num_colors = attr.ib(init=False)
    contig_color = attr.ib(init=False)
    non_contig_colors = attr.ib(init=False)
    colors = attr.ib(init=False)
    empty_kmer_builder = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.graph_parser = parser.RandomAccess(self.graph_handle)
        self.num_colors = self.graph_parser.header.num_colors + 1
        self.contig_color = self.num_colors - 1
        self.non_contig_colors = list(range(self.num_colors - 1))
        self.colors = self.non_contig_colors + [self.contig_color]
        self.empty_kmer_builder = EmptyKmerBuilder(num_colors=self.num_colors)

    def get_kmers(self, contig):
        kmer_size = self.graph_parser.header.kmer_size
        assert len(contig) >= kmer_size
        kmers = []
        for kmer_start in range(len(contig) - kmer_size + 1):
            kmer_string = contig[kmer_start:(kmer_start + kmer_size)]
            try:
                kmer = self.graph_parser.get_kmer_for_string(kmer_string)
                kmer.append_color(1)
            except KeyError:
                kmer = self.empty_kmer_builder.build_or_get(kmer_string)
                if kmer.coverage[-1] == 0:
                    kmer.coverage[-1] = 1
            kmers.append((kmer, kmer_string))
        for kmer_idx in range(len(kmers) - 1):
            this_kmer = kmers[kmer_idx][0]
            next_kmer = kmers[kmer_idx + 1][0]

            parser.kmer.connect_kmers(this_kmer, next_kmer, self.contig_color)
        return kmers

    def get_kmer_graph(self, contig):
        kmer_graph = nx.MultiDiGraph()
        kmers = list(self.get_kmers(contig))
        for kmer, kmer_string in kmers:
            kmer_graph.add_node(kmer_string, kmer=kmer)
        for kmer_idx, (kmer, kmer_string) in enumerate(kmers):
            incoming_kmers, outgoing_kmers = get_incoming_and_outgoing_kmers(kmer)
            is_revcomp = kmer.kmer != kmer_string
            if is_revcomp:
                incoming_kmers, outgoing_kmers = outgoing_kmers, incoming_kmers
            if 0 < kmer_idx:
                prev_kmer, prev_kmer_string = kmers[kmer_idx - 1]
                if prev_kmer.kmer in incoming_kmers[self.contig_color]:
                    kmer_graph.add_edge(prev_kmer_string, kmer_string, key=self.contig_color)
            for is_incoming, neighbor_kmers in [(True, incoming_kmers), (False, outgoing_kmers)]:
                for color in self.non_contig_colors:
                    for neighbor_kmer_string in neighbor_kmers[color]:
                        neighbor_kmer_string, _ = flip_kmer_string_to_match(
                            neighbor_kmer_string,
                            kmer_string,
                            flip_is_after_reference_kmer=not is_incoming
                        )
                        if neighbor_kmer_string not in kmer_graph:
                            kmer_graph.add_node(
                                neighbor_kmer_string,
                                kmer=self.empty_kmer_builder.build_or_get(neighbor_kmer_string)
                            )
                        if is_incoming:
                            kmer_graph.add_edge(neighbor_kmer_string, kmer_string, key=color)
                        else:
                            kmer_graph.add_edge(kmer_string, neighbor_kmer_string, key=color)
        for kmer_node in kmer_graph:
            kmer_graph.nodes[kmer_node]['repr'] = kmer_node
        return kmer_graph


def get_incoming_and_outgoing_kmers(kmer):
    incoming_kmers = defaultdict(set)
    outgoing_kmers = defaultdict(set)
    for color, edge_set in enumerate(kmer.edges):
        for incoming_kmer in edge_set.get_incoming_kmers(kmer.kmer):
            incoming_kmers[color].add(incoming_kmer)
        for outgoing_kmer in edge_set.get_outgoing_kmers(kmer.kmer):
            outgoing_kmers[color].add(outgoing_kmer)
    return incoming_kmers, outgoing_kmers
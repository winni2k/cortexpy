from collections import defaultdict

import attr
import networkx as nx

from cortexpy.utils import lexlo
from . import parser as parser
from .parser.kmer import EmptyKmerBuilder, revcomp_target_to_match_ref

RETRIEVED_CONTIG_NAME = 'retrieved_contig'


@attr.s(slots=True)
class ContigRetriever(object):
    graph_handle = attr.ib()
    seen_kmer_strings = attr.ib(attr.Factory(dict))
    graph_parser = attr.ib(init=False)
    num_colors = attr.ib(init=False)
    contig_color = attr.ib(init=False)
    non_contig_colors = attr.ib(init=False)
    colors = attr.ib(init=False)
    empty_kmer_builder = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.graph_parser = parser.RandomAccess(self.graph_handle)
        self.num_colors = self.graph_parser.num_colors + 1
        self.contig_color = self.num_colors - 1
        self.non_contig_colors = list(range(self.num_colors - 1))
        self.colors = self.non_contig_colors + [self.contig_color]
        self.empty_kmer_builder = EmptyKmerBuilder(num_colors=self.num_colors, default_coverage=0)

    def get_kmers(self, contig):
        kmer_size = self.graph_parser.kmer_size
        assert len(contig) >= kmer_size
        kmers = []
        for kmer_start in range(len(contig) - kmer_size + 1):
            kmer_string = contig[kmer_start:(kmer_start + kmer_size)]
            lexlo_kmer_string = lexlo(kmer_string)
            if lexlo_kmer_string in self.seen_kmer_strings:
                kmer = self.seen_kmer_strings[lexlo_kmer_string]
            else:
                try:
                    kmer = self.graph_parser.get_kmer_for_string(kmer_string)
                except KeyError:
                    kmer = self.empty_kmer_builder.build(kmer_string)
                self.seen_kmer_strings[lexlo_kmer_string] = kmer
            kmer.increment_color_coverage(self.num_colors - 1)
            kmers.append((kmer, kmer_string))
        for kmer_idx in range(len(kmers) - 1):
            this_kmer = kmers[kmer_idx][0]
            next_kmer = kmers[kmer_idx + 1][0]

            parser.kmer.connect_kmers(this_kmer, next_kmer, self.contig_color,
                                      identical_kmer_check=False)
        return kmers

    def get_kmer_graph(self, contig):
        kmer_graph = nx.MultiDiGraph()
        kmers = list(self.get_kmers(contig))
        for kmer, kmer_string in kmers:
            kmer_graph.add_node(kmer_string, kmer=kmer)
        for kmer_idx, (kmer, kmer_string) in enumerate(kmers):
            incoming_kmers, outgoing_kmers = get_incoming_and_outgoing_kmers(kmer)
            not_lexlo = kmer.kmer != kmer_string
            if not_lexlo:
                incoming_kmers, outgoing_kmers = outgoing_kmers, incoming_kmers
            if 0 < kmer_idx:
                prev_kmer, prev_kmer_string = kmers[kmer_idx - 1]
                if prev_kmer.kmer in incoming_kmers[self.contig_color]:
                    kmer_graph.add_edge(prev_kmer_string, kmer_string, key=self.contig_color)
            for is_incoming, neighbor_kmers in [(True, incoming_kmers), (False, outgoing_kmers)]:
                for color in self.non_contig_colors:
                    for neighbor_kmer_string in neighbor_kmers[color]:
                        neighbor_kmer_string, _ = revcomp_target_to_match_ref(
                            neighbor_kmer_string,
                            kmer_string,
                            rc_is_after_reference_kmer=not is_incoming
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

        for kmer_node in list(kmer_graph):
            kmer_graph.nodes[kmer_node]['repr'] = kmer_node

        return self._add_metadata_to_graph(kmer_graph)

    def _add_metadata_to_graph(self, kmer_graph):
        kmer_graph.graph['colors'] = tuple(self.colors)
        sample_names = [n.decode() for n in self.graph_parser.sample_names]
        sample_names.append(RETRIEVED_CONTIG_NAME)
        kmer_graph.graph['sample_names'] = sample_names
        kmer_graph.graph['is_kmer_graph'] = True
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

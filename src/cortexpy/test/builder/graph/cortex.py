import io
import json
from collections import defaultdict

import attr
from delegation import SingleDelegated

import cortexpy.graph.cortex as cortex
from cortexpy.graph.interactor import Interactor
from cortexpy.graph.parser.kmer import EmptyKmerBuilder
from cortexpy.graph.parser.random_access import load_ra_cortex_graph
from cortexpy.links import Links
from cortexpy.test.builder import Graph
from cortexpy.utils import lexlo


def get_cortex_builder():
    return CortexBuilder(Graph())


def get_cortex_graph_mapping_builder():
    return CortexGraphMappingBuilder(Graph())


class CortexBuilder(SingleDelegated):

    def build(self):
        return load_ra_cortex_graph(self.delegate.build())


@attr.s()
class CortexGraphMappingBuilder(SingleDelegated):
    delegate = attr.ib()
    ra_parser_args = attr.ib(attr.Factory(dict))

    def build(self):
        return load_ra_cortex_graph(self.delegate.build(),
                                    ra_parser_args=self.ra_parser_args)._kmer_mapping

    def with_kmer_cache_size(self, n):
        self.ra_parser_args['kmer_cache_size'] = n
        return self


@attr.s(slots=True)
class CortexGraphBuilder(object):
    graph = attr.ib(attr.Factory(cortex.CortexDiGraph))
    colors = attr.ib(attr.Factory(set))
    kmer_builder = attr.ib(attr.Factory(EmptyKmerBuilder))
    consistent_seeds = attr.ib(None)

    def __attrs_post_init__(self):
        self.with_colors(0)

    def with_kmer_size(self, k):
        if 'kmer_size' in self.graph.graph:
            assert self.graph.graph['kmer_size'] == k
        else:
            self.graph.graph['kmer_size'] = k
        return self

    def with_node_coverage(self, node, *coverage):
        if node not in self.graph:
            self.add_node(node)
        self.graph.node[node].coverage = tuple(coverage)
        return self

    def with_node_kmer(self, node, kmer):
        self.with_kmer_size(len(kmer.kmer))
        self.graph.add_node(node, kmer=kmer)
        return self

    def with_node(self, node):
        return self.add_node(node)

    def add_node(self, node):
        if node not in self.graph:
            kmer = self.kmer_builder.build_or_get(node)
            self.with_kmer_size(len(kmer.kmer))
            self.with_node_kmer(node, kmer)
        return self

    def with_colors(self, *colors):
        for color in colors:
            self.with_color(color)
        return self

    def with_color(self, color):
        self.colors.add(color)
        self.kmer_builder.num_colors = len(self.colors)
        self.graph.graph['colors'] = sorted(list(self.colors))
        return self

    def make_consistent(self, *seeds):
        self.consistent_seeds = seeds
        return self

    def add_edge(self, u, v, color=0, coverage=1):
        self.add_edge_with_color_and_coverage(u, v, color, color_coverage=coverage)
        return self

    def add_edge_with_color_and_coverage(self, u, v, color, color_coverage):
        for node in [u, v]:
            if node in self.graph:
                coverages = list(self.graph[node].kmer.coverage)
            else:
                coverages = [0 for _ in range(len(self.colors))]
            coverages[color] = color_coverage
            self.with_node_coverage(node, *coverages)
        self.graph.add_edge(u, v, key=color)

    def add_edge_with_coverage(self, u, v, color_coverage):
        self.add_edge(u, v, 0, color_coverage)
        return self

    def add_edge_with_color(self, u, v, color):
        self.add_edge(u, v, color, 1)
        return self

    def add_path(self, *k_strings, color=0, coverage=1):
        if len(k_strings) == 1 and isinstance(k_strings[0], list):
            k_strings = k_strings[0]
        kmer = self.kmer_builder.build_or_get(k_strings[0])

        num_colors = len(self.colors)
        self.graph.add_node(k_strings[0], kmer=kmer)
        self.with_node_coverage(kmer.kmer, *[coverage for _ in range(num_colors)])
        if len(k_strings) > 1:
            for k_string1, k_string2 in zip(k_strings[:-1], k_strings[1:]):
                kmer = self.kmer_builder.build_or_get(k_string2)
                self.graph.add_node(k_string2, kmer=kmer)
                self.with_node_coverage(kmer.kmer, *[coverage for _ in range(num_colors)])
                self.add_edge_with_color(k_string1, k_string2, color)

    def build(self):
        if self.consistent_seeds:
            self.graph = Interactor.from_graph(self.graph) \
                .make_graph_nodes_consistent(self.consistent_seeds) \
                .graph
        return self.graph


@attr.s(slots=True)
class LinksBuilder:
    header = attr.ib(None)
    links = attr.ib(None)

    def __attrs_post_init__(self):
        if self.header is None:
            self.header = {'graph': {'num_colours': 1}}
        if self.links is None:
            self.links = defaultdict(list)

    def with_link_for_kmer(self, link, kmer):
        "link: <F|R> <num_juncs> <counts0,counts1,...> <junctions>"
        assert lexlo(kmer) == kmer
        self.links[kmer].append(link)
        return self

    def build(self):
        header = json.dumps(self.header)
        body = []
        for kmer in sorted(self.links.keys()):
            group = self.links[kmer]
            body.append(f'{kmer} {len(group)}')
            body += group
        out = '\n'.join([header] + body)
        return Links.from_binary_stream(io.BufferedReader(io.BytesIO(out.encode())))

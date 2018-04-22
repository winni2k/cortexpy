import attr
import networkx as nx
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging

from networkx.algorithms.coloring import strategy_connected_sequential_dfs

from cortexpy.constants import EdgeTraversalOrientation, NodeEdgeDirection
from .colored_de_bruijn import ColoredDeBruijn

logger = logging.getLogger(__name__)


def make_multi_graph(graph):
    if isinstance(graph, ColoredDeBruijn):
        return graph
    return nx.MultiGraph(graph)


@attr.s(slots=True)
class Interactor(object):
    graph = attr.ib()
    colors = attr.ib()

    def add_edge_to_graph_for_kmer_pair(self, kmer1, kmer2, kmer1_string, kmer2_string):
        if isinstance(self.graph, ColoredDeBruijn):
            return self
        first_is_revcomp = bool(kmer1.kmer != kmer1_string)
        for color in self.colors:
            if first_is_revcomp:
                add_edge = kmer1.has_incoming_edge_from_kmer_in_color(kmer2, color)
            else:
                add_edge = kmer1.has_outgoing_edge_to_kmer_in_color(kmer2, color)
            if add_edge:
                self.graph.add_edge(kmer1_string, kmer2_string, key=color)
        return self

    def prune_tips_less_than(self, n):
        nodes_to_prune = set()
        assert self.colors is None
        graph = make_multi_graph(self.graph)
        num_tips = 0
        num_tips_to_prune = 0
        for node, direction in edge_nodes_of(self.graph):
            tip = find_tip_from(
                n=n,
                graph=self.graph,
                next_node_generator=nx.dfs_preorder_nodes(graph, source=node)
            )
            num_tips += 1
            if tip:
                nodes_to_prune |= tip
                num_tips_to_prune += 1
        logger.info('Found %s of %s tips shorter than %s.', num_tips_to_prune, num_tips, n)
        logger.info('Pruning %s nodes', len(nodes_to_prune))
        self.graph.remove_nodes_from(nodes_to_prune)
        return self


def find_tip_from(*, n, graph, next_node_generator):
    tip = set()
    query = next(next_node_generator)
    for _ in range(n):
        if len(set(graph.in_edges(query))) > 1 or len(set(graph.out_edges(query))) > 1:
            return tip
        else:
            tip.add(query)
            try:
                query = next(next_node_generator)
            except StopIteration:
                return tip
    return None


def make_copy_of_color(graph, color, include_self_refs=False):
    """Makes a copy of graph, but only copies over links with key=color.
    Only copies over nodes that are linked by a link with key=color.
    """
    out_graph = graph.fresh_copy()
    for u, v, key, data in graph.edges(keys=True, data=True):
        if key == color:
            if u == v and not include_self_refs:
                continue
            out_graph.add_edge(u, v, key=key, **data)
    for node in list(out_graph):
        out_graph.add_node(node, **graph.node[node])
    return out_graph


def convert_kmer_path_to_contig(path):
    if len(path) == 0:
        return ''
    contig = [path[0]]
    for kmer in path[1:]:
        contig.append(kmer[-1])
    return ''.join(contig)


def in_nodes_of(graph):
    for source in graph.nodes():
        if graph.in_degree(source) == 0:
            yield source


def out_nodes_of(graph):
    for target in graph.nodes():
        if graph.out_degree(target) == 0:
            yield target


def edge_nodes_of(graph):
    for node in graph.nodes():
        if graph.out_degree(node) == 0:
            yield (node, NodeEdgeDirection.outgoing)
        if graph.in_degree(node) == 0:
            yield (node, NodeEdgeDirection.incoming)


@attr.s(slots=True)
class Contigs(object):
    graph = attr.ib()
    color = attr.ib(None)

    def all_simple_paths(self):
        if self.color is not None:
            graph = make_copy_of_color(self.graph, self.color, include_self_refs=False)
        else:
            graph = self.graph
        idx = 0
        for source in in_nodes_of(graph):
            for target in out_nodes_of(graph):
                if source == target:
                    continue
                for path in nx.all_simple_paths(graph, source=source, target=target):
                    yield SeqRecord(Seq(convert_kmer_path_to_contig(path)), id=str(idx),
                                    description='')
                    idx += 1

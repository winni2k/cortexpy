import attr
import networkx as nx
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging

from cortexpy.constants import EdgeTraversalOrientation
from cortexpy.graph.parser.kmer import revcomp_target_to_match_ref
from cortexpy.utils import lexlo, revcomp
from .cortex import CortexDiGraph, ConsistentCortexDiGraph

logger = logging.getLogger(__name__)


def make_multi_graph(graph):
    if isinstance(graph, nx.Graph):
        return nx.MultiGraph(graph)
    if isinstance(graph, nx.DiGraph):
        return nx.MultiDiGraph(graph)
    return graph


@attr.s(slots=True)
class Interactor(object):
    graph = attr.ib()
    colors = attr.ib()

    def add_edge_to_graph_for_kmer_pair(self, kmer1, kmer2, kmer1_string, kmer2_string):
        if isinstance(self.graph, CortexDiGraph):
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
        for node, direction_of_end in edge_nodes_of(self.graph):
            tip = find_tip_from(
                n=n,
                graph=self.graph,
                start=node,
                next_node_generator=node_generator_from_edges(
                    nx.edge_dfs(
                        graph,
                        source=node,
                        orientation=EdgeTraversalOrientation.other(direction_of_end).name
                    )
                )
            )
            num_tips += 1
            if tip:
                nodes_to_prune |= tip
                num_tips_to_prune += 1
        logger.info('Found %s of %s tips shorter than %s.', num_tips_to_prune, num_tips, n)
        logger.info('Pruning %s nodes', len(nodes_to_prune))
        self.graph.remove_nodes_from(nodes_to_prune)
        return self

    def make_graph_nodes_consistent(self, seed_kmer_strings=None):
        """
        Take a CDBG and make all nodes have kmer_strings that are consistent with each other.
        If a seed kmer string is provided, then start with that seed kmer.
        """
        graph = CortexDiGraph(self.graph)
        new_graph = ConsistentCortexDiGraph(graph=self.graph.graph)

        seeds = SeedKmerStringIterator(self.graph.nodes(), seed_kmer_strings)

        for seed, lexlo_seed in seeds:
            new_graph.add_node(seed, kmer=self.graph.node[lexlo_seed])
            seeds.remove(lexlo_seed)
            for source, sink, key, direction in nx.edge_dfs(graph, lexlo_seed, 'ignore'):
                if direction == 'forward':
                    rc_after_ref_kmer = True
                    ref, target = source, sink
                elif direction == 'reverse':
                    ref, target = sink, source
                    rc_after_ref_kmer = False
                else:
                    raise Exception("unknown direction: {}".format(direction))
                if ref not in new_graph.node:
                    ref = revcomp(ref)
                    rc_after_ref_kmer = not rc_after_ref_kmer
                matched_target, _ = revcomp_target_to_match_ref(target, ref, rc_after_ref_kmer)
                new_graph.add_node(matched_target, kmer=graph.node[matched_target])
                seeds.remove(lexlo(matched_target))
        self.graph = new_graph
        return self


@attr.s(slots=True)
class SeedKmerStringIterator(object):
    all_lexlo_kmer_strings = attr.ib()
    seed_kmer_strings = attr.ib(attr.Factory(list))
    _unseen_lexlo_kmer_strings = attr.ib(attr.Factory(list))
    _seen_lexlo_kmer_strings = attr.ib(attr.Factory(set))

    def __attrs_post_init__(self):
        self._unseen_lexlo_kmer_strings = {lexlo(k_string) for k_string in
                                           self.all_lexlo_kmer_strings}

    def __iter__(self):
        return self

    def __next__(self):
        while self.seed_kmer_strings:
            seed = self.seed_kmer_strings.pop()
            lexlo_seed = lexlo(seed)
            if lexlo_seed not in self._unseen_lexlo_kmer_strings:
                continue
            self._seen_lexlo_kmer_strings.add(lexlo_seed)
            return seed, lexlo_seed
        while self._unseen_lexlo_kmer_strings:
            unseen = self._unseen_lexlo_kmer_strings.pop()
            if unseen in self._seen_lexlo_kmer_strings:
                continue
            self._seen_lexlo_kmer_strings.add(unseen)
            return unseen, unseen
        raise StopIteration

    def remove(self, lexlo_kmer_string):
        assert lexlo_kmer_string == lexlo(lexlo_kmer_string)
        self._seen_lexlo_kmer_strings.add(lexlo_kmer_string)


def seed_kmer_string_generator(seed_kmer_strings, unseen_lexlo_kmer_strings):
    while seed_kmer_strings:
        seed = seed_kmer_strings.pop()
        lexlo_seed = lexlo(seed)
        if lexlo_seed not in unseen_lexlo_kmer_strings:
            continue
        yield seed, lexlo_seed


def node_generator_from_edges(edge_generator):
    for edge in edge_generator:
        if len(edge) == 2:
            yield edge[1]
        elif len(edge) == 3:
            if edge[2] == EdgeTraversalOrientation.reverse.name:
                yield edge[0]
            else:
                yield edge[1]
        elif len(edge) == 4:
            if edge[3] == EdgeTraversalOrientation.reverse.name:
                yield edge[0]
            elif edge[3] == EdgeTraversalOrientation.original.name:
                yield edge[1]
            else:
                ValueError
        else:
            raise ValueError


def find_tip_from(*, n, start, graph, next_node_generator):
    tip = set()
    query = start
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
    """Find all edge nodes of a graph
    second return value is direction of edge
    """
    for node in graph.nodes():
        if graph.out_degree(node) == 0:
            yield (node, EdgeTraversalOrientation.original)
        if graph.in_degree(node) == 0:
            yield (node, EdgeTraversalOrientation.reverse)


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

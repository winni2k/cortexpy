import collections
import logging
from collections import OrderedDict

import attr
import networkx as nx
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from cortexpy.constants import EdgeTraversalOrientation
from cortexpy.graph.cortex import CortexDiGraph, ConsistentCortexDiGraph
from cortexpy.graph.parser.kmer import revcomp_target_to_match_ref
from cortexpy.graph.serializer.unitig import UnitigCollapser
from cortexpy.links import UnitigLinkWalker, LinkedGraphTraverser
from cortexpy.utils import lexlo, revcomp

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

    @classmethod
    def from_graph(cls, graph):
        return cls(graph=graph)

    def compose_in_graph(self, graph_to_add):
        """Compose two graphs and store the result in the interactor graph"""
        for g in [graph_to_add, self.graph]:
            assert g.is_multigraph()
            assert g.is_directed()
        self.graph.add_nodes_from(graph_to_add.nodes(data=True))
        if isinstance(graph_to_add, nx.Graph):
            self.graph.add_edges_from(graph_to_add.edges(keys=True))

    def find_nodes_of_tips_less_than(self, n):
        nodes_to_prune = set()
        graph = make_multi_graph(self.graph)
        num_tips = 0
        num_tips_to_prune = 0
        for node, direction_of_end in list(edge_nodes_of(self.graph)):
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
        return nodes_to_prune

    def prune_tips_less_than(self, n):
        logger.info(f'Removing tips shorter than {n} k-mers')
        nodes_to_prune = self.find_nodes_of_tips_less_than(n)
        logger.info('Pruning %s nodes', len(nodes_to_prune))
        self.graph.remove_nodes_from(nodes_to_prune)
        return self

    def make_graph_nodes_consistent(self, seed_kmer_strings=None):
        """
        Take a Cortex graph and make all nodes have kmer_strings that are consistent with each
        other. If a seed kmer string is provided, then start with that seed kmer.
        """
        if self.graph.is_consistent():
            return self
        if seed_kmer_strings is None:
            seed_kmer_strings = []
        graph = CortexDiGraph(self.graph)
        new_graph = ConsistentCortexDiGraph(graph=self.graph.graph)

        seeds = SeedKmerStringIterator.from_all_kmer_strings_and_seeds(self.graph.nodes(),
                                                                       seed_kmer_strings)

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

    def keep_color(self, color):
        self.graph = make_copy_of_color_for_kmer_graph(self.graph, color, include_self_refs=False)
        return self

    def all_simple_paths(self, extra_incoming_node=None, links=None):
        if not isinstance(self.graph, nx.Graph):
            assert self.graph.is_consistent()
        if extra_incoming_node:
            for neighbor in self.graph.pred[extra_incoming_node]:
                for color in self.graph.graph['colors']:
                    self.graph.remove_edge(neighbor, extra_incoming_node, color)
        unitig_graph = UnitigCollapser(self.graph) \
            .collapse_kmer_unitigs() \
            .unitig_graph
        unitig_graph = nx.DiGraph(unitig_graph)
        unitig_graph = nx.convert_node_labels_to_integers(unitig_graph)
        path_converter = UnitigGraphPathConverter.from_unitig_graph(unitig_graph)

        record_idx = 0
        in_nodes = sorted(list(in_nodes_of(unitig_graph)))
        logger.info(f"Found {len(in_nodes)} incoming tip nodes")
        out_nodes = set(sorted(list(out_nodes_of(unitig_graph))))
        logger.info(f"Found {len(out_nodes)} outgoing tip nodes")
        for sidx, source in enumerate(in_nodes):
            if source in out_nodes:
                logger.info('Incoming node %s; %s outgoing nodes; Path number %s',
                            sidx,
                            len(out_nodes),
                            record_idx)
                yield SeqRecord(Seq(unitig_graph.node[source]['unitig']), id=str(record_idx),
                                description='')
                record_idx += 1
                continue
            if links is not None:
                wrapped_graph = LinkedGraphTraverser.from_graph_and_link_walker(
                    unitig_graph,
                    UnitigLinkWalker.from_links_unitigs_kmer_size_unitig(
                        links,
                        unitig_graph,
                        unitig_graph.graph['kmer_size'],
                        source
                    )
                )
                paths = _all_simple_paths_graph(wrapped_graph,
                                                source,
                                                out_nodes,
                                                cutoff=len(self.graph) - 1)
            else:
                paths = _all_simple_paths_graph(unitig_graph,
                                                source,
                                                out_nodes,
                                                cutoff=len(self.graph) - 1)

            for pidx, path in enumerate(paths):
                if pidx % 100000 == 0:
                    logger.info('Incoming node %s; %s outgoing nodes; Path number %s', sidx,
                                len(out_nodes), record_idx)
                yield SeqRecord(
                    Seq(path_converter.to_contig(path)),
                    id=str(record_idx),
                    description='')
                record_idx += 1


@attr.s(slots=True)
class SeedKmerStringIterator(object):
    """
    Iterates seeds and their lexlo representations that exist in the supplied all_kmers:
    >>> list(SeedKmerStringIterator.from_all_kmer_strings_and_seeds(['AAC'], ['GTT']))
    [('GTT', 'AAC')]

    Kmers that are not in the seed list are return after that:
    >>> list(SeedKmerStringIterator.from_all_kmer_strings_and_seeds(['AAA', 'AAC'], ['GTT']))
    [('GTT', 'AAC'), ('AAA', 'AAA')]

    Seeds that do not exist in the all_kmers are not returned.
    >>> list(SeedKmerStringIterator.from_all_kmer_strings_and_seeds([], ['CCC']))
    []

    Returned kmers from all_kmers list are returned in order.
    >>> list(SeedKmerStringIterator.from_all_kmer_strings_and_seeds(['AAA', 'AAG', 'AAC'], []))
    [('AAA', 'AAA'), ('AAG', 'AAG'), ('AAC', 'AAC')]
    """
    seed_kmer_strings = attr.ib()
    _unseen_lexlo_kmer_strings = attr.ib()
    _seen_lexlo_kmer_strings = attr.ib(attr.Factory(set))

    @classmethod
    def from_all_kmer_strings_and_seeds(cls, all_kmers, seeds):
        return cls(list(seeds),
                   OrderedDict.fromkeys(lexlo(k_string) for k_string in all_kmers))

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
            unseen, _ = self._unseen_lexlo_kmer_strings.popitem(last=False)
            if unseen in self._seen_lexlo_kmer_strings:
                continue
            self._seen_lexlo_kmer_strings.add(unseen)
            return unseen, unseen
        raise StopIteration

    def remove(self, lexlo_kmer_string):
        assert lexlo_kmer_string == lexlo(lexlo_kmer_string)
        self._seen_lexlo_kmer_strings.add(lexlo_kmer_string)


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


def make_copy_of_color_for_kmer_graph(graph, color, include_self_refs=False):
    """Makes a copy of graph, but only copies over links with key=color.
    Only copies over nodes that are linked by a link with key=color.
    """
    out_graph = nx.MultiDiGraph(colors=[color])
    for u, v, key, data in graph.edges(keys=True, data=True):
        if key == color:
            if u == v and not include_self_refs:
                continue
            out_graph.add_edge(u, v, key=key, **data)
    for node in list(out_graph):
        out_graph.add_node(node, kmer=graph.node[node])
    return out_graph


@attr.s(slots=True)
class UnitigGraphPathConverter:
    reprs = attr.ib()

    @classmethod
    def from_unitig_graph(cls, graph):
        reprs = {u: rpr for u, rpr in graph.nodes.data('repr')}
        return cls(reprs)

    def to_contig(self, path):
        return ''.join([self.reprs[u] for u in path])


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


def _all_simple_paths_with_links(G, source, target, links, cutoff=None):
    """This function was copied from Networkx before being edited by Warren Kretzschmar"""
    walker = UnitigLinkWalker.from_links_unitigs_kmer_size(links, G, G.graph['kmer_size'])
    walker.load_unitig(source)
    assert cutoff is not None
    if target in G:
        targets = {target}
    else:
        targets = set(target)
    visited = collections.OrderedDict.fromkeys([source])
    stack = [walker.link_successors()]
    while stack:
        children = stack[-1]
        child = next(children, None)
        if child is None:
            stack.pop()
            visited.popitem()
        elif len(visited) < cutoff:
            if child in targets:
                yield list(visited) + [child]
            elif child not in visited:
                visited[child] = None
                stack.append(iter(G[child]))
        else:  # len(visited) == cutoff:
            if child in targets or len(targets & children) != 0:
                yield list(visited) + [child]
            stack.pop()
            visited.popitem()


def _all_simple_paths_graph(G, source, target, cutoff=None):
    """This function was copied from Networkx before being edited by Warren Kretzschmar
    todo: switch back to nx.all_simple_paths once Networkx 2.2 is released"""
    assert cutoff is not None
    if target in G:
        targets = {target}
    else:
        targets = set(target)
    visited = collections.OrderedDict.fromkeys([source])
    stack = [iter(G[source])]
    while stack:
        children = stack[-1]
        child = next(children, None)
        if child is None:
            stack.pop()
            visited.popitem()
        elif len(visited) < cutoff:
            if child in targets:
                yield list(visited) + [child]
            elif child not in visited:
                visited[child] = None
                stack.append(iter(G[child]))
        else:  # len(visited) == cutoff:
            if child in targets or len(targets & children) != 0:
                yield list(visited) + [child]
            stack.pop()
            visited.popitem()

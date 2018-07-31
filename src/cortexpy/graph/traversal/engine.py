import copy
import attr
import collections
import networkx as nx

from ..interactor import Interactor
from cortexpy.graph.cortex import build_empty_cortex_graph_from_ra_parser
from cortexpy.graph.parser.kmer import EmptyKmerBuilder
from cortexpy.utils import lexlo, IntervalLogger, kmerize_contig, kmerize_fasta
from . import branch
from cortexpy.constants import EdgeTraversalOrientation, EngineTraversalOrientation
import logging

logger = logging.getLogger(__name__)


def add_graph_to(graph, graph_to_add):
    """Compose two graphs and store the result in the first graph"""
    for g in [graph_to_add, graph]:
        assert g.is_multigraph()
        assert g.is_directed()
    graph.add_nodes_from(graph_to_add.nodes(data=True))
    if isinstance(graph_to_add, nx.Graph):
        graph.add_edges_from(graph_to_add.edges(keys=True))


@attr.s(slots=True)
class Engine(object):
    """This engine creates subgraphs of Cortex graphs"""
    ra_parser = attr.ib()
    traversal_colors = attr.ib((0,))
    orientation = attr.ib(EngineTraversalOrientation.original)
    max_nodes = attr.ib(None)
    branch_queue = attr.ib(attr.Factory(collections.deque))
    graph = attr.ib(init=False)
    last_graph_size = attr.ib(0)
    logging_interval = attr.ib(0)
    queuer = attr.ib(init=False)
    branch_traverser = attr.ib(init=False)
    logger = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.graph = build_empty_cortex_graph_from_ra_parser(self.ra_parser)
        self._add_graph_metadata()
        self.logger = IntervalLogger(logger, min_log_interval_seconds=self.logging_interval)

    def traverse_from_each_kmer_in_fasta(self, fasta):
        kmer_generator = kmerize_fasta(fasta, self.ra_parser.kmer_size)
        self._traverse_from_each_kmer_in(kmer_generator)
        self._post_process_graph()
        return self

    def traverse_from_each_kmer_in(self, contig):
        self._traverse_from_each_kmer_in(kmerize_contig(contig, self.ra_parser.kmer_size))
        self._post_process_graph()
        return self

    def _traverse_from_each_kmer_in(self, kmer_generator):
        for start_kmer in kmer_generator:
            try:
                add_graph_to(self.graph, self._traverse_from(start_kmer).graph)
                self.log_graph_size()
            except KeyError:
                pass
            if self.max_nodes and len(self.graph) > self.max_nodes:
                raise Exception(("Terminating contig traversal after kmer {}"
                                 " because max node limit is reached").format(start_kmer))
        return self

    def traverse_from_each_kmer_in_iterable(self, iterable):
        for kmer in iterable:
            self._traverse_from(kmer)
        self._post_process_graph()
        return self

    def traverse_from(self, start_string):
        self._traverse_from(start_string)
        self._post_process_graph()
        return self

    def _traverse_from(self, start_string):
        assert len(start_string) == self.ra_parser.kmer_size
        self.branch_traverser = {
            color: branch.Traverser(self.ra_parser,
                                    traversal_color=color,
                                    other_stopping_colors=set(self.traversal_colors) - {color})
            for color in self.traversal_colors
        }
        self.queuer = branch.Queuer(self.branch_queue,
                                    traversal_colors=self.traversal_colors,
                                    engine_orientation=self.orientation)

        self._process_initial_branch(start_string)
        while 0 < len(self.branch_queue) and (
            self.max_nodes is None or len(self.graph) < self.max_nodes
        ):
            self._traverse_a_branch_from_queue()
        if self.max_nodes and len(self.graph) > self.max_nodes:
            raise Exception("Max nodes ({}) exceeded: {} nodes found".format(self.max_nodes,
                                                                             len(self.graph)))
        return self

    def _post_process_graph(self):
        self.graph = annotate_kmer_graph_edges(self.graph)

    def _add_graph_metadata(self):
        self.graph.graph['colors'] = self.ra_parser.colors
        self.graph.graph['sample_names'] = [n.decode() for n in self.ra_parser.sample_names]

    def _process_initial_branch(self, start_string):
        if self.orientation == EngineTraversalOrientation.both:
            first_traversal_orientation = EdgeTraversalOrientation.original
        else:
            first_traversal_orientation = EdgeTraversalOrientation[self.orientation.name]
        self.queuer.add_from(start_string,
                             first_traversal_orientation,
                             connecting_node=None,
                             traversal_color=self.traversal_colors[0])
        self._traverse_a_branch_from_queue()
        start_kmer = self.ra_parser.get_kmer_for_string(start_string)
        if self.orientation == EngineTraversalOrientation.both:
            for color in self.traversal_colors:
                oriented_edge_set = start_kmer.edges[color].oriented(
                    EdgeTraversalOrientation.reverse)
                kmer_strings = list(oriented_edge_set.neighbor_kmer_strings(start_string))
                if len(kmer_strings) == 1:
                    self.queuer.add_from(start_string=kmer_strings[0],
                                         orientation=EdgeTraversalOrientation.reverse,
                                         connecting_node=start_string,
                                         traversal_color=color)
        for color in self.traversal_colors[1:]:
            oriented_edge_set = start_kmer.edges[color].oriented(first_traversal_orientation)
            for kmer_string in oriented_edge_set.neighbor_kmer_strings(start_string):
                self.queuer.add_from(start_string=kmer_string,
                                     orientation=first_traversal_orientation,
                                     connecting_node=start_string,
                                     traversal_color=color)

    def _traverse_a_branch_from_queue(self):
        setup = self.branch_queue.popleft()
        color_branch_traverser = self.branch_traverser[setup.traversal_color]
        branch = color_branch_traverser.traverse_from(setup.start_string,
                                                      orientation=setup.orientation,
                                                      parent_graph=self.graph)
        add_graph_to(self.graph, branch.graph)
        self._connect_branch_to_parent_graph(branch, setup)
        self._link_branch_and_queue_neighbor_traversals(branch)

    def _connect_branch_to_parent_graph(self, branch, setup):
        if setup.connecting_node is not None and not branch.is_empty():
            self._add_edge_in_orientation(setup.connecting_node, branch.first_kmer_string,
                                          setup.orientation)

    def _link_branch_and_queue_neighbor_traversals(self, branch):
        branch = self._add_neighbors_from_other_colors_to_branch(branch)
        orientations_and_kmer_strings = [(branch.orientation, branch.neighbor_kmer_strings)]
        if self.orientation == EngineTraversalOrientation.both:
            orientations_and_kmer_strings.append(
                (EdgeTraversalOrientation.other(branch.orientation),
                 branch.reverse_neighbor_kmer_strings)
            )
        for orientation, kmer_strings in orientations_and_kmer_strings:
            for neighbor_string in kmer_strings:
                if neighbor_string in self.graph:
                    self._add_edge_in_orientation(branch.last_kmer_string,
                                                  neighbor_string,
                                                  orientation)
                else:
                    self.queuer.add_from_branch(branch)

    def _add_neighbors_from_other_colors_to_branch(self, branch):
        if branch.is_empty():
            return branch
        last_kmer_string = branch.last_kmer_string
        last_kmer = self.ra_parser.get_kmer_for_string(last_kmer_string)
        neighbor_kmer_strings = set(branch.neighbor_kmer_strings)
        reverse_neighbor_kmer_strings = set(branch.reverse_neighbor_kmer_strings)
        for traversal_color in self.traversal_colors:
            oriented_edge_set = last_kmer.edges[traversal_color].oriented(branch.orientation)
            neighbor_kmer_strings |= set(oriented_edge_set.neighbor_kmer_strings(last_kmer_string))
            reverse_neighbor_kmer_strings |= set(
                oriented_edge_set.other_orientation().neighbor_kmer_strings(last_kmer_string))

        branch = copy.copy(branch)
        branch.neighbor_kmer_strings = list(neighbor_kmer_strings)
        branch.reverse_neighbor_kmer_strings = list(reverse_neighbor_kmer_strings)
        return branch

    def _add_edge_in_orientation(self, kmer1_string, kmer2_string, orientation):
        if orientation == EdgeTraversalOrientation.reverse:
            kmer1_string, kmer2_string = kmer2_string, kmer1_string
        kmer1 = self.ra_parser.get_kmer_for_string(kmer1_string)
        kmer2 = self.ra_parser.get_kmer_for_string(kmer2_string)
        interactor = Interactor.from_graph(self.graph)
        interactor.add_edge_to_graph_for_kmer_pair(kmer1, kmer2, kmer1_string, kmer2_string,
                                                   kmer1.colors)

    def log_graph_size(self):
        if len(self.graph) > self.last_graph_size:
            self.last_graph_size = len(self.graph)
            self.logger.info('current graph size: {}'.format(self.last_graph_size))


def annotate_kmer_graph_edges(graph):
    """Adds nodes to graph for kmer_strings that only exist as edges in a node's kmer."""
    colors = graph.graph['colors']
    kmer_builder = EmptyKmerBuilder(num_colors=len(colors), default_coverage=1)
    for kmer_string, kmer in list(graph.nodes(data='kmer')):
        is_lexlo = bool(kmer_string == lexlo(kmer_string))
        for color in colors:
            for new_kmer_string in kmer.edges[color].get_outgoing_kmer_strings(kmer_string,
                                                                               is_lexlo=is_lexlo):
                if new_kmer_string not in graph.nodes:
                    graph.add_node(new_kmer_string, kmer=kmer_builder.build_or_get(new_kmer_string))
                    graph.add_edge(kmer_string, new_kmer_string, key=color)
            for new_kmer_string in kmer.edges[color].get_incoming_kmer_strings(kmer_string,
                                                                               is_lexlo=is_lexlo):
                if new_kmer_string not in graph.nodes:
                    graph.add_node(new_kmer_string, kmer=kmer_builder.build_or_get(new_kmer_string))
                    graph.add_edge(new_kmer_string, kmer_string, key=color)
    return graph

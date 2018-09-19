import itertools
from collections.abc import Collection, Mapping, MutableMapping

import attr
import networkx as nx

from cortexpy.constants import EdgeTraversalOrientation
from cortexpy.graph.parser.kmer import find_all_neighbors, disconnect_kmers
from cortexpy.utils import lexlo


def build_cortex_graph_from_header(header, **kwargs):
    return build_cortex_graph(sample_names=header.sample_names,
                              kmer_size=header.kmer_size,
                              num_colors=header.num_colors,
                              colors=header.colors,
                              **kwargs)


def build_empty_cortex_graph_from_ra_parser(ra_parser):
    return build_cortex_graph(sample_names=ra_parser.sample_names,
                              kmer_size=ra_parser.kmer_size,
                              num_colors=ra_parser.num_colors,
                              colors=ra_parser.colors)


def build_cortex_graph(*, sample_names, kmer_size, num_colors, colors, kmer_generator=None,
                       kmer_mapping=None):
    """Colored de Bruijn graph constructor"""
    assert (kmer_generator is None) or (kmer_mapping is None)
    if kmer_generator is not None:
        graph = CortexDiGraph({k.kmer: k for k in kmer_generator})
    elif kmer_mapping is not None:
        graph = CortexDiGraph(kmer_mapping)
    else:
        graph = CortexDiGraph()
    graph.graph['sample_names'] = sample_names
    graph.graph['kmer_size'] = kmer_size
    graph.graph['num_colors'] = num_colors
    graph.graph['colors'] = colors
    return graph


@attr.s(slots=True)
class CortexGraphMapping(MutableMapping):
    """Create a dict-like kmer mapping from a RandomAccess parser (ra_parser)

    The exclusion set tracks kmers deleted from the ra_parser.
    The new_kmers track kmers that have been added to the mapping.
    Kmers that exist in both new_kmers and ra_parser are considered overwritten. The kmers in
    new_kmers have precedence.
    """
    ra_parser = attr.ib()
    _exclusion_set = attr.ib(attr.Factory(set))
    _new_kmers = attr.ib(attr.Factory(dict))
    _n_duplicates = attr.ib(0)

    def __attrs_post_init__(self):
        if isinstance(self.ra_parser, type(self)):
            self._exclusion_set = self.ra_parser._exclusion_set
            self.ra_parser = self.ra_parser.ra_parser

    def __getitem__(self, key):
        lexlo_key = lexlo(key)
        if lexlo_key in self._exclusion_set:
            raise KeyError
        if lexlo_key in self._new_kmers:
            return self._new_kmers[lexlo_key]
        return self.ra_parser[lexlo_key]

    def __setitem__(self, key, value):
        lexlo_key = lexlo(key)
        if lexlo_key in self._exclusion_set:
            self._exclusion_set.discard(lexlo_key)
        if lexlo_key in self.ra_parser and lexlo_key not in self._new_kmers:
            self._n_duplicates += 1
        self._new_kmers[lexlo_key] = value

    def __delitem__(self, item):
        lexlo_string = lexlo(item)
        if lexlo_string in self._exclusion_set:
            raise KeyError
        in_new_kmers = lexlo_string in self._new_kmers
        in_ra_parser = lexlo_string in self.ra_parser
        if in_new_kmers:
            del self._new_kmers[lexlo_string]
        if in_ra_parser:
            self._exclusion_set.add(lexlo_string)
        if in_new_kmers and in_ra_parser:
            self._n_duplicates -= 1

    def __iter__(self):
        for kmer_string in self._new_kmers:
            yield kmer_string
        for kmer_string in self.ra_parser:
            if kmer_string not in self._exclusion_set and kmer_string not in self._new_kmers:
                yield kmer_string

    def __len__(self):
        return len(self.ra_parser) + len(self._new_kmers) - len(
            self._exclusion_set) - self._n_duplicates

    def disconnect_kmers(self, first, second, colors):
        """Disconnect two kmers"""
        are_neighbors = False
        for ref_kmer, flip_kmer, ref_letter, flip_letter in find_all_neighbors(first, second):
            are_neighbors = True
            if colors:
                self[ref_kmer.kmer] = ref_kmer
                self[flip_kmer.kmer] = flip_kmer
            for color in colors:
                ref_kmer.edges[color].remove_edge(ref_letter)
                flip_kmer.edges[color].remove_edge(flip_letter)
        if not are_neighbors:
            raise ValueError(
                'first kmer ({}) cannot be connected to second kmer ({})'.format(first.kmer,
                                                                                 second.kmer)
            )

    def connect_kmers(self, first, second, color, identical_kmer_check=True):
        """Connect two kmers"""
        if identical_kmer_check and first == second and first is not second:
            raise ValueError('Kmers are equal, but not the same object')
        are_neighbors = False
        for ref_kmer, flip_kmer, ref_letter, flip_letter in find_all_neighbors(first, second):
            are_neighbors = True
            self[ref_kmer.kmer] = ref_kmer
            self[flip_kmer.kmer] = flip_kmer
            ref_kmer.edges[color].add_edge(ref_letter)
            flip_kmer.edges[color].add_edge(flip_letter)
        if not are_neighbors:
            raise ValueError(
                'first kmer ({}) cannot be connected to second kmer ({})'.format(first.kmer,
                                                                                 second.kmer)
            )


@attr.s(slots=True)
class CortexDiGraph(Collection):
    """Stores cortex k-mers and conforms to parts of the interface of networkx.MultiDiGraph"""
    _kmer_mapping = attr.ib(attr.Factory(dict))
    graph = attr.ib(attr.Factory(dict))

    def __attrs_post_init__(self):
        # todo: implement .from_ra_parser class method as described in
        # http://www.attrs.org/en/stable/init.html
        if isinstance(self._kmer_mapping, type(self)):
            self.graph = self._kmer_mapping.graph
            self._kmer_mapping = self._kmer_mapping._kmer_mapping
        else:
            self._kmer_mapping = CortexGraphMapping(self._kmer_mapping)

    @property
    def nodes(self):
        return NodeView(self._kmer_mapping)

    @property
    def node(self):
        return NodeView(self._kmer_mapping)

    def __len__(self):
        return len(self._kmer_mapping)

    def __iter__(self):
        return iter(self._kmer_mapping)

    def __contains__(self, item):
        return item in self._kmer_mapping.keys()

    def __getitem__(self, item):
        return self.succ[item]

    def add_node(self, kmer_string, *, kmer):
        self._kmer_mapping[kmer_string] = kmer

    def add_nodes_from(self, node_iterable):
        for node in node_iterable:
            self.add_node(node[0], kmer=node[1])

    def is_multigraph(self):
        return True

    def is_directed(self):
        return True

    def is_consistent(self):
        return False

    @property
    def edges(self):
        return EdgeView(self)

    def out_edges(self, node, keys=False, default=None, data=None):
        kmer = self.node[node]
        is_lexlo = kmer.kmer == node
        for color in kmer.colors:
            for out_node in kmer.edges[color].get_outgoing_kmer_strings(node, is_lexlo=is_lexlo):
                if keys:
                    yield (node, out_node, color)
                else:
                    yield (node, out_node)

    def in_edges(self, node, keys=False, default=None, data=None):
        kmer = self.node[node]
        is_lexlo = kmer.kmer == node
        for color in kmer.colors:
            for in_node in kmer.edges[color].get_incoming_kmer_strings(node, is_lexlo=is_lexlo):
                if keys:
                    yield (in_node, node, color)
                else:
                    yield (in_node, node)

    def out_degree(self, node):
        return len(list(self.out_edges(node)))

    def in_degree(self, node):
        return len(list(self.in_edges(node)))

    def add_edge(self, first, second, *, key):
        """Note: edges can only be added to existing nodes"""
        first_kmer = self.node[first]
        second_kmer = self.node[second]
        self._kmer_mapping.connect_kmers(first_kmer, second_kmer, color=key)

    @property
    def succ(self):
        return MultiAdjacencyView(self._kmer_mapping, EdgeTraversalOrientation.original)

    @property
    def pred(self):
        return MultiAdjacencyView(self._kmer_mapping, EdgeTraversalOrientation.reverse)

    def remove_node(self, node):
        try:
            node_kmer = self.node[node]
        except KeyError:
            return
        for succ in self.succ[node]:
            succ_kmer = self.node[succ]
            self._kmer_mapping.disconnect_kmers(node_kmer, succ_kmer, node_kmer.colors)
        for pred in self.pred[node]:
            pred_kmer = self.node[pred]
            self._kmer_mapping.disconnect_kmers(node_kmer, pred_kmer, node_kmer.colors)
        del self._kmer_mapping[node]

    def remove_nodes_from(self, nodes):
        for node in nodes:
            self.remove_node(node)

    def nbunch_iter(self, nbunch=None):
        """Return an iterator over nodes contained in nbunch that are
        also in the graph.

        The nodes in nbunch are checked for membership in the graph
        and if not are silently ignored.

        Parameters
        ----------
        nbunch : single node, container, or all nodes (default= all nodes)
            The view will only report edges incident to these nodes.

        Returns
        -------
        niter : iterator
            An iterator over nodes in nbunch that are also in the graph.
            If nbunch is None, iterate over all nodes in the graph.

        Raises
        ------
        NetworkXError
            If nbunch is not a node or or sequence of nodes.
            If a node in nbunch is not hashable.

        See Also
        --------
        Graph.__iter__

        Notes
        -----
        When nbunch is an iterator, the returned iterator yields values
        directly from nbunch, becoming exhausted when nbunch is exhausted.

        To test whether nbunch is a single node, one can use
        "if nbunch in self:", even after processing with this routine.

        If nbunch is not a node or a (possibly empty) sequence/iterator
        or None, a :exc:`NetworkXError` is raised.  Also, if any object in
        nbunch is not hashable, a :exc:`NetworkXError` is raised.

        Licence
        -------
        This method was copied from Networkx version 2.1 and then modified
        """
        if nbunch is None:  # include all nodes via iterator
            bunch = iter(self._adj)
        elif nbunch in self:  # if nbunch is a single node
            bunch = iter([nbunch])
        else:  # if nbunch is a sequence of nodes
            def bunch_iter(nlist, adj):
                try:
                    for n in nlist:
                        if n in adj:
                            yield n
                except TypeError as e:
                    message = e.args[0]
                    # capture error for non-sequence/iterator nbunch.
                    if 'iter' in message:
                        msg = "nbunch is not a node or a sequence of nodes."
                        raise nx.NetworkXError(msg)
                    # capture error for unhashable node.
                    elif 'hashable' in message:
                        msg = "Node {} in sequence nbunch is not a valid node."
                        raise nx.NetworkXError(msg.format(n))
                    else:
                        raise

            bunch = bunch_iter(nbunch, self._adj)
        return bunch


@attr.s(slots=True)
class ConsistentCortexDiGraph(Collection):
    """Graph that stores kmer strings that are consistent with each other"""
    _kmer_mapping = attr.ib(attr.Factory(dict))
    graph = attr.ib(attr.Factory(dict))  # refers to graph attribute of nx.Graph

    def __attrs_post_init__(self):
        if isinstance(self._kmer_mapping, (CortexDiGraph, ConsistentCortexDiGraph)):
            self.graph = self._kmer_mapping.graph
        if isinstance(self._kmer_mapping, ConsistentCortexDiGraph):
            self._kmer_mapping = self._kmer_mapping._kmer_mapping

    def is_consistent(self):
        return True

    def is_directed(self):
        return True

    def is_multigraph(self):
        return True

    def __contains__(self, item):
        return item in self._kmer_mapping

    def __iter__(self):
        yield from self._kmer_mapping

    def __len__(self):
        return len(self._kmer_mapping)

    def __getitem__(self, item):
        return self.succ[item]

    @property
    def succ(self):
        return MultiAdjacencyView(self._kmer_mapping,
                                  EdgeTraversalOrientation.original,
                                  return_lexlo_kmers=False)

    @property
    def pred(self):
        return MultiAdjacencyView(self._kmer_mapping,
                                  EdgeTraversalOrientation.reverse,
                                  return_lexlo_kmers=False)

    @property
    def nodes(self):
        return NodeView(self._kmer_mapping)

    @property
    def node(self):
        return NodeView(self._kmer_mapping)

    def add_node(self, kmer_string, *, kmer):
        self._kmer_mapping[kmer_string] = kmer

    def remove_edge(self, k, v, key):
        disconnect_kmers(self._kmer_mapping[k], self._kmer_mapping[v], [key])

    @property
    def edges(self):
        return EdgeView(self)

    def out_degree(self, node):
        return len(list(self.out_edges(node)))

    def in_degree(self, node):
        return len(list(self.in_edges(node)))

    def out_edges(self, node, keys=False, default=None, data=None):
        kmer = self._kmer_mapping[node]
        is_lexlo = kmer.kmer == node
        for color in kmer.colors:
            for out_node in kmer.edges[color].get_outgoing_kmer_strings(node, is_lexlo=is_lexlo):
                if keys:
                    yield (node, out_node, color)
                else:
                    yield (node, out_node)

    def in_edges(self, node, keys=False, default=None, data=None):
        kmer = self.node[node]
        is_lexlo = kmer.kmer == node
        for color in kmer.colors:
            for in_node in kmer.edges[color].get_incoming_kmer_strings(node, is_lexlo=is_lexlo):
                if keys:
                    yield (in_node, node, color)
                else:
                    yield (in_node, node)


@attr.s(slots=True)
class NodeView(Collection):
    _nodes = attr.ib()

    def __call__(self, *args, data=False, **kwargs):
        if data:
            yield from self._nodes.items()
        else:
            yield from self._nodes.keys()

    def __len__(self):
        return len(self._nodes)

    def __iter__(self):
        return iter(self._nodes)

    def __contains__(self, item):
        return item in self._nodes

    def __getitem__(self, item):
        return self._nodes[item]


def get_canonical_edge(first, second):
    """Get canonical edge.

    Canonical edges are between lexlo kmers and are ordered lexicographically

    Return canonical edge, if the first and second nodes were lexlo"""
    lexlo_first = lexlo(first)
    lexlo_second = lexlo(second)
    flip_first_second = lexlo_second < lexlo_first
    if flip_first_second:
        lexlo_second, lexlo_first = lexlo_first, lexlo_second
    return lexlo_first, lexlo_second, flip_first_second


@attr.s(slots=True)
class EdgeView(object):
    graph = attr.ib()

    def __call__(self, *args, data=False, keys=False, **kwargs):
        if not args:
            if self.graph.is_consistent():
                yield from self._edge_iter_consistent(data=data, keys=keys)
            else:
                yield from self._edge_iter(data=data, keys=keys)
        else:
            if self.graph.is_directed():
                yield from self.graph.out_edges(args[0], data=data, keys=keys)
            else:
                yield from self._edge_iter(data=data, keys=keys)

    def __len__(self):
        return len(list(self()))

    def __iter__(self):
        return self()

    def __contains__(self, item):
        return item in self()

    def _edge_iter(self, data=False, keys=False):
        """Non-consistent edges of graph. Only returns canonical edges"""
        seen_self_edges = set()
        for node, kmer in self.graph.nodes(data=True):
            for color in kmer.colors:
                edge = kmer.edges[color]
                neighbors = itertools.chain(
                    edge.get_outgoing_kmers(kmer.kmer),
                    edge.get_incoming_kmers(kmer.kmer)
                )
                for neighbor in neighbors:
                    if neighbor not in self.graph:
                        continue
                    if neighbor < kmer.kmer:
                        continue
                    if neighbor == kmer.kmer:
                        if neighbor in seen_self_edges:
                            continue
                        else:
                            seen_self_edges.add(neighbor)

                    ret = [node, neighbor]
                    if keys:
                        ret.append(color)
                    if data:
                        ret.append({})
                    yield tuple(ret)

    def _edge_iter_consistent(self, data=False, keys=False):
        for node, kmer in self.graph.nodes(data=True):
            is_lexlo = node == kmer.kmer
            for color in kmer.colors:
                for kmer_string in kmer.edges[color].get_outgoing_kmer_strings(node,
                                                                               is_lexlo=is_lexlo):
                    if kmer_string in self.graph.nodes:
                        ret = [node, kmer_string]
                        if keys:
                            ret.append(color)
                        if data:
                            ret.append({})
                        yield tuple(ret)


@attr.s(slots=True)
class AdjancencyView(object):
    _nodes = attr.ib()
    kmer = attr.ib()
    query = attr.ib()
    orientation = attr.ib()
    return_lexlo_kmers = attr.ib(True)

    @property
    def colors(self):
        return self.kmer.colors

    def __iter__(self):
        edge_kmers = set()
        query_is_lexlo = self.kmer.kmer == self.query
        for color in self.colors:
            if self.orientation == EdgeTraversalOrientation.original:
                node_iter = self.kmer.edges[color].get_outgoing_kmer_strings(
                    self.query, is_lexlo=query_is_lexlo
                )
            else:
                node_iter = self.kmer.edges[color].get_incoming_kmer_strings(
                    self.query, is_lexlo=query_is_lexlo
                )
            for out_node in node_iter:
                if self.return_lexlo_kmers:
                    out_node = lexlo(out_node)
                edge_kmers.add(out_node)
        return iter(edge_kmers)

    def __getitem__(self, item):
        other = self._nodes[item]
        edge_colors = set()
        if self.orientation == EdgeTraversalOrientation.original:
            edge_func = self.kmer.has_outgoing_edge_to_kmer_in_color
        else:
            edge_func = self.kmer.has_incoming_edge_from_kmer_in_color
        try:
            for color in self.colors:
                if edge_func(other, color):
                    edge_colors.add(color)
        except ValueError:
            pass
        return edge_colors


@attr.s(slots=True)
class MultiAdjacencyView(object):
    _nodes = attr.ib()
    orientation = attr.ib()
    return_lexlo_kmers = attr.ib(True)

    def __getitem__(self, item):
        kmer = self._nodes[item]
        return AdjancencyView(self._nodes,
                              kmer,
                              query=item,
                              orientation=self.orientation,
                              return_lexlo_kmers=self.return_lexlo_kmers)


@attr.s(slots=True)
class DictView(Mapping):
    base_dict = attr.ib()
    allowed_keys = attr.ib(None)

    def __getitem__(self, item):
        if item in self.allowed_keys:
            return self.base_dict[item]

    def __iter__(self):
        for key, val in self.base_dict.items():
            if key in self.allowed_keys:
                yield key, val

    def __len__(self):
        num_overlap = len(self.allowed_keys & self.base_dict.keys())
        return len(self.base_dict) - num_overlap

    def __copy__(self):
        return {k: v for k, v in self}

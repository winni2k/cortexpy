import json
import attr
import networkx as nx
from networkx.readwrite import json_graph


@attr.s(slots=True)
class Serializer(object):
    """A class for serializing kmer graphs."""
    graph = attr.ib()
    collapse_kmer_unitigs = attr.ib(False)

    def to_json_serializable(self):
        graph = self.graph.copy()
        if self.collapse_kmer_unitigs:
            graph = collapse_kmer_unitigs(graph)
        return make_graph_json_representable(graph)

    def to_json(self):
        graph = self.to_json_serializable()
        serializable = json_graph.node_link_data(graph, attrs={'link': 'edges'})
        return json.dumps(serializable)


def make_graph_json_representable(graph):
    graph = graph.copy()
    for node, node_data in graph.nodes.items():
        kmer = node_data.pop('kmer', None)
        if kmer is None:
            node_data['is_missing'] = True
        else:
            node_data['coverage'] = list(kmer.coverage)
        if isinstance(node_data.get('node_object'), nx.Graph):
            node_data['node_object'] = repr(node_data['node_object'])
    return graph


def traverse_to_next_junction(graph, start_node, unvisited_nodes, orientation):
    path = nx.DiGraph()
    path.add_node(start_node)
    for source, target in nx.edge_dfs(graph, start_node, orientation=orientation):
        if target in unvisited_nodes:
            path.add_edge(source, target)
            if graph.degree[target] > 2:
                break
        else:
            break
    return path


def find_unitig_from(start_node, graph):
    """Find a unitig that contains :param start_node:

    # from https://github.com/mcveanlab/mccortex/wiki/unitig
    Unitig definition: A maximal run of de Bruijn nodes (kmers, each with an orientation)
    where all nodes must have in-degree and out-degree exactly 1, with the exception of in-degree
    of the first node and out-degree of the last node (they may have any number).
    Additionally, a unitig must not visit a kmer more than once.
    """
    unitig_graph = nx.DiGraph()
    left_node = right_node = start_node
    seen = {start_node}
    for source, target in nx.edge_dfs(graph, start_node, orientation='forward'):
        if graph.in_degree(target) != 1 or target in seen:
            break
        right_node = target
        unitig_graph.add_edge(source, target)
        seen.add(target)
        if graph.out_degree(target) != 1:
            break
    for source, target, orientation in nx.edge_dfs(graph, start_node, orientation='reverse'):
        if graph.out_degree(source) != 1 or source in seen:
            break
        left_node = source
        unitig_graph.add_edge(source, target)
        seen.add(source)
        if graph.in_degree(source) != 1:
            break
    return unitig_graph, left_node, right_node


def find_unitigs(graph, in_place=False):
    """Finds unitigs of 2 or more nodes and replaces them by a subgraph containing the
    unitig nodes.

    Unitigs are described in :py:func:`find_unitig_from`.

    returns a graph containing unitig subgraphs
    """

    unitigable_nodes = set()
    if not in_place:
        graph = graph.copy()
    for node in graph:
        if graph.in_degree(node) < 2 and graph.out_degree(node) < 2:
            unitigable_nodes.add(node)

    visited_nodes = set()
    for start_node in unitigable_nodes:
        if start_node in visited_nodes:
            continue
        unitig_graph, left_node, right_node = find_unitig_from(start_node, graph)
        unitig_graph_set = set(unitig_graph)
        assert visited_nodes & unitig_graph_set == set()
        visited_nodes |= unitig_graph_set
        if len(unitig_graph) > 1:
            graph.add_node(unitig_graph, is_unitig=True, left_node=left_node, right_node=right_node)
            for source, _ in graph.in_edges(left_node):
                graph.add_edge(source, unitig_graph)
                if source == right_node:
                    graph.add_edge(unitig_graph, unitig_graph)
            for _, target in graph.out_edges(right_node):
                graph.add_edge(unitig_graph, target)
                if target == left_node:
                    graph.add_edge(unitig_graph, unitig_graph)
            graph.remove_nodes_from(unitig_graph)
    return graph


def collapse_kmer_unitigs(graph):
    """Collapses unitig kmers into a single graph node.

    All nodes have an attribute `name` added which is the node label with the left k-1 letters
    removed unless the node has no incoming edges.

    Unitigs are described in :py:func:`find_unitig_from`.

    :param graph:
    :return graph:
    """
    unitig_graph = find_unitigs(graph)
    out_graph = unitig_graph.copy()
    for node, data in unitig_graph.nodes.data():
        if data.get('is_unitig'):
            left_node = data['left_node']
            if graph.in_degree(left_node) > 0:
                short_kmer_name = [left_node[-1]]
            else:
                short_kmer_name = [left_node]
            for source, target in nx.edge_dfs(node, left_node):
                short_kmer_name.append(target[-1])
            short_kmer_name = ''.join(short_kmer_name)
        else:
            if graph.in_degree(node) == 0:
                short_kmer_name = node
            else:
                short_kmer_name = node[-1]
        out_graph.nodes[node]['repr'] = short_kmer_name
    out_graph = nx.convert_node_labels_to_integers(out_graph, label_attribute='node_object')
    return out_graph

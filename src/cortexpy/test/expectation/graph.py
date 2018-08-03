import logging

import attr

logger = logging.getLogger(__name__)


@attr.s(slots=True)
class KmerGraphsExpectation(object):
    graph_list = attr.ib()

    def has_n_graphs(self, n):
        assert n == len(self.graph_list)
        return self

    def has_nodes(self, *expected_nodes):
        expected_nodes = set(expected_nodes)
        nodes = set()
        for graph in self.graph_list:
            nodes |= graph.nodes
        assert expected_nodes == nodes
        return self


@attr.s(slots=True)
class KmerGraphExpectation(object):
    graph = attr.ib()
    sort_edges = attr.ib(False)

    def __attrs_post_init__(self):
        # logger.info(nx.node_link_data(self.graph))
        pass

    def has_node_coverages(self, *node_coverage_strings):
        for node_coverage_string in node_coverage_strings:
            self.has_node_coverage(node_coverage_string)
        self.has_n_nodes(len(node_coverage_strings))

    def has_node_coverage(self, node_coverage_string):
        node, *coverages = node_coverage_string.split()
        coverages = [int(c) for c in coverages]
        self.has_node(node).has_coverages(*coverages)
        return self

    def has_node(self, node):
        assert node in self.graph
        return KmerNodeExpectation(self.graph.node[node])

    def has_n_nodes(self, n):
        assert n == len(self.graph)
        return self

    def has_n_edges(self, n):
        assert n == len(self.graph.edges)
        return self

    def has_nodes(self, *expected_nodes):
        assert sorted(expected_nodes) == sorted(self.graph.nodes)
        return self

    def has_edges(self, *edges):
        expected_edges = []
        for edge in edges:
            if isinstance(edge, str):
                edge = edge.split(' ')
                edge[2] = int(edge[2])
            if self.sort_edges:
                edge[0:2] = sorted(edge[0:2])
            expected_edges.append(tuple(edge))

        assert sorted(expected_edges) == sorted(self.graph.edges(keys=True))
        return self

    def has_edge(self, source, target, color):
        assert (source, target, color) in self.graph.edges
        return self


@attr.s
class KmerNodeExpectation(object):
    kmer_node = attr.ib()

    def has_coverages(self, *coverages):
        if len(coverages) == 1 and isinstance(coverages[0], str):
            coverages = tuple([int(c) for c in coverages[0].split(' ')])
        assert (coverages) == self.kmer_node['kmer'].coverage
        return self

    def is_missing(self):
        print(self.kmer_node['kmer'].coverage)
        assert all(c == 0 for c in self.kmer_node['kmer'].coverage)
        return self

    def is_not_missing(self):
        assert not self.kmer_node.get('is_missing', False)
        return self

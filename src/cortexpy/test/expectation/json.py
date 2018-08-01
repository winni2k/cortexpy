import json
import logging
import sys
from collections import defaultdict

import attr

logger = logging.getLogger(__name__)


@attr.s(slots=True)
class JsonGraphs(object):
    graph_file_names = attr.ib()
    graphs = attr.ib(attr.Factory(list))

    def __attrs_post_init__(self):
        for file_name in self.graph_file_names:
            with open(str(file_name)) as fh:
                self.graphs.append(json.load(fh))

    def has_graph_with_node_reprs(self, *reprs):
        reprs = set(reprs)
        graphs_with_node_reprs = [g for g in self.graphs if
                                  reprs <= {n['repr'] for n in g['nodes']}]
        assert 0 < len(graphs_with_node_reprs)
        return self

    def has_n_graphs(self, n):
        assert n == len(self.graphs)
        return self


@attr.s(slots=True)
class JsonGraph(object):
    graph = attr.ib()
    colors = attr.ib(init=False)
    n_colors = attr.ib(init=False)
    node_id_by_repr = attr.ib(attr.Factory(lambda: defaultdict(list)), init=False)
    node_id_by_unitig = attr.ib(attr.Factory(lambda: defaultdict(list)), init=False)

    def __attrs_post_init__(self):
        logger.info('With JSON graph: {}'.format(json.dumps(self.graph, indent=2, sort_keys=True)))
        assert sum(['is_missing' in n for n in self.graph['nodes']]) == 0
        assert sum(['is_missing' in e for e in self.graph['edges']]) == 0
        self.colors = self.graph['graph']['colors']
        self.n_colors = len(self.colors)
        for node_id, node in enumerate(self.graph['nodes']):
            self.node_id_by_repr[node['repr']].append(node_id)
            self.node_id_by_unitig[node['unitig']].append(node_id)

    @classmethod
    def from_string(cls, string):
        return cls(json.loads(string))

    def is_directed(self):
        assert self.graph['directed']
        return self

    def has_colors(self, colors):
        assert self.colors == colors
        return self

    def has_n_nodes(self, n):
        assert len(self.graph['nodes']) == n
        return self

    def has_node_unitig(self, unitig):
        node_ids = self.node_id_by_unitig[unitig]
        assert len(node_ids) == 1
        node = self.graph['nodes'][node_ids[0]]
        return JsonNodeExpectation(node, self.n_colors)

    def has_node_reprs(self, *reprs):
        assert set(reprs) == set(self.node_id_by_repr.keys())
        return self

    def has_node_repr(self, repr):
        node_ids = self.node_id_by_repr[repr]
        assert len(node_ids) > 0
        if len(node_ids) != 1:
            raise AssertionError('Found more than one node matching repr: ' + repr)
        node = self.graph['nodes'][node_ids[0]]
        return JsonNodeExpectation(node, self.n_colors)

    def has_n_edges(self, n):
        assert len(self.graph['edges']) == n
        return self

    def has_repr_edge(self, source_repr, target_repr, color):
        source_id_list = self.node_id_by_repr[source_repr]
        target_id_list = self.node_id_by_repr[target_repr]
        assert len(source_id_list) == 1
        assert len(target_id_list) == 1
        print('source_repr={}, target_repr={}'.format(source_repr, target_repr))
        return self.has_edge(source_id_list[0], target_id_list[0], color)

    def has_edge(self, source, target, color):
        num_matching_edges = 0
        matching_edge = None
        for e in self.graph['edges']:
            edge = (e['source'], e['target'], e['key'])
            if edge == (source, target, color):
                num_matching_edges += 1
                matching_edge = edge
        assert matching_edge == (source, target, color)
        assert num_matching_edges == 1
        return self


@attr.s(slots=True)
class JsonNodeExpectation(object):
    node = attr.ib()
    n_colors = attr.ib(None)

    def has_coverages_by_color(self, *coverages):
        coverages = list(listify_elements(coverages))
        assert coverages == self.node['coverage']
        return self

    def has_coverages_by_kmer(self, *coverages):
        coverages = list(listify_elements(coverages))
        return self.has_coverages_by_color(*[list(i) for i in zip(*coverages)])


def expect_zero_return_code(completed_process):
    stdout = completed_process.stdout
    if completed_process.returncode != 0:
        print(stdout)
        print(completed_process.stderr, file=sys.stderr)
        assert completed_process.returncode == 0


def listify_elements(iterable):
    return ([[e], e][int(isinstance(e, list))] for e in iterable)


@attr.s(slots=True)
class TraversalExpectation(object):
    completed_process = attr.ib()

    def __attrs_post_init__(self):
        expect_zero_return_code(self.completed_process)

    def returns_json_graph(self):
        return JsonGraph(json.loads(self.completed_process.stdout))

import os
import json
import sys
from collections import defaultdict
import attr

from cortexpy.test import builder as builder, runner as runner

if os.environ.get('CI'):
    SPAWN_PROCESS = True
else:
    SPAWN_PROCESS = False


def expect_zero_return_code(completed_process):
    stdout = completed_process.stdout
    if completed_process.returncode != 0:
        print(stdout)
        print(completed_process.stderr, file=sys.stderr)
        assert completed_process.returncode == 0


def listify_elements(iterable):
    return ([[e], e][int(isinstance(e, list))] for e in iterable)


@attr.s(slots=True)
class JsonNodeExpectation(object):
    node = attr.ib()
    n_colors = attr.ib(None)

    def has_coverages(self, *coverages):
        coverages = list(listify_elements(coverages))
        for coverage in coverages:
            assert len(coverage) == self.n_colors
        assert self.node['coverage'] == coverages
        return self


@attr.s(slots=True)
class JsonGraphExpectation(object):
    graph = attr.ib()
    colors = attr.ib(init=False)
    n_colors = attr.ib(init=False)
    node_id_by_repr = attr.ib(init=False)

    def __attrs_post_init__(self):
        print('With JSON graph: {}'.format(self.graph))
        assert sum(['is_missing' in n for n in self.graph['nodes']]) == 0
        assert sum(['is_missing' in e for e in self.graph['edges']]) == 0
        self.colors = self.graph['graph']['colors']
        self.n_colors = len(self.colors)
        self.node_id_by_repr = defaultdict(list)
        for node_id, node in enumerate(self.graph['nodes']):
            self.node_id_by_repr[node['repr']].append(node_id)

    def is_directed(self):
        assert self.graph['directed']
        return self

    def has_colors(self, colors):
        assert self.colors == colors
        return self

    def has_n_nodes(self, n):
        assert len(self.graph['nodes']) == n
        return self

    def has_node_repr(self, repr):
        nodes = list(filter(lambda n: n['repr'] == repr, self.graph['nodes']))
        assert len(nodes) > 0
        return JsonNodeExpectation(nodes[0], self.n_colors)

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


class TestContig(object):
    def test_outputs_json(self, tmpdir):
        # given
        record = 'ACCTT'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        # when
        completed_process = (runner
                             .Cortexpy(SPAWN_PROCESS)
                             .view_contig(output_type='json', graph=output_graph, contig=record))
        stdout = completed_process.stdout

        # then
        assert completed_process.returncode == 0, completed_process
        expect = JsonGraphExpectation(json.loads(stdout))
        expect.is_directed()

    def test_collapse_kmer_unitigs_option(self, tmpdir):
        # given
        record1 = 'AAACCCGAA'
        record2 = 'ACCG'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record1)
                        .with_dna_sequence(record2)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))
        runner.Mccortex(kmer_size).view(output_graph)

        # when
        completed_process = (runner
                             .Cortexpy(SPAWN_PROCESS)
                             .view_contig(contig=record1,
                                          output_type='json',
                                          graph=output_graph))

        # then
        expect_zero_return_code(completed_process)

        stdout = completed_process.stdout
        expect = JsonGraphExpectation(json.loads(stdout))

        expect.has_colors([0, 1])
        expect.has_n_nodes(3)
        expect.has_node_repr('AAACC').has_coverages([1, 1], [1, 1], [2, 1])
        expect.has_node_repr('C').has_coverages([1, 1])
        expect.has_node_repr('GAA').has_coverages([2, 1], [1, 1], [1, 1])

        for edge in [('AAACC', 'C', 0), ('AAACC', 'C', 1), ('C', 'GAA', 0), ('C', 'GAA', 1),
                     ('AAACC', 'GAA', 0)]:
            expect.has_repr_edge(*edge)
        expect.has_n_edges(5)

    def test_collapse_kmer_unitigs_option_with_missing_kmers(self, tmpdir):
        # given
        record1 = 'AAACCCGAA'
        record2 = 'ACCG'
        query_record = record1 + 'G'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record1)
                        .with_dna_sequence(record2)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))
        runner.Mccortex(kmer_size).view(output_graph)

        # when
        completed_process = (runner
                             .Cortexpy(SPAWN_PROCESS)
                             .view_contig(contig=query_record,
                                          output_type='json',
                                          graph=output_graph))

        # then
        expect_zero_return_code(completed_process)

        stdout = completed_process.stdout
        expect = JsonGraphExpectation(json.loads(stdout))

        expect.has_n_nodes(4)
        expect.has_node_repr('AAACC').has_coverages([1, 1], [1, 1], [2, 1])
        expect.has_node_repr('C').has_coverages([1, 1])
        expect.has_node_repr('GAA').has_coverages([2, 1], [1, 1], [1, 1])
        expect.has_node_repr('G').has_coverages([0, 1])

        for color in [0, 1]:
            for edge in [['AAACC', 'C'], ['C', 'GAA']]:
                expect.has_repr_edge(edge[0], edge[1], color)
        expect.has_repr_edge('GAA', 'G', 1)
        expect.has_repr_edge('AAACC', 'GAA', 0)
        expect.has_n_edges(6)


class TestTraversal(object):
    def test_with_initial_kmers(self, tmpdir):
        # given
        record1 = 'AAACCCGAA'
        record2 = 'ACCG'
        kmer_size = 3
        for query_string_idx, query_string in enumerate(['CCC', record1, record2]):
            query_dir = tmpdir.join(str(query_string_idx))
            query_dir.ensure(dir=True)

            output_graph = (builder.Mccortex()
                            .with_dna_sequence(record1)
                            .with_dna_sequence(record2)
                            .with_kmer_size(kmer_size)
                            .build(query_dir))
            runner.Mccortex(kmer_size).view(output_graph)

            # when
            completed_process = (runner
                                 .Cortexpy(SPAWN_PROCESS)
                                 .view_traversal(contig=query_string, graph=output_graph, color=0))

            # then
            expect_zero_return_code(completed_process)

            stdout = completed_process.stdout
            expect = JsonGraphExpectation(json.loads(stdout))

            expect.has_n_nodes(3)
            expect.has_node_repr('C').has_coverages([1])
            expect.has_node_repr('AAACC').has_coverages([1], [1], [2])
            expect.has_node_repr('GAA').has_coverages([2], [1], [1])

            for color in [0]:
                for edge in [['AAACC', 'C'], ['C', 'GAA']]:
                    expect.has_repr_edge(edge[0], edge[1], color)
            expect.has_repr_edge('AAACC', 'GAA', 0)
            expect.has_n_edges(3)

    def test_with_peripheral_edges_creates_unitigs_(self, tmpdir):
        # given
        record1 = 'AAACCCGAA'
        record2 = 'ACCG'
        record3 = record1 + 'G'
        kmer_size = 3

        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record1)
                        .with_dna_sequence(record2)
                        .with_dna_sequence(record3, name='sample_1')
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))
        runner.Mccortex(kmer_size).view(output_graph)

        # when
        completed_process = (runner
                             .Cortexpy(SPAWN_PROCESS)
                             .view_traversal(contig='CCC', graph=output_graph, color=0))

        # then
        expect_zero_return_code(completed_process)

        stdout = completed_process.stdout
        expect = JsonGraphExpectation(json.loads(stdout))

        expect.has_node_repr('C').has_coverages([1, 1])
        expect.has_node_repr('AAACC').has_coverages([1, 1], [1, 1], [2, 1])
        expect.has_node_repr('GAA').has_coverages([2, 1], [1, 1], [1, 1])
        expect.has_node_repr('G').has_coverages([0, 0])
        expect.has_n_nodes(4)

        for color in [0, 1]:
            for edge in [['AAACC', 'C'], ['C', 'GAA']]:
                expect.has_repr_edge(edge[0], edge[1], color)
        expect.has_repr_edge('AAACC', 'GAA', 0)
        expect.has_repr_edge('GAA', 'G', 1)
        expect.has_n_edges(6)
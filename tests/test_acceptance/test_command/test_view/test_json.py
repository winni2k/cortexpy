import json
import os

from cortexpy.test import builder as builder, runner as runner
from cortexpy.test.expectation.json import JsonGraph, expect_zero_return_code

if os.environ.get('CI'):
    SPAWN_PROCESS = True
else:
    SPAWN_PROCESS = False


class TestContig(object):
    def test_outputs_json(self, tmpdir):
        # given
        record = 'ACCTT'
        kmer_size = 3
        output_graph = builder.Mccortex() \
            .with_dna_sequence(record) \
            .with_kmer_size(kmer_size) \
            .build(tmpdir)

        # when
        completed_process = runner \
            .Cortexpy(SPAWN_PROCESS) \
            .view_contig(to_json=True, graph=output_graph, contig=record)
        stdout = completed_process.stdout

        # then
        assert completed_process.returncode == 0, completed_process
        expect = JsonGraph(json.loads(stdout))
        expect.is_directed()

    def test_collapse_kmer_unitigs_option(self, tmpdir):
        # given
        record1 = 'AAACCCGAA'
        record2 = 'ACCG'
        record3 = 'TTCGGGTTT'
        kmer_size = 3
        output_graph = builder.Mccortex() \
            .with_dna_sequence(record1) \
            .with_dna_sequence(record2) \
            .with_dna_sequence(record3) \
            .with_kmer_size(kmer_size) \
            .build(tmpdir)
        runner.Mccortex(kmer_size).view(output_graph)

        # when
        completed_process = runner \
            .Cortexpy(SPAWN_PROCESS) \
            .view_contig(contig=record1,
                         to_json=True,
                         graph=output_graph)

        # then
        expect_zero_return_code(completed_process)

        stdout = completed_process.stdout
        expect = JsonGraph(json.loads(stdout))

        expect.has_colors([0, 1])
        expect.has_n_nodes(3)
        expect.has_node_repr('AAACC').has_coverages_by_kmer([2, 1], [2, 1], [3, 1])
        expect.has_node_repr('C').has_coverages_by_kmer([2, 1])
        expect.has_node_repr('GAA').has_coverages_by_kmer([3, 1], [2, 1], [2, 1])

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
        output_graph = builder.Mccortex() \
            .with_dna_sequence(record1) \
            .with_dna_sequence(record2) \
            .with_kmer_size(kmer_size) \
            .build(tmpdir)
        runner.Mccortex(kmer_size).view(output_graph)

        # when
        completed_process = runner \
            .Cortexpy(SPAWN_PROCESS) \
            .view_contig(contig=query_record,
                         to_json=True,
                         graph=output_graph)

        # then
        expect_zero_return_code(completed_process)

        stdout = completed_process.stdout
        expect = JsonGraph(json.loads(stdout))

        expect.has_n_nodes(4)
        expect.has_node_repr('AAACC').has_coverages_by_kmer([1, 1], [1, 1], [2, 1])
        expect.has_node_repr('C').has_coverages_by_kmer([1, 1])
        expect.has_node_repr('GAA').has_coverages_by_kmer([2, 1], [1, 1], [1, 1])
        expect.has_node_repr('G').has_coverages_by_kmer([0, 1])

        for color in [0, 1]:
            for edge in [['AAACC', 'C'], ['C', 'GAA']]:
                expect.has_repr_edge(edge[0], edge[1], color)
        expect.has_repr_edge('GAA', 'G', 1)
        expect.has_repr_edge('AAACC', 'GAA', 0)
        expect.has_n_edges(6)

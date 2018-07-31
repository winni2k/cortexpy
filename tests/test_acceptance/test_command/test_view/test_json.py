import os
import json

from cortexpy.test import builder as builder, runner as runner
from cortexpy.test.driver import command
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


class TestTraversal(object):
    def test_with_initial_kmers(self, tmpdir):
        # given
        record1 = 'AAACCCGAA'
        record2 = 'ACCG'
        for query_string_idx, query_string in enumerate(['CCC', record1, record2]):
            query_dir = tmpdir / str(query_string_idx)
            query_dir.ensure(dir=True)
            d = command.ViewTraversal(query_dir)
            d.with_records(record1, record2)
            d.with_kmer_size(3)
            d.with_json_output()

            # when
            expect = d.run()

            expect.has_n_nodes(3)
            expect.has_node_reprs('C', 'AAACC', 'GAA')
            expect.has_node_repr('C').has_coverages_by_kmer([1])
            expect.has_node_repr('AAACC').has_coverages_by_kmer([1], [1], [2])
            expect.has_node_repr('GAA').has_coverages_by_kmer([2], [1], [1])

            for color in [0]:
                for edge in [['AAACC', 'C'], ['C', 'GAA']]:
                    expect.has_repr_edge(edge[0], edge[1], color)
            expect.has_repr_edge('AAACC', 'GAA', 0)
            expect.has_n_edges(3)

    def test_with_revcomp_seed_kmer(self, tmpdir):
        # given
        record1 = 'AAACCCGAG'
        record2 = 'ACCG'
        query_string = 'CTC'
        d = command.ViewTraversal(tmpdir)
        d.with_records(record1, record2)
        d.with_kmer_size(3)
        d.with_seed_strings(query_string)
        d.with_json_output()

        # when
        expect = d.run()

        expect.has_node_reprs('G', 'TTT', 'CTCGG')
        expect.has_node_repr('G').has_coverages_by_kmer([1])
        expect.has_node_repr('TTT').has_coverages_by_kmer([2], [1], [1])
        expect.has_node_repr('CTCGG').has_coverages_by_kmer([1], [1], [2])

        expect.has_n_edges(3)
        expect.has_repr_edge('G', 'TTT', 0)
        expect.has_repr_edge('CTCGG', 'G', 0)
        expect.has_repr_edge('CTCGG', 'TTT', 0)

    def test_with_peripheral_edges_creates_unitigs(self, tmpdir):
        # given
        record1 = 'AAACCCGAA'
        d = command.ViewTraversal(tmpdir)
        d.with_json_output()
        d.with_kmer_size(3)
        d.with_record('ACCG')
        d.with_record(record1)
        d.with_record(record1 + 'G', name='sample_1')
        d.with_initial_contigs('CCC')
        d.with_traversal_colors(0)

        # when
        expect = d.run()

        expect.has_node_repr('C').has_coverages_by_kmer([1, 1])
        expect.has_node_repr('AAACC').has_coverages_by_kmer([1, 1], [1, 1], [2, 1])
        expect.has_node_repr('GAA').has_coverages_by_kmer([2, 1], [1, 1], [1, 1])
        expect.has_node_repr('G').has_coverages_by_kmer([1, 1])
        expect.has_n_nodes(4)

        for color in [0, 1]:
            for edge in [['AAACC', 'C'], ['C', 'GAA']]:
                expect.has_repr_edge(edge[0], edge[1], color)
        expect.has_repr_edge('AAACC', 'GAA', 0)
        expect.has_repr_edge('GAA', 'G', 1)
        expect.has_n_edges(6)

    def test_with_non_matching_start_string_returns_empty_json(self, tmpdir):
        # given
        d = command.ViewTraversal(tmpdir)
        d.with_json_output()
        d.with_records('AAACCCGAA')
        d.with_kmer_size(3)
        d.with_initial_contigs('GCGC')

        # when
        expect = d.run()

        # then
        expect.has_n_nodes(0)

    def test_of_three_colors_returns_complete_graph(self, tmpdir):
        # given
        d = command.ViewTraversal(tmpdir)
        d.with_sample_records(
            'AAAT',
            'AAAGG',
            'CAAA',
        )
        d.with_kmer_size(3)
        d.with_json_output()

        # when
        expect = d.run()

        # then
        expect.has_node_repr('A').has_coverages_by_color(1, 1, 1)
        expect.has_node_repr('CAA').has_coverages_by_color(0, 0, 1)
        expect.has_node_repr('T').has_coverages_by_color(1, 0, 0)
        expect.has_node_repr('GG').has_coverages_by_color([0, 0], [1, 1], [0, 0])
        expect.has_n_nodes(4)

        expect.has_node_unitig('AAA')
        expect.has_node_unitig('CAA')
        expect.has_node_unitig('AAT')
        expect.has_node_unitig('AAGG')

        expect.has_repr_edge('A', 'T', 0)
        expect.has_repr_edge('A', 'GG', 1)
        expect.has_repr_edge('CAA', 'A', 2)
        expect.has_n_edges(3)

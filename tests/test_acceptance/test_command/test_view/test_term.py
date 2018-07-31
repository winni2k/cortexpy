import attr

import cortexpy.test.builder as builder
from cortexpy.test import runner


@attr.s(slots=True)
class CortexpyPrintOutputParser(object):
    output = attr.ib()

    def get_kmer_strings(self):
        return self.output.rstrip().split('\n')


@attr.s(slots=True)
class ViewExpectation(object):
    output = attr.ib()
    parser = attr.ib(init=False)
    kmer_strings = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.parser = CortexpyPrintOutputParser(self.output)
        self.kmer_strings = self.parser.get_kmer_strings()

    def has_kmer(self, kmer_string):
        assert kmer_string in self.kmer_strings
        return self

    def has_n_kmers(self, n):
        assert len(self.kmer_strings) == n
        return self


class Test(object):
    def test_prints_three_kmers_including_one_revcomp(self, tmpdir):
        # given
        record = 'ACCTT'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected_kmers = [
            'AAG 1 ......G.',
            'ACC 1 .......T',
            'AGG 1 a......T',
        ]

        # when
        stdout = runner.Cortexpy().view_graph(output_graph).stdout

        # then
        assert expected_kmers == CortexpyPrintOutputParser(stdout).get_kmer_strings()


class TestTermWithRecord(object):
    def test_prints_single_kmer(self, tmpdir):
        # given
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence('ACCAA')
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected_kmer = 'CAA: CAA 1 1 .c...... ........'

        # when
        stdout = runner.Cortexpy().view_contig(graph=output_graph, contig='CAA').stdout

        # then
        assert [expected_kmer] == CortexpyPrintOutputParser(stdout).get_kmer_strings()

    def test_prints_one_missing_kmer(self, tmpdir):
        # given
        kmer_size = 3
        record = 'GGG'
        output_graph = (builder.Mccortex()
                        .with_dna_sequence('AAAA')
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected_kmer = 'CCC: GGG 0 1 ........ ........'

        # when
        stdout = runner.Cortexpy().view_contig(graph=output_graph, contig=record).stdout

        # then
        assert [expected_kmer] == CortexpyPrintOutputParser(stdout).get_kmer_strings()

    def test_prints_three_kmers(self, tmpdir):
        # given
        record = 'ACCAA'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record)
                        .with_kmer_size(kmer_size).build(tmpdir))

        expected_kmers = [
            'ACC: ACC 1 1 ....A... ....A...',
            'CCA: CCA 1 1 a...A... a...A...',
            'CAA: CAA 1 1 .c...... .c......',
        ]

        # when
        stdout = runner.Cortexpy().view_contig(graph=output_graph, contig=record).stdout

        # then
        assert expected_kmers == CortexpyPrintOutputParser(stdout).get_kmer_strings()

    def test_prints_three_kmers_including_one_revcomp(self, tmpdir):
        # given
        record = 'ACCTT'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected_kmers = [
            'ACC: ACC 1 1 .......T .......T',
            'AGG: CCT 1 1 A......t A......t',
            'AAG: CTT 1 1 .C...... .C......',
        ]

        # when
        stdout = runner.Cortexpy().view_contig(graph=output_graph, contig=record).stdout

        # then
        assert expected_kmers == CortexpyPrintOutputParser(stdout).get_kmer_strings()

    def test_prints_one_missing_and_one_revcomp_kmer(self, tmpdir):
        # given
        dna_sequence = 'ACCTT'
        search_record = 'ACTT'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(dna_sequence)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected_kmers = [
            'ACT: ACT 0 1 ........ .......T',
            'AAG: CTT 1 1 .C...... A.......',
        ]

        # when
        stdout = runner.Cortexpy().view_contig(graph=output_graph, contig=search_record).stdout
        expect = ViewExpectation(stdout)

        # then
        (expect
         .has_kmer(expected_kmers[0])
         .has_kmer(expected_kmers[1])
         .has_n_kmers(2))

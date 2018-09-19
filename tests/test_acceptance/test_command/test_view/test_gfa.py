import attr
import pytest

from cortexpy.test import runner
from cortexpy.test.expectation import Fasta


@attr.s(slots=True)
class GfaBuilder:
    headers = attr.ib(attr.Factory(list))
    segments = attr.ib(attr.Factory(list))
    links = attr.ib(attr.Factory(list))

    def build_lines(self):
        self.headers.insert(0, 'H\tVN:Z:1.0')
        lines = []
        lines += self.headers
        lines += self.segments
        lines += self.links
        return lines

    def build(self):
        return '\n'.join(self.build_lines())


class TestTraversal:
    @pytest.mark.skip
    def test_with_empty_input_outputs_empty_fasta(self, tmpdir):
        # given
        input_graph = str(tmpdir / 'input.gfa')
        output_graph = tmpdir / 'output.fast'
        with open(input_graph, 'w') as fh:
            fh.write(GfaBuilder().build())

        # when
        runner.Cortexpy().traverse(input_graph, out=output_graph, input_gfa=True)

        # then
        expect = Fasta(output_graph.read())
        expect.has_no_records()

    def test_with_single_unitig_outputs_single_record(self):
        pass

    def test_with_x_shaped_unitig_graph_outputs_four_records(self):
        pass

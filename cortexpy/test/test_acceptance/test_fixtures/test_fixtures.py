import os
import cortexpy.graph.parser as parser


class TestManyColorsHeaderOnlyFixture(object):
    def test_does_not_raise(self):
        # given
        graph_header_fixture = os.path.join(os.path.dirname(__file__),
                                            'many_colors_header_only.ctx')

        # when
        with open(graph_header_fixture, 'rb') as header_handle:
            header = parser.header.from_stream(header_handle)

        # then
        assert header.version == 6
        assert header.kmer_size == 47
        assert header.kmer_container_size == 2
        assert header.num_colors == 25

import contextlib
import io
import os
import subprocess
from pathlib import Path

import attr
import sys
import logging
from cortexpy.__main__ import main

logger = logging.getLogger(__name__)


@attr.s(slots=True)
class Mccortex(object):
    kmer_size = attr.ib()
    mccortex_bin = attr.ib(None)

    def __attrs_post_init__(self):
        if self.mccortex_bin is None:
            bin_dir = os.environ['BIN_DIR']
            self.mccortex_bin = os.path.join(bin_dir, 'mccortex')

    def view(self, graph):
        return self.run(['view', '-k', graph])

    def clean(self, graphs, out, *, clip_tips_shorter_than=0,
              unitigs_with_mean_coverage_less_than=1):
        if isinstance(graphs, str):
            graphs = [graphs]
        return self.run(['clean', '--tips={}'.format(str(clip_tips_shorter_than)),
                         '--unitigs={}'.format(str(unitigs_with_mean_coverage_less_than)), '--out',
                         out, *graphs])

    def run(self, mccortex_args):
        if self.kmer_size <= 31:
            mccortex_k = 31
        elif self.kmer_size <= 63:
            mccortex_k = 63
        else:
            raise Exception(
                'Kmer size ({}) too big for available mccortex binaries'.format(self.kmer_size))
        command = [self.mccortex_bin, str(mccortex_k)]
        command += mccortex_args
        command = [str(arg) for arg in command]

        return subprocess.run(command, check=True)


@attr.s(slots=True)
class Cortexpy(object):
    spawn_process = attr.ib(False)

    def view_graph(self, graph):
        return self.run(['view', 'graph', graph])

    def view_contig(self, contig, graph, to_json=False, other_args=()):
        run_args = []
        if to_json:
            run_args.append('--to-json')
        run_args += list(other_args)
        return self.run(['view', 'contig', str(graph), contig] + run_args)

    def view_traversal(self, contig, graph,
                       to_json=None, kmers=None,
                       color=0, max_nodes=None, colors=None,
                       # deprecated options
                       output_format=None, output_type=None,):
        assert output_format is None or to_json is None
        assert output_type is None or kmers is None
        if (output_format is None and to_json is None) or output_format == 'json':
            to_json = True
        if (output_type is None and kmers is None) or output_format == 'kmers':
            kmers = True

        if colors is None:
            colors = [color]
        if isinstance(colors, int):
            colors = [colors]
        intermediate_graph = str(Path(graph).with_suffix('.traverse.pickle'))
        command1 = [
            'traverse', str(graph),
            '--initial-contig',
            contig, '--colors',
            ','.join(str(c) for c in colors),
            '--out', intermediate_graph
        ]
        if max_nodes is not None:
            command1 += ['--max-nodes', str(max_nodes)]
        command1_ret = self.run(command1)

        command2 = ['view', 'traversal', intermediate_graph, ]
        if to_json:
            command2.append('--to-json')
        if kmers:
            command2.append('--kmers')
        command2_ret = self.run(command2)

        stdout = command1_ret.stdout + command2_ret.stdout
        stderr = command1_ret.stderr + command2_ret.stderr
        returncode = command2_ret.returncode

        return subprocess.CompletedProcess(command1 + command2,
                                           returncode,
                                           stdout=stdout,
                                           stderr=stderr)

    def assemble(self, *, graph, initial_seqs):
        command = ['assemble', graph, initial_seqs]
        return self.run([str(c) for c in command])

    def run(self, args):
        logger.info('Running with args: {}'.format(args))
        if self.spawn_process:
            command = [sys.executable, '-m', 'cortexpy'] + args
            return subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  universal_newlines=True)
        else:
            stdout = io.StringIO()
            stderr = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                with contextlib.redirect_stderr(stderr):
                    main(args)
            return subprocess.CompletedProcess(args,
                                               0,
                                               stdout=stdout.getvalue(),
                                               stderr=stderr.getvalue())

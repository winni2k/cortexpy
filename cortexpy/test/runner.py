import contextlib
import io
import os
import subprocess

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

        return subprocess.check_call(command)


@attr.s(slots=True)
class Cortexpy(object):
    spawn_process = attr.ib(False)

    def view_graph(self, graph):
        return self.run(['view', 'graph', graph])

    def view_contig(self, contig, graph, output_format=None, other_args=()):
        run_args = []
        if output_format is not None:
            run_args.extend(['--output-format', output_format])
        run_args += list(other_args)
        return self.run(['view', 'contig', graph, contig] + run_args)

    def view_traversal(self, contig, graph, output_format='json', output_type='kmers',
                       orientation='both', color=0, max_nodes=None, colors=None):
        if colors is None:
            colors = [color]
        command = ['view', 'traversal',
                   graph, contig,
                   '--output-format', output_format,
                   '--colors', ','.join(str(c) for c in colors),
                   '--output-type', output_type]
        if max_nodes is not None:
            command.extend(['--max-nodes', str(max_nodes)])
        return self.run(command)

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

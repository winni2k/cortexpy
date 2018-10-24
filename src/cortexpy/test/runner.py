import contextlib
import io
import logging
import os
import subprocess
import sys

import attr

from cortexpy.__main__ import main

logger = logging.getLogger(__name__)


@attr.s(slots=True)
class Mccortex(object):
    kmer_size = attr.ib()
    mccortex_bin = attr.ib('mccortex')

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

        return subprocess.run(command, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


@attr.s(slots=True)
class Cortexpy(object):
    spawn_process = attr.ib(False)

    def view_graph(self, graph, kmers=False, out=None):
        args = ['view', 'graph', graph]
        if kmers:
            args.append('--kmers')
        if out is not None:
            args += ['--out', out]
        return self.run(args)

    def view_contig(self, contig, graph, to_json=False, other_args=()):
        run_args = []
        if to_json:
            run_args.append('--to-json')
        run_args += list(other_args)
        return self.run(['view', 'contig', str(graph), contig] + run_args)

    def traverse(self, graph,
                 contig=None,
                 to_json=False,
                 out=None,
                 max_paths=None,
                 input_gfa=False,
                 kmers=False,
                 graph_index=None,
                 extra_start_kmer=None,
                 links_file=None):
        cmd = ['traverse', graph]
        if to_json:
            cmd.append('--to-json')
        if contig is not None:
            cmd += ['--seed-strings', contig]
        if max_paths:
            cmd += ['--max-paths', max_paths]
        if input_gfa:
            cmd.append('--from-gfa')
        if kmers:
            cmd.append('--kmers')
        if graph_index is not None:
            cmd += ['--graph-index', graph_index]
        if extra_start_kmer is not None:
            cmd += ['--extra-start-kmer', extra_start_kmer]
        if out:
            cmd += ['--out', out]
        if links_file:
            cmd += ['--links-file', links_file]
        cmd += ['-v']
        return self.run(cmd)

    def subgraph(self, *, graphs, contig, out='/dev/null',
                 contig_fasta=False, colors=None,
                 max_nodes=None, verbose=False,
                 silent=False, logging_interval=None):
        cmd = ['subgraph', contig, '--out', out]
        assert len(graphs) > 0
        cmd.append('--graphs')
        cmd += graphs
        if colors:
            cmd.append('--colors')
            cmd += colors
        if contig_fasta:
            cmd.append('--initial-fasta')
        if max_nodes:
            cmd.extend(['--max-nodes', max_nodes])
        if verbose:
            cmd.append('--verbose')
        if silent:
            cmd.append('--silent')
        if logging_interval is not None:
            cmd += ['--logging-interval', logging_interval]
        return self.run(cmd)

    def prune(self, *, graph, out, remove_tips=None, verbose=None):
        cmd = ['prune', graph, '--out', out]
        if remove_tips:
            cmd.extend(['--remove-tips', remove_tips])
        if verbose:
            cmd.append('--verbose')
        cmd = [str(c) for c in cmd]
        return self.run(cmd)

    def assemble(self, *, graph, initial_seqs, out='/dev/null'):
        command = ['assemble', graph, initial_seqs, '--out', out]
        return self.run(command)

    def run(self, args):
        args = [str(a) for a in args]
        logger.info('Running with args: "{}"'.format(' '.join(args)))
        if self.spawn_process:
            command = [sys.executable, '-m', 'cortexpy'] + args
            return subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  universal_newlines=True)
        else:
            args = ['cortexpy'] + args
            stdout = io.StringIO()
            stderr = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                with contextlib.redirect_stderr(stderr):
                    main(args)
            return subprocess.CompletedProcess(args,
                                               0,
                                               stdout=stdout.getvalue(),
                                               stderr=stderr.getvalue())

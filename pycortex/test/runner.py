import contextlib
import io
import os
import subprocess

import attr

from pycortex.__main__ import main


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


def run_in_process(args):
    pycortex_stdout = io.StringIO()
    pycortex_stderr = io.StringIO()
    exit_status = 0
    with contextlib.redirect_stdout(pycortex_stdout):
        with contextlib.redirect_stdout(pycortex_stderr):
            try:
                main(args)
            except Exception as e:
                print(e, file=pycortex_stderr)
                exit_status = 1
    return pycortex_stdout.getvalue(), pycortex_stderr.getvalue(), exit_status


@attr.s(slots=True)
class Pycortex(object):
    spawn_process = attr.ib(False)

    def view(self, args):
        return self.run(['view'] + args)

    def run(self, args):
        if self.spawn_process:
            raise NotImplementedError
        else:
            return run_in_process(args)

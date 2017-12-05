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
class cortexpy(object):
    spawn_process = attr.ib(False)

    def view(self, args):
        return self.run(['view'] + args)

    def run(self, args):
        if self.spawn_process:
            command = [sys.executable, '-m', 'cortexpy'] + args
            print('Running: ' + ' '.join(command))
            return subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            main(args)
            return subprocess.CompletedProcess(args, 0)

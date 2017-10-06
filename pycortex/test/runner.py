import os
import subprocess

import attr

BIN_DIR = os.environ['BIN_DIR']
MCCORTEX = os.path.join(BIN_DIR, 'mccortex')


@attr.s(slots=True)
class Mccortex(object):
    kmer_size = attr.ib()

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
        command = [MCCORTEX, str(mccortex_k)]
        command += mccortex_args

        return subprocess.check_call(command)


@attr.s(slots=True)
class Pycortex(object):
    spawn_process = attr.ib(False)

    def run(self):
        pass

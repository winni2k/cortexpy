from cortexpy.test import driver
import cortexpy.test.driver.command


class Test(object):
    def test_creates_two_transcripts_from_four_records(self, tmpdir):
        # given
        d = driver.command.Assemble(tmpdir)
        d.with_records(
            'CCCC',
            'CCCGG',
            'AAAA',
            'AAATAA',
        )
        d.with_initial_sequences('AAA', 'CCCC')
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_record('CCCCGGG')
        expect.has_record('AAAATAA')
        expect.has_n_records(2)

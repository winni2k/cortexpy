import cortexpy.test.driver.command as command


class Test(object):
    def test_creates_two_transcripts_from_four_records(self, tmpdir):
        # given
        d = command.Assemble(tmpdir)
        d.with_records(
            'CCCG',
            'CCGC',
            'AAAT',
            'AATA',
        )
        d.with_initial_sequences('AAA', 'CCC')
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_record('AAATA')
        expect.has_record('CCCGC')
        expect.has_n_records(2)

    def test_creates_two_transcripts_from_four_records_in_four_colors(self, tmpdir):
        # given
        d = command.Assemble(tmpdir)
        for seq in ['CCCG', 'CCGC', 'AAAT', 'AATA']:
            d.with_dna_sequence(seq, name=seq)
        d.with_initial_sequences('AAA', 'CCC')
        d.with_kmer_size(3)

        # when
        expect = d.run()

        # then
        expect.has_record('AAATA')
        expect.has_record('CCCGC')
        expect.has_n_records(2)

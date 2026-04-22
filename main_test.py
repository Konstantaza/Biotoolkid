from main import DNASequence, NucleicAcidSequence, filter_fastq, AminoAcidSequence, RNASequence
import pytest


def test_dna_complement():
    seq = 'ATGCTAGCTGAC'
    test = DNASequence(seq)
    assert test.complement().sequence == 'TACGATCGACTG'


def test_class_error():
    seq = 'ATGCTAGCTGAC'
    with pytest.raises(TypeError):
        test = NucleicAcidSequence(seq)


def test_filter_fastq_file_creation(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    input_file = tmp_path / "test_input.fastq"
    input_file.write_text("@read1\nATGC\n+\nIIII\n")
    filter_fastq("test_input.fastq", "output.fastq")
    expected_file = tmp_path / "filtered" / "output.fastq"
    assert expected_file.exists()


def test_count_hydrophobic():
    seq = 'MAVW'
    test = AminoAcidSequence(seq)
    assert test.count_hydrophobic() == 4


def test_transcribe():
    seq = 'ATGCTAGCTGAC'
    dna = DNASequence(seq)
    test = dna.transcribe()
    assert isinstance(test, RNASequence)
    assert test.sequence == 'AUGCUAGCUGAC'


def test_check_alphabet_is_true():
    seq = 'ATGCTAGCTGAC'
    dna = DNASequence(seq)
    rna = dna.transcribe()
    aa = AminoAcidSequence('MAVW')
    assert dna.check_alphabet() == True
    assert rna.check_alphabet() == True
    assert aa.check_alphabet() == True
    

def test_check_alphabet_is_false():
    seq = 'ATGCTAGCTGACZ'
    dna = DNASequence(seq)
    rna = RNASequence('AUGCUAGCUGACZ')
    aa = AminoAcidSequence('MAVWZ')
    assert dna.check_alphabet() == False
    assert rna.check_alphabet() == False
    assert aa.check_alphabet() == False


def test_seq_length():
    seq = 'ATGCTAGCTGAC'
    dna = DNASequence(seq)
    assert len(dna) == 12
    




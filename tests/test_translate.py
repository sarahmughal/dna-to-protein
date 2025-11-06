from src.translate import dna_to_rna, translate_dna

def test_dna_to_rna_basic():
    assert dna_to_rna('ACTT') == 'ACUU'

def test_translate_simple():
    # ATG -> AUG -> M
    assert translate_dna('ATG') == 'M'

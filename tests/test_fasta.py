import io, os
from src.io_utils import read_fasta, write_fasta

def test_fasta_roundtrip(tmp_path):
    recs = {'a':'ACGT', 'b':'GGGTTTAAA'}
    p = tmp_path/'x.fasta'
    write_fasta(str(p), recs, width=4)
    out = read_fasta(str(p))
    assert out == recs

from src.align import needleman_wunsch, smith_waterman

def test_global_alignment():
    a,b,s = needleman_wunsch('GATTACA','GCATGCU')
    assert isinstance(s, int)
    assert len(a)==len(b)

def test_local_alignment():
    a,b,s = smith_waterman('GATTACA','GCATGCU')
    assert s >= 0

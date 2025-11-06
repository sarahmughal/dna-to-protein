from src.tree import upgma, to_newick
from src.distance import distance_matrix

def test_upgma_small():
    seqs = ['AAAA','AAAT','AATT']
    names = ['s1','s2','s3']
    D = distance_matrix(seqs, model='p')
    root = upgma(names, D)
    nwk = to_newick(root) + ';'
    assert nwk.endswith(';')

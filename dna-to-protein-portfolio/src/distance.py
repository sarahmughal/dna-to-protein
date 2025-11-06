import math
from typing import List
from .align import needleman_wunsch

def p_distance(a: str, b: str) -> float:
    """Proportion of differing sites after global alignment."""
    al_a, al_b, _ = needleman_wunsch(a, b)
    matches = sum(1 for x,y in zip(al_a, al_b) if x==y and x!='-' and y!='-')
    comps = sum(1 for x,y in zip(al_a, al_b) if x!='-' and y!='-')
    if comps == 0:
        return 1.0
    return 1 - matches/comps

def jukes_cantor(p: float) -> float:
    """Estimate substitutions/site under JC69 from p-distance."""
    if p >= 0.75:
        return float('inf')
    return -3/4 * math.log(1 - 4*p/3)

def distance_matrix(seqs: List[str], model: str = 'p') -> List[List[float]]:
    n = len(seqs)
    D = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            p = p_distance(seqs[i], seqs[j])
            d = p if model=='p' else jukes_cantor(p)
            D[i][j] = D[j][i] = d
    return D

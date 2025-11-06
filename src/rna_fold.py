# Nussinov algorithm for RNA secondary structure (max base pairs)
from typing import List, Tuple

def can_pair(a: str, b: str) -> bool:
    pairs = {('A','U'),('U','A'),('G','C'),('C','G'),('G','U'),('U','G')}
    return (a,b) in pairs

def nussinov(rna: str, min_loop: int = 0) -> str:
    n = len(rna)
    dp = [[0]*n for _ in range(n)]
    bt = [[None]*n for _ in range(n)]
    for k in range(1, n):
        for i in range(n-k):
            j = i + k
            best = dp[i+1][j]
            choice = ('skip_i', i+1, j)
            if dp[i][j-1] > best:
                best = dp[i][j-1]; choice = ('skip_j', i, j-1)
            if can_pair(rna[i], rna[j]) and (j - i - 1) >= min_loop:
                score = dp[i+1][j-1] + 1
                if score > best:
                    best = score; choice = ('pair', i+1, j-1)
            # bifurcation
            for t in range(i+1, j):
                score = dp[i][t] + dp[t+1][j]
                if score > best:
                    best = score; choice = ('split', i, t, t+1, j)
            dp[i][j] = best
            bt[i][j] = choice
    # traceback to dot-bracket
    res = ['.']*n
    def tb(i,j):
        if i>=j or bt[i][j] is None:
            return
        c = bt[i][j]
        if c[0]=='skip_i':
            tb(c[1], c[2])
        elif c[0]=='skip_j':
            tb(c[1], c[2])
        elif c[0]=='pair':
            res[i] = '('
            res[j] = ')'
            tb(c[1], c[2])
        else: # split
            _, i1, t, t1, j1 = c
            tb(i1, t); tb(t1, j1)
    tb(0, n-1)
    return ''.join(res)

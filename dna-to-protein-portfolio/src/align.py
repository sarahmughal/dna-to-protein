from typing import Tuple, List
from dataclasses import dataclass

@dataclass
class Scoring:
    match: int = 1
    mismatch: int = -1
    gap: int = -1

def needleman_wunsch(a: str, b: str, scoring: Scoring = Scoring()) -> Tuple[str,str,int]:
    """Global alignment."""
    m, n = len(a), len(b)
    # DP matrices
    score = [[0]*(n+1) for _ in range(m+1)]
    ptr = [[None]*(n+1) for _ in range(m+1)]
    for i in range(1, m+1):
        score[i][0] = i*scoring.gap
        ptr[i][0] = 'U'
    for j in range(1, n+1):
        score[0][j] = j*scoring.gap
        ptr[0][j] = 'L'
    for i in range(1, m+1):
        for j in range(1, n+1):
            diag = score[i-1][j-1] + (scoring.match if a[i-1]==b[j-1] else scoring.mismatch)
            up = score[i-1][j] + scoring.gap
            left = score[i][j-1] + scoring.gap
            best = max(diag, up, left)
            score[i][j] = best
            if best == diag:
                ptr[i][j] = 'D'
            elif best == up:
                ptr[i][j] = 'U'
            else:
                ptr[i][j] = 'L'
    # traceback
    i, j = m, n
    al_a, al_b = [], []
    while i>0 or j>0:
        p = ptr[i][j]
        if p == 'D':
            al_a.append(a[i-1]); al_b.append(b[j-1])
            i -= 1; j -= 1
        elif p == 'U':
            al_a.append(a[i-1]); al_b.append('-')
            i -= 1
        else: # 'L'
            al_a.append('-'); al_b.append(b[j-1])
            j -= 1
    return ''.join(reversed(al_a)), ''.join(reversed(al_b)), score[m][n]

def smith_waterman(a: str, b: str, scoring: Scoring = Scoring()) -> Tuple[str,str,int]:
    """Local alignment."""
    m, n = len(a), len(b)
    score = [[0]*(n+1) for _ in range(m+1)]
    ptr = [[None]*(n+1) for _ in range(m+1)]
    best_i, best_j, best_score = 0, 0, 0
    for i in range(1, m+1):
        for j in range(1, n+1):
            diag = score[i-1][j-1] + (scoring.match if a[i-1]==b[j-1] else scoring.mismatch)
            up = score[i-1][j] + scoring.gap
            left = score[i][j-1] + scoring.gap
            val = max(0, diag, up, left)
            score[i][j] = val
            if val == 0:
                ptr[i][j] = None
            elif val == diag:
                ptr[i][j] = 'D'
            elif val == up:
                ptr[i][j] = 'U'
            else:
                ptr[i][j] = 'L'
            if val > best_score:
                best_score, best_i, best_j = val, i, j
    # traceback from best
    i, j = best_i, best_j
    al_a, al_b = [], []
    while i>0 and j>0 and score[i][j]>0:
        p = ptr[i][j]
        if p == 'D':
            al_a.append(a[i-1]); al_b.append(b[j-1])
            i -= 1; j -= 1
        elif p == 'U':
            al_a.append(a[i-1]); al_b.append('-')
            i -= 1
        elif p == 'L':
            al_a.append('-'); al_b.append(b[j-1])
            j -= 1
        else:
            break
    return ''.join(reversed(al_a)), ''.join(reversed(al_b)), best_score

def center_star_msa(seqs: List[str]) -> List[str]:
    """Very simple MSA: choose center that maximizes sum of pairwise global alignment scores,
    then align others to the center and merge columns.
    """
    if len(seqs) == 1:
        return [seqs[0]]
    # compute scores
    scores = [[0]*len(seqs) for _ in range(len(seqs))]
    for i in range(len(seqs)):
        for j in range(i+1, len(seqs)):
            _,_,s = needleman_wunsch(seqs[i], seqs[j])
            scores[i][j] = scores[j][i] = s
    center = max(range(len(seqs)), key=lambda i: sum(scores[i]))
    # build MSA by aligning to center
    aligned = [None]*len(seqs)
    aligned[center] = seqs[center]
    for i in range(len(seqs)):
        if i == center:
            continue
        ac, si, _ = needleman_wunsch(aligned[center], seqs[i])
        # merge gaps across all already-aligned sequences
        # expand existing aligned sequences wherever 'ac' has gaps
        def expand(s, template):
            k = 0
            out = []
            for ch in template:
                if ch == '-':
                    out.append('-')
                else:
                    out.append(s[k]); k += 1
            return ''.join(out)
        # expand current center across template ac (which includes gaps)
        if aligned[center] == seqs[center]:
            aligned[center] = ac
        else:
            aligned[center] = expand(aligned[center].replace('-', ''), ac)
        # now add the new sequence aligning with 'ac' vs 'si'
        aligned[i] = si
        # for every previously added sequence, expand to match the current center length
        for j in range(len(seqs)):
            if aligned[j] is not None and len(aligned[j]) < len(aligned[center]):
                aligned[j] = expand(aligned[j].replace('-', ''), aligned[center])
    return aligned

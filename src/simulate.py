import random
from typing import Dict
random.seed(42)

NUCS = ['A','C','G','T']

def random_dna(length: int, gc: float = 0.5) -> str:
    p_gc = gc/2
    p_at = (1-gc)/2
    choices = ['A','C','G','T']
    weights = [p_at, p_gc, p_gc, p_at]
    return ''.join(random.choices(choices, weights=weights, k=length))

def mutate(seq: str, mu: float) -> str:
    out = []
    for ch in seq:
        if random.random() < mu:
            out.append(random.choice([n for n in NUCS if n!=ch]))
        else:
            out.append(ch)
    return ''.join(out)

def simulate_family(n: int = 5, length: int = 120, gc: float = 0.5, mu: float = 0.05) -> Dict[str,str]:
    root = random_dna(length, gc=gc)
    fam = {"root": root}
    for i in range(1, n+1):
        fam[f"seq{i}"] = mutate(root, mu)
    return fam

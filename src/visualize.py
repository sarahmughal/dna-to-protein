
import math
from typing import Dict, List, Tuple
import matplotlib.pyplot as plt

def _savefig(path: str, title: str = None):
    if title:
        plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()

def plot_gc_content(records: Dict[str, str], out_path: str):
    names = list(records.keys())
    vals = []
    for s in records.values():
        s = s.upper()
        vals.append((s.count('G') + s.count('C')) / len(s) * 100.0 if s else 0.0)
    plt.figure(figsize=(8,4))
    plt.bar(range(len(names)), vals)
    plt.xticks(range(len(names)), names, rotation=45, ha='right')
    plt.ylabel('GC%')
    _savefig(out_path, title='GC Content per Sequence')

def codon_usage(seq: str):
    s = seq.upper().replace('T','U')
    counts = {}
    for i in range(0, len(s)-2, 3):
        c = s[i:i+3]
        if len(c)==3 and all(ch in 'AUCG' for ch in c):
            counts[c] = counts.get(c,0)+1
    return counts

def plot_codon_usage(seq: str, out_path: str):
    cu = codon_usage(seq)
    items = sorted(cu.items())
    plt.figure(figsize=(10,4))
    plt.bar(range(len(items)), [v for _,v in items])
    plt.xticks(range(len(items)), [k for k,_ in items], rotation=90)
    plt.ylabel('Count')
    _savefig(out_path, title='Codon Usage (frame 0)')

def plot_distance_heatmap(D: List[List[float]], names: List[str], out_path: str, title: str='Distance Heatmap'):
    plt.figure(figsize=(5,4))
    plt.imshow(D, aspect='auto')
    plt.colorbar(label='Distance')
    plt.xticks(range(len(names)), names, rotation=45, ha='right')
    plt.yticks(range(len(names)), names)
    _savefig(out_path, title=title)

def plot_alignment(aligned: List[str], names: List[str], out_path: str):
    if not aligned: 
        return
    cols = len(aligned[0])
    consensus = []
    for j in range(cols):
        col = [s[j] for s in aligned]
        bases = [b for b in col if b!='-']
        consensus.append(max(set(bases), key=bases.count) if bases else '-')
    mat = []
    for s in aligned:
        row = [1 if ch==consensus[j] and ch!='-' else 0 for j,ch in enumerate(s)]
        mat.append(row)
    plt.figure(figsize=(max(6, cols*0.15), max(3, len(aligned)*0.3)))
    plt.imshow(mat, aspect='auto')
    plt.colorbar(label='Match to consensus')
    plt.yticks(range(len(names)), names)
    plt.xlabel('Alignment column')
    _savefig(out_path, title='Multiple Sequence Alignment (consensus match)')

def plot_rna_arcs(rna: str, dot_bracket: str, out_path: str):
    pairs = []
    stack = []
    for i,ch in enumerate(dot_bracket):
        if ch=='(':
            stack.append(i)
        elif ch==')' and stack:
            j = stack.pop()
            pairs.append((j,i))
    x = list(range(len(rna)))
    y = [0]*len(rna)
    plt.figure(figsize=(max(8, len(rna)*0.2), 3))
    plt.plot(x, y)
    for i,j in pairs:
        cx = (i+j)/2
        r = (j-i)/2
        pts = 60
        xs = [cx + r*math.cos(math.pi*t/pts) for t in range(pts+1)]
        ys = [r*math.sin(math.pi*t/pts) for t in range(pts+1)]
        plt.plot(xs, ys)
    plt.xticks(range(len(rna)), list(rna), rotation=90)
    plt.yticks([])
    _savefig(out_path, title='RNA Base-Pair Arcs')

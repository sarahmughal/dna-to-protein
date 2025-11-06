
from typing import Dict, List, Tuple, Optional
import math
import matplotlib.pyplot as plt

# ---- Theme helpers ----------------------------------------------------------
_THEME = {
    "gc_bar": {"color": "#6C8AE4"},            # periwinkle
    "codon_bar": {"color": "#F7A072"},         # coral
    "heatmap_cmap": "magma",                   # pretty dark-to-bright
    "msa_cmap": "viridis",                     # perceptually uniform
    "arc_color": "#8BD3DD",                    # teal
    "accent": "#7DCE82",                       # green accent (unused placeholder)
    "font_size": 11
}

def set_theme(font_size: int = 11):
    """Apply a clean Matplotlib style and font sizes."""
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": "#222",
        "axes.labelcolor": "#222",
        "xtick.color": "#222",
        "ytick.color": "#222",
        "font.size": font_size,
        "axes.titlesize": font_size + 1,
        "axes.labelsize": font_size,
        "xtick.labelsize": font_size - 1,
        "ytick.labelsize": font_size - 1,
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
    })

def _savefig(path: str, title: Optional[str] = None, xlabel: Optional[str] = None, ylabel: Optional[str] = None):
    if title:
        plt.title(title, pad=10)
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()

# ---- Plots ------------------------------------------------------------------
def plot_gc_content(records: Dict[str, str], out_path: str):
    """Bar chart of GC% per sequence (how G/C-heavy each is)."""
    set_theme(_THEME["font_size"])
    names = list(records.keys())
    vals = []
    for s in records.values():
        s = s.upper()
        vals.append((s.count('G') + s.count('C')) / len(s) * 100.0 if s else 0.0)
    plt.figure(figsize=(9, 4))
    plt.bar(range(len(names)), vals, **_THEME["gc_bar"])
    plt.xticks(range(len(names)), names, rotation=35, ha='right')
    _savefig(out_path, title="GC Content per Sequence", ylabel="GC (%)")

def codon_usage(seq: str):
    """Return dict {codon: count} for frame-0 RNA (T→U)."""
    s = seq.upper().replace('T','U')
    counts = {}
    for i in range(0, len(s)-2, 3):
        c = s[i:i+3]
        if len(c)==3 and all(ch in 'AUCG' for ch in c):
            counts[c] = counts.get(c,0)+1
    return counts

def plot_codon_usage(seq: str, out_path: str):
    """Bar chart of codon counts (frame 0)."""
    set_theme(_THEME["font_size"])
    cu = codon_usage(seq)
    items = sorted(cu.items())
    xs = list(range(len(items)))
    counts = [v for _, v in items]
    labels = [k for k, _ in items]
    plt.figure(figsize=(12, 4.8))
    plt.bar(xs, counts, **_THEME["codon_bar"])
    plt.xticks(xs, labels, rotation=90)
    _savefig(out_path, title="Codon Usage (frame 0)", ylabel="Count")

def plot_distance_heatmap(D: List[List[float]], names: List[str], out_path: str, title: str='Distance Heatmap'):
    """Heatmap of pairwise distances (p or JC69)."""
    set_theme(_THEME["font_size"])
    plt.figure(figsize=(6.5, 5.2))
    im = plt.imshow(D, aspect='auto', cmap=_THEME["heatmap_cmap"])
    plt.colorbar(im, label="Distance")
    plt.xticks(range(len(names)), names, rotation=45, ha='right')
    plt.yticks(range(len(names)), names)
    _savefig(out_path, title=title)

def plot_alignment(aligned: List[str], names: List[str], out_path: str):
    """Consensus-match visualization of a multiple alignment (bright = match)."""
    if not aligned:
        return
    set_theme(_THEME["font_size"])
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
    plt.figure(figsize=(max(7, cols*0.2), max(3.2, len(aligned)*0.38)))
    im = plt.imshow(mat, aspect='auto', cmap=_THEME["msa_cmap"], vmin=0, vmax=1)
    plt.colorbar(im, label="Match to consensus")
    plt.yticks(range(len(names)), names)
    plt.xlabel('Alignment column')
    _savefig(out_path, title="Multiple Sequence Alignment (Consensus Match)")

def plot_rna_arcs(rna: str, dot_bracket: str, out_path: str):
    """Arc diagram from dot-bracket (paired bases draw semicircles)."""
    set_theme(_THEME["font_size"])
    pairs = []
    stack = []
    for i, ch in enumerate(dot_bracket):
        if ch == '(':
            stack.append(i)
        elif ch == ')' and stack:
            j = stack.pop()
            pairs.append((j, i))
    x = list(range(len(rna)))
    y = [0] * len(rna)
    plt.figure(figsize=(max(8, len(rna)*0.22), 3.2))
    # baseline
    plt.plot(x, y, color="#999", linewidth=1)
    # arcs
    for i, j in pairs:
        cx = (i + j) / 2
        r = (j - i) / 2
        pts = 64
        xs = [cx + r * math.cos(math.pi * t / pts) for t in range(pts + 1)]
        ys = [r * math.sin(math.pi * t / pts) for t in range(pts + 1)]
        plt.plot(xs, ys, color=_THEME["arc_color"], linewidth=1.6)
    plt.xticks(range(len(rna)), list(rna), rotation=90)
    plt.yticks([])
    _savefig(out_path, title="RNA Base-Pair Arcs")

import streamlit as st
from pathlib import Path

from src.io_utils import read_fasta
from src.translate import dna_to_rna, translate_dna
from src.align import needleman_wunsch, smith_waterman, center_star_msa
from src.distance import distance_matrix
from src.tree import upgma, to_newick
from src.rna_fold import nussinov
from src.visualize_pretty import (
    plot_gc_content,
    plot_codon_usage,
    plot_distance_heatmap,
    plot_alignment,
    plot_rna_arcs,
)

# ⭐ If you enabled the pretty theme earlier, keep this import:
from src.visualize_pretty import set_theme as pretty_set_theme  # optional

st.set_page_config(page_title="Genetic Visualization App", layout="wide")
st.title("🧬 Genetic Visualization App")

with st.sidebar:
    uploaded = st.file_uploader("Upload FASTA", type=["fasta", "fa", "fna"])
    use_demo = st.checkbox("Use demo dataset", value=not bool(uploaded))
    model = st.radio("Distance model", ["p", "jc"], index=1, horizontal=True)
    mode = st.radio("Pairwise", ["global", "local"], index=0, horizontal=True)
    use_pretty = st.checkbox("Use pretty theme", value=True)
    # ⭐ NEW: a run button so we don’t recompute on every small change
    run = st.button("Run analysis")

work_dir = Path("streamlit_outputs")
work_dir.mkdir(exist_ok=True)

if use_pretty:
    pretty_set_theme(font_size=12)


# ---------------------------
# ⭐ NEW: Cached heavy helpers
# ---------------------------
@st.cache_data(show_spinner=False)
def cached_msa(seqs_tuple):
    # seqs_tuple: tuple of strings (hashable for cache)
    return center_star_msa(list(seqs_tuple))


@st.cache_data(show_spinner=False)
def cached_distance(seqs_tuple, model):
    return distance_matrix(list(seqs_tuple), model=model)


@st.cache_data(show_spinner=False)
def cached_tree(names_tuple, D):
    return to_newick(upgma(list(names_tuple), D)) + ";"


# ---------------------------
# Load records (same as before)
# ---------------------------
def load_records():
    if uploaded is not None and not use_demo:
        data = uploaded.read().decode("utf-8")
        path = work_dir / "uploaded.fasta"
        path.write_text(data, encoding="utf-8")
        return read_fasta(str(path)), str(path)
    else:
        demo_path = Path("data") / "example_group.fasta"
        if not demo_path.exists():
            demo_path = work_dir / "demo.fasta"
            demo_path.write_text(
                ">seqA\nGATTACA\n>seqB\nGCATGCU\n>seqC\nGACTATA\n", encoding="utf-8"
            )
        return read_fasta(str(demo_path)), str(demo_path)


recs, fasta_path = load_records()
names = list(recs.keys())
seqs = list(recs.values())

st.write(f"Loaded **{len(seqs)}** sequences from `{Path(fasta_path).name}`")
st.code(
    "\n".join(
        [f">{h}\n{recs[h][:60]}{'...' if len(recs[h]) > 60 else ''}" for h in names]
    ),
    language="text",
)

# ⭐ NEW: size guardrails (tune as you like)
MAX_SEQS_FOR_MSA = 200  # cap number of sequences
MAX_LEN_FOR_MSA = 3000  # cap per-sequence length
too_many = len(seqs) > MAX_SEQS_FOR_MSA
too_long = any(len(s) > MAX_LEN_FOR_MSA for s in seqs)

if too_many or too_long:
    st.warning(
        f"Large input detected (n={len(seqs)}, max_len={max(len(s) for s in seqs)}). "
        f"MSA/Distance/Tree steps are disabled to keep the app responsive. "
        f"Consider subsetting sequences or lowering limits."
    )

# ⭐ NEW: require click to run heavy analyses
if not run:
    st.info("Adjust settings, then click **Run analysis** in the sidebar.")
    st.stop()

# ---------------------------
# GC content (lightweight)
# ---------------------------
gc_png = work_dir / "gc.png"
with st.spinner("Plotting GC%…"):
    plot_gc_content(recs, str(gc_png))
st.image(str(gc_png), caption="GC% per sequence")

# ---------------------------
# Pairwise alignment (fast)
# ---------------------------
if len(seqs) >= 2:
    s1 = st.selectbox("Sequence 1", names, index=0)
    s2 = st.selectbox("Sequence 2", names, index=min(1, len(names) - 1))
    with st.spinner("Aligning selected pair…"):
        if mode == "global":
            a, b, score = needleman_wunsch(recs[s1], recs[s2])
        else:
            a, b, score = smith_waterman(recs[s1], recs[s2])
    st.write(f"**Score:** {score}")
    st.code(a + "\n" + b, language="text")
    pf = work_dir / "pairwise.png"
    plot_alignment([a, b], [s1, s2], str(pf))
    st.image(str(pf), caption="Pairwise alignment (consensus match)")

# ---------------------------
# MSA + distance + tree (heavy)
# ---------------------------
if len(seqs) >= 2 and not (too_many or too_long):
    seqs_tuple = tuple(seqs)
    names_tuple = tuple(names)

    with st.spinner("Building MSA…"):
        aln = cached_msa(seqs_tuple)
    mf = work_dir / "msa.png"
    plot_alignment(aln, names, str(mf))
    st.image(str(mf), caption="MSA consensus view")

    with st.spinner(f"Computing distances ({model})…"):
        D = cached_distance(seqs_tuple, model)
    df = work_dir / "dist.png"
    plot_distance_heatmap(D, names, str(df), title=f"{model}-distance Heatmap")
    st.image(str(df), caption="Distance matrix")

    with st.spinner("Clustering tree…"):
        newick = cached_tree(names_tuple, D)
    (work_dir / "tree.newick").write_text(newick, encoding="utf-8")
    st.code(newick, language="text")

# ---------------------------
# Translation + codon usage (lightweight)
# ---------------------------
if seqs:
    prot = translate_dna(seqs[0], frame=0, stop_behavior="truncate")
    st.write("Protein (first sequence, frame 0):")
    st.code(prot, language="text")
    cf = work_dir / "codon.png"
    with st.spinner("Plotting codon usage…"):
        plot_codon_usage(seqs[0], str(cf))
    st.image(str(cf), caption="Codon usage (frame 0)")

# ---------------------------
# RNA folding (mid-weight)
# ---------------------------
if seqs:
    with st.spinner("Folding RNA…"):
        rna = dna_to_rna(seqs[0]).upper()
        dot = nussinov(rna, min_loop=0)
    st.write("RNA fold (dot-bracket):")
    st.code(dot, language="text")
    rf = work_dir / "rna_arcs.png"
    plot_rna_arcs(rna, dot, str(rf))
    st.image(str(rf), caption="RNA base-pair arcs")

st.markdown("---")
st.caption(
    "Tip: If it still feels slow, subset your FASTA or raise the MSA limits incrementally."
)

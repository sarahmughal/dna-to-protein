
# Beginner Bioinformatics Cheat Sheet (FASTA → Alignments → Trees → RNA)

## What is a FASTA file?
A FASTA is a simple text format for biological sequences. It has:
- A header line that starts with `>` followed by the sequence name/description.
- One or more lines of letters for the sequence (DNA uses A, C, G, T; RNA uses A, U, C, G; proteins are amino-acid letters).

**Example:**
```
>my_gene
ATGGCCATTGTAATGGGCCGCTGA
```

We read these into Python as `{header: sequence}` to process and analyze.

---

## GC Content (bar chart)
**What it shows:** The percentage of `G` and `C` letters in each DNA sequence.  
**Why it matters:** GC-rich DNA tends to be more thermally stable. Differences in GC% can hint at species differences, genome regions with bias, or sequencing artifacts.

---

## Pairwise Alignment & MSA (consensus match view)
**Pairwise alignment:** lines up two sequences to maximize matches and minimize penalties for mismatches and gaps.  
- **Global (Needleman–Wunsch):** aligns end-to-end—good when sequences are similar along their full length.  
- **Local (Smith–Waterman):** finds the best matching *subregion*—useful when only parts are similar.

**Multiple Sequence Alignment (MSA):** lines up *many* sequences together.  
**Consensus match plot:** bright cells = positions that match the consensus base for that column; dark = mismatches or gaps.  
**Why it matters:** Highly conserved (bright) columns may be functionally important; variable columns might be evolving faster.

---

## Distance Matrices & Heatmaps (p-distance / JC69)
**p-distance:** the fraction of differing aligned positions between two sequences.  
**Jukes–Cantor (JC69):** adjusts for unseen multiple substitutions at the same site.  
**Why it matters:** Distance matrices summarize “how different” sequences are pairwise. Heatmaps make the patterns obvious at a glance.

---

## Gene Trees (UPGMA, Newick format)
**UPGMA:** a simple clustering method that builds an ultrametric tree from distances.  
**Newick format:** a parenthesis-based text format to store trees, e.g., `(A:0.1,(B:0.05,C:0.05):0.05);`  
**Why it matters:** Trees visualize hypothesized evolutionary relationships. UPGMA is quick; for unequal evolutionary rates, try Neighbor-Joining later.

---

## Translation & Codon Usage
**DNA → RNA → Protein:** DNA with `T` is transcribed to RNA with `U`, then every 3 letters (codon) translate to an amino acid.  
**Codon usage plot:** counts of each codon in a sequence (frame 0).  
**Why it matters:** Organisms use some codons more than others; bias can reflect gene expression levels or tRNA availability.

---

## RNA Folding (dot-bracket + arc diagram)
**Nussinov algorithm:** predicts secondary structure by maximizing base pairs (A-U, G-C, and wobble G-U).  
**Dot-bracket:** paired positions show as `(` and `)`, unpaired as `.`  
**Arc diagram:** semicircles connect paired positions along the sequence.  
**Why it matters:** Secondary structure affects stability and function (e.g., stem-loops can regulate translation).

---

## General tips
- Normalize inputs (uppercase letters, remove whitespace).  
- Always check if sequences are related before interpreting trees.  
- Visuals help you see biological *signals* (conservation, divergence, structure) quickly.

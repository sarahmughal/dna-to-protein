# DNA → RNA → Protein: Alignments, Gene Trees, and RNA Folding (Beginner Friendly)

This beginner-level Python project showcases core bioinformatics tasks:

- Translate DNA → RNA → Protein
- Read and write FASTA
- Pairwise alignment (Needleman–Wunsch + Smith–Waterman)
- Simple multiple sequence alignment (center-star, progressive guide)
- Distance matrices (p-distance, Jukes–Cantor)
- UPGMA gene trees with Newick export
- RNA secondary structure (Nussinov algorithm)
- Dataset generation (random + simulated homologs)

Everything is pure-Python (no external deps required), with optional Biopython noted in comments.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
pytest
python -m src.cli --help
```

## Examples

Translate a DNA sequence and write protein:
```bash
python -m src.cli translate --dna ACTGATGGCTGA --frame 0 --stop truncate
```

Pairwise global alignment (Needleman–Wunsch):
```bash
python -m src.cli align --seq1 GATTACA --seq2 GCATGCU --mode global
```

Build a simple gene tree from a FASTA of related sequences:
```bash
python -m src.cli tree --fasta data/simulated_family/family4.fasta --model jc
```

Fold an RNA sequence (dot-bracket notation):
```bash
python -m src.cli fold --rna GCGCGAUUCGCG
```

Generate your own simulated homologous family:
```bash
python -m src.cli simulate --out data/my_family.fasta --n 6 --length 120 --gc 0.45 --mu 0.05
```

## Project layout

```
src/
  cli.py           # Command-line interface (argparse)
  io_utils.py      # FASTA read/write
  translate.py     # DNA→RNA→Protein
  align.py         # Needleman–Wunsch (global), Smith–Waterman (local), simple MSA
  distance.py      # p-distance, Jukes–Cantor
  tree.py          # UPGMA and Newick export
  rna_fold.py      # Nussinov secondary structure
data/
  example_single.fasta
  example_group.fasta
  simulated_family/...
tests/
  test_translate.py
  test_fasta.py
  test_align.py
  test_tree.py
.github/workflows/python-tests.yml  # CI with pytest
```

## License

MIT (see `LICENSE`).

---

*Generated on 2025-11-06T03:17:27.619636Z*

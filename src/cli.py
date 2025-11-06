import argparse, sys, os
from .translate import dna_to_rna, translate_dna
from .io_utils import read_fasta, write_fasta
from .align import needleman_wunsch, smith_waterman, center_star_msa
from .distance import distance_matrix
from .tree import upgma, to_newick
from .rna_fold import nussinov
from .simulate import simulate_family

def cmd_translate(args):
    prot = translate_dna(args.dna, frame=args.frame, stop_behavior=args.stop)
    rna = dna_to_rna(args.dna)
    print("RNA:", rna)
    print("Protein:", prot)

def cmd_fasta(args):
    recs = read_fasta(args.fasta)
    print(f"Loaded {len(recs)} records:")
    for h, s in recs.items():
        print(f"- {h}: {len(s)} bp")

def cmd_align(args):
    if args.mode == 'global':
        a,b,s = needleman_wunsch(args.seq1, args.seq2)
    else:
        a,b,s = smith_waterman(args.seq1, args.seq2)
    print(a)
    print(b)
    print("score:", s)

def cmd_msa(args):
    recs = read_fasta(args.fasta)
    seqs = list(recs.values())
    aln = center_star_msa(seqs)
    out = {h: aln[i] for i, h in enumerate(recs.keys())}
    write_fasta(args.out, out)
    print(f"Wrote MSA to {args.out}")

def cmd_tree(args):
    recs = read_fasta(args.fasta)
    names = list(recs.keys())
    seqs = list(recs.values())
    D = distance_matrix(seqs, model=args.model)
    root = upgma(names, D)
    newick = to_newick(root) + ";"
    if args.out:
        with open(args.out, 'w', encoding='utf-8') as f:
            f.write(newick + "\n")
        print(f"Wrote Newick to {args.out}")
    else:
        print(newick)

def cmd_fold(args):
    structure = nussinov(args.rna, min_loop=args.min_loop)
    print(args.rna)
    print(structure)

def cmd_simulate(args):
    fam = simulate_family(n=args.n, length=args.length, gc=args.gc, mu=args.mu)
    write_fasta(args.out, fam)
    print(f"Wrote simulated family to {args.out}")

def main(argv=None):
    p = argparse.ArgumentParser(prog='bio-portfolio', description='Beginner bioinformatics toolkit')
    sub = p.add_subparsers(required=True)

    t = sub.add_parser('translate', help='Translate DNA → RNA → protein')
    t.add_argument('--dna', required=True)
    t.add_argument('--frame', type=int, default=0, choices=[0,1,2])
    t.add_argument('--stop', default='truncate', choices=['truncate','keep','ignore'])
    t.set_defaults(func=cmd_translate)

    f = sub.add_parser('fasta', help='Read/inspect FASTA')
    f.add_argument('--fasta', required=True)
    f.set_defaults(func=cmd_fasta)

    a = sub.add_parser('align', help='Pairwise alignment')
    a.add_argument('--seq1', required=True)
    a.add_argument('--seq2', required=True)
    a.add_argument('--mode', choices=['global','local'], default='global')
    a.set_defaults(func=cmd_align)

    m = sub.add_parser('msa', help='Simple MSA (center-star)')
    m.add_argument('--fasta', required=True)
    m.add_argument('--out', required=True)
    m.set_defaults(func=cmd_msa)

    tr = sub.add_parser('tree', help='Build UPGMA tree and output Newick')
    tr.add_argument('--fasta', required=True)
    tr.add_argument('--model', choices=['p','jc'], default='jc')
    tr.add_argument('--out')
    tr.set_defaults(func=cmd_tree)

    rf = sub.add_parser('fold', help='RNA folding (Nussinov)')
    rf.add_argument('--rna', required=True)
    rf.add_argument('--min_loop', type=int, default=0)
    rf.set_defaults(func=cmd_fold)

    sim = sub.add_parser('simulate', help='Simulate homologous sequences into FASTA')
    sim.add_argument('--out', required=True)
    sim.add_argument('--n', type=int, default=5)
    sim.add_argument('--length', type=int, default=120)
    sim.add_argument('--gc', type=float, default=0.5)
    sim.add_argument('--mu', type=float, default=0.05)
    sim.set_defaults(func=cmd_simulate)

    args = p.parse_args(argv)
    args.func(args)

if __name__ == '__main__':
    main()

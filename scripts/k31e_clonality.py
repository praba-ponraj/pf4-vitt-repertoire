#!/usr/bin/env python3
"""
k31e_clonality.py  --  are K31E sequences independent, or one expanded clone?

A 2.13% sequence-level K31E frequency means something very different if the
6,522 sequences carrying it collapse into a handful of lineages.  This script
compares CDR-L3 diversity among K31E-bearing IGLV3-21 sequences against the
K31-germline IGLV3-21 background from the same runs.

Same germline anchoring as analyze_bulk_iglv321.py / k31e_null.py v3.

Reports, for K31E and for background separately:
  - unique CDR3 count, and unique CDR3 per 1000 sequences
  - share of sequences in the top 1 / top 10 CDR3s
  - Simpson / Shannon diversity
  - the top 15 CDR3s with counts
If K31E diversity is far below background, the percentage is a clone
frequency, not a mutation frequency.

Usage:
  python3 k31e_clonality.py Subject-HIP2/*.csv.gz --label HIP2
"""

import argparse
import math
import sys
from collections import Counter

import pandas as pd

ANCHOR = "NIGSKS"
K_OFF = 4
CHUNK = 200_000


def diversity(c):
    n = sum(c.values())
    if n == 0:
        return {}
    uniq = len(c)
    top1 = c.most_common(1)[0][1] / n
    top10 = sum(v for _, v in c.most_common(10)) / n
    simpson = sum((v / n) ** 2 for v in c.values())
    shannon = -sum((v / n) * math.log(v / n) for v in c.values())
    return dict(n=n, uniq=uniq, per1k=1000 * uniq / n, top1=top1,
                top10=top10, simpson=simpson, shannon=shannon,
                eff=math.exp(shannon))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--label", default="sample")
    ap.add_argument("--vgene", default="IGLV3-21")
    ap.add_argument("--cdr3col", default=None,
                    help="override CDR3 column (default: auto)")
    args = ap.parse_args()

    k31e, bg = Counter(), Counter()
    k31e_allele, bg_allele = Counter(), Counter()
    k31e_rows = []
    n_iglv = n_anch = 0
    cdr3col = args.cdr3col

    for path in args.files:
        try:
            reader = pd.read_csv(path, skiprows=1, chunksize=CHUNK,
                                 low_memory=False)
        except Exception as e:
            print(f"[skip] {path}: {e}", file=sys.stderr)
            continue
        got = 0
        for chunk in reader:
            if cdr3col is None:
                for c in ("cdr3_aa", "junction_aa", "cdr3"):
                    if c in chunk.columns:
                        cdr3col = c
                        break
                if cdr3col is None:
                    sys.exit(f"no CDR3 column; have: "
                             f"{[c for c in chunk.columns if 'cdr' in c.lower() or 'junc' in c.lower()]}")
                print(f"[info] using CDR3 column: {cdr3col}", file=sys.stderr)

            m = (chunk["v_call"].astype(str).str.split(",").str[0]
                 .str.split("*").str[0].eq(args.vgene))
            sub = chunk[m]
            if not len(sub):
                continue
            n_iglv += len(sub)
            got += len(sub)

            qs = sub["v_sequence_alignment_aa"].astype(str).tolist()
            gs = sub["v_germline_alignment_aa"].astype(str).tolist()
            cs = sub[cdr3col].astype(str).tolist()
            vs = sub["v_call"].astype(str).tolist()

            for q, g, c3, vc in zip(qs, gs, cs, vs):
                if q in ("nan", "") or g in ("nan", ""):
                    continue
                p = g.find(ANCHOR)
                if p < 0 or p + K_OFF >= min(len(g), len(q)):
                    continue
                if g[p + K_OFF] != "K":
                    continue
                n_anch += 1
                aa = q[p + K_OFF]
                al = vc.split(",")[0].split("*")[-1][:2]
                if c3 in ("nan", ""):
                    continue
                if aa == "E":
                    k31e[c3] += 1
                    k31e_allele[al] += 1
                    k31e_rows.append((c3, vc.split(",")[0]))
                elif aa == "K":
                    bg[c3] += 1
                    bg_allele[al] += 1
        print(f"[read] {path}: {got:,} {args.vgene}", file=sys.stderr)

    print(f"\n=== {args.label} ===")
    print(f"{args.vgene} sequences     : {n_iglv:,}")
    print(f"anchored at K31          : {n_anch:,}")

    dk, db = diversity(k31e), diversity(bg)
    if not dk or not db:
        sys.exit("insufficient data")

    print(f"\n{'':<28}{'K31E':>14}{'K31 germline':>16}")
    rows = [
        ("sequences", f"{dk['n']:,}", f"{db['n']:,}"),
        ("unique CDR3", f"{dk['uniq']:,}", f"{db['uniq']:,}"),
        ("unique CDR3 per 1000 seq", f"{dk['per1k']:.1f}", f"{db['per1k']:.1f}"),
        ("top 1 CDR3 share", f"{100*dk['top1']:.2f}%", f"{100*db['top1']:.2f}%"),
        ("top 10 CDR3 share", f"{100*dk['top10']:.2f}%", f"{100*db['top10']:.2f}%"),
        ("Simpson index", f"{dk['simpson']:.5f}", f"{db['simpson']:.5f}"),
        ("Shannon H", f"{dk['shannon']:.3f}", f"{db['shannon']:.3f}"),
        ("effective # clones", f"{dk['eff']:.0f}", f"{db['eff']:.0f}"),
    ]
    for a, b, c in rows:
        print(f"{a:<28}{b:>14}{c:>16}")

    print(f"\nratio of unique-CDR3-per-1000 (K31E / background): "
          f"{dk['per1k']/db['per1k']:.3f}")
    print("  ~1.0  -> K31E sequences are as diverse as background "
          "(independent events)")
    print("  <<1.0 -> K31E is concentrated in few lineages "
          "(clone frequency, not mutation frequency)")

    print(f"\n--- top 15 CDR3s among K31E sequences ---")
    for c3, n in k31e.most_common(15):
        print(f"  {n:>7,}  {100*n/dk['n']:>6.2f}%  {c3}")

    print(f"\n--- top 5 CDR3s in K31-germline background (for scale) ---")
    for c3, n in bg.most_common(5):
        print(f"  {n:>7,}  {100*n/db['n']:>6.2f}%  {c3}")

    print(f"\n--- allele assignment of K31E vs background ---")
    for a in sorted(set(k31e_allele) | set(bg_allele)):
        ke, be = k31e_allele.get(a, 0), bg_allele.get(a, 0)
        print(f"  *{a}: K31E {ke:>7,} ({100*ke/max(1,dk['n']):5.2f}%)   "
              f"background {be:>8,} ({100*be/max(1,db['n']):5.2f}%)")

    out = f"{args.label}_K31E_cdr3.tsv"
    with open(out, "w") as fh:
        fh.write("cdr3_aa\tcount\n")
        for c3, n in k31e.most_common():
            fh.write(f"{c3}\t{n}\n")
    print(f"\nwrote CDR3 table -> {out}")


if __name__ == "__main__":
    main()

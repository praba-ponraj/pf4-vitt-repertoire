#!/usr/bin/env python3
"""
k31e_allele_runs.py  --  final light-chain QC.

Two jobs, one pass:

 (A) EXACT top-N sensitivity analysis.  Reports K31E prevalence after
     removing the 1 / 5 / 10 / 25 most abundant K31E-associated CDR-L3
     sequences, computed from raw counts against the original anchored
     denominator -- not back-computed from rounded percentages.

 (B) RUN-STRATIFIED allele association.  Per run and pooled, the
     *02 vs *04 odds ratio for K31E versus K31-germline, plus a
     Mantel-Haenszel pooled OR that controls for run.  This is the test
     of whether the HIP2 *02 enrichment is biological or an artefact of
     the three high-yield runs.

Same germline anchoring as analyze_bulk_iglv321.py (NIGSKS, K at +4).

Usage:
  python3 k31e_allele_runs.py Subject-HIP2/*.csv.gz --label HIP2_all9
  python3 k31e_allele_runs.py Subject-HIP1/*.csv.gz --label HIP1
"""

import argparse
import math
import os
import sys
from collections import Counter

import pandas as pd

ANCHOR = "NIGSKS"
K_OFF = 4
CHUNK = 200_000


def woolf_ci(a, b, c, d):
    """OR with Woolf 95% CI for 2x2: [[a,b],[c,d]]."""
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    orr = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return orr, math.exp(math.log(orr) - 1.96 * se), \
        math.exp(math.log(orr) + 1.96 * se)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--label", default="sample")
    ap.add_argument("--vgene", default="IGLV3-21")
    ap.add_argument("--topn", type=int, nargs="+", default=[1, 5, 10, 25])
    args = ap.parse_args()

    per_run = []          # (run, n_iglv, n_anch, k31e, strata dict)
    all_cdr3 = Counter()  # pooled K31E CDR-L3 counts
    tot_anch = tot_k31e = tot_iglv = tot_noanchor = 0
    pooled = Counter()    # ('k31e'|'bg', allele) -> count

    for path in args.files:
        run = os.path.basename(path).split("_")[0]
        n_iglv = n_anch = n_k31e = n_noanchor = 0
        strata = Counter()
        try:
            reader = pd.read_csv(path, skiprows=1, chunksize=CHUNK,
                                 low_memory=False)
        except Exception as e:
            print(f"[skip] {path}: {e}", file=sys.stderr)
            continue

        cdr3col = None
        for chunk in reader:
            if cdr3col is None:
                for c in ("cdr3_aa", "junction_aa"):
                    if c in chunk.columns:
                        cdr3col = c
                        break
            m = (chunk["v_call"].astype(str).str.split(",").str[0]
                 .str.split("*").str[0].eq(args.vgene))
            sub = chunk[m]
            if not len(sub):
                continue
            n_iglv += len(sub)

            qs = sub["v_sequence_alignment_aa"].astype(str).tolist()
            gs = sub["v_germline_alignment_aa"].astype(str).tolist()
            vs = sub["v_call"].astype(str).tolist()
            cs = (sub[cdr3col].astype(str).tolist() if cdr3col
                  else [""] * len(sub))

            for q, g, vc, c3 in zip(qs, gs, vs, cs):
                if q in ("nan", "") or g in ("nan", ""):
                    n_noanchor += 1
                    continue
                p = g.find(ANCHOR)
                if p < 0 or p + K_OFF >= min(len(g), len(q)) \
                        or g[p + K_OFF] != "K":
                    n_noanchor += 1
                    continue
                n_anch += 1
                aa = q[p + K_OFF]
                al = vc.split(",")[0].split("*")[-1][:2]
                if aa == "E":
                    n_k31e += 1
                    strata[("k31e", al)] += 1
                    pooled[("k31e", al)] += 1
                    if c3 not in ("nan", ""):
                        all_cdr3[c3] += 1
                elif aa == "K":
                    strata[("bg", al)] += 1
                    pooled[("bg", al)] += 1

        per_run.append((run, n_iglv, n_anch, n_k31e, strata))
        tot_iglv += n_iglv
        tot_anch += n_anch
        tot_k31e += n_k31e
        tot_noanchor += n_noanchor
        print(f"[read] {run}: {n_iglv:,} {args.vgene}, "
              f"{n_k31e:,} K31E", file=sys.stderr)

    if not tot_anch:
        sys.exit("nothing anchored")

    print(f"\n=== {args.label} ===")
    print(f"{args.vgene} sequences : {tot_iglv:,}")
    print(f"anchored at K31      : {tot_anch:,}")
    print(f"anchor failures      : {tot_noanchor:,} "
          f"({100*tot_noanchor/tot_iglv:.4f}%)")
    print(f"K31E                 : {tot_k31e:,} "
          f"({100*tot_k31e/tot_anch:.4f}% of anchored)")

    # ---------- (A) exact top-N sensitivity -------------------------------
    print(f"\n--- (A) exact top-N CDR-L3 removal sensitivity ---")
    print(f"distinct K31E CDR-L3 sequences: {len(all_cdr3):,}")
    print(f"{'removed':>10} {'seqs removed':>14} {'K31E left':>12} "
          f"{'prevalence':>12}")
    print(f"{'none':>10} {0:>14,} {tot_k31e:>12,} "
          f"{100*tot_k31e/tot_anch:>11.4f}%")
    top = all_cdr3.most_common(max(args.topn))
    for n in sorted(args.topn):
        removed = sum(v for _, v in top[:n])
        left = tot_k31e - removed
        print(f"{'top ' + str(n):>10} {removed:>14,} {left:>12,} "
              f"{100*left/tot_anch:>11.4f}%")
    print("\nReport as a sensitivity analysis, not a replacement estimate.")

    # ---------- (B) run-stratified allele association ---------------------
    print(f"\n--- (B) *02 vs *04 association with K31E, by run ---")
    print(f"{'run':>12} {'K31E 02':>9} {'K31E 04':>9} {'bg 02':>10} "
          f"{'bg 04':>10} {'OR':>8} {'95% CI':>20}")
    mh_num = mh_den = 0.0
    usable = 0
    for run, n_iglv, n_anch, n_k31e, s in per_run:
        a = s[("k31e", "02")] + s[("k31e", "03")]
        b = s[("k31e", "04")] + s[("k31e", "01")]
        c = s[("bg", "02")] + s[("bg", "03")]
        d = s[("bg", "04")] + s[("bg", "01")]
        if min(a + b, c + d) == 0 or (a + b) < 10:
            print(f"{run:>12} {a:>9,} {b:>9,} {c:>10,} {d:>10,} "
                  f"{'--':>8} {'(too few K31E)':>20}")
            continue
        orr, lo, hi = woolf_ci(a, b, c, d)
        n = a + b + c + d
        mh_num += a * d / n
        mh_den += b * c / n
        usable += 1
        print(f"{run:>12} {a:>9,} {b:>9,} {c:>10,} {d:>10,} "
              f"{orr:>8.3f} {f'{lo:.2f}-{hi:.2f}':>20}")

    A = pooled[("k31e", "02")] + pooled[("k31e", "03")]
    B = pooled[("k31e", "04")] + pooled[("k31e", "01")]
    C = pooled[("bg", "02")] + pooled[("bg", "03")]
    D = pooled[("bg", "04")] + pooled[("bg", "01")]
    if min(A, B, C, D) > 0:
        orr, lo, hi = woolf_ci(A, B, C, D)
        print(f"\n{'CRUDE pooled':>12} {A:>9,} {B:>9,} {C:>10,} {D:>10,} "
              f"{orr:>8.3f} {f'{lo:.2f}-{hi:.2f}':>20}")
    if mh_den > 0 and usable > 1:
        print(f"{'Mantel-Haenszel (run-stratified)':>44}: "
              f"OR = {mh_num/mh_den:.3f}   ({usable} runs contributing)")
        print("\nIf the MH OR tracks the crude OR, the association is not "
              "explained by run;\nif it collapses toward 1.0, it was "
              "confounded by run composition.")
    elif usable <= 1:
        print("\nOnly one run has enough K31E events -- run stratification "
              "is not possible,\nand the association cannot be separated "
              "from that run.")

    print(f"\n--- allele composition detail (pooled) ---")
    for al in ("01", "02", "03", "04"):
        k, bgc = pooled[("k31e", al)], pooled[("bg", al)]
        kt = sum(pooled[("k31e", x)] for x in ("01", "02", "03", "04"))
        bt = sum(pooled[("bg", x)] for x in ("01", "02", "03", "04"))
        print(f"  *{al}: K31E {k:>8,} ({100*k/max(1,kt):6.2f}%)   "
              f"background {bgc:>9,} ({100*bgc/max(1,bt):6.2f}%)")


if __name__ == "__main__":
    main()

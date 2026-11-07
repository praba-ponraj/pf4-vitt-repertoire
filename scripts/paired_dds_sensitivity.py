#!/usr/bin/env python3
"""
paired_dds_sensitivity.py  --  sequence-only sensitivity analysis for Table 1C.

The DDS(D) proxy is a disjunction:
    IGLV3-21  AND  (allele in {02,03}  OR  cdr2_aa_light == "DDS")
Of the 9,394 proxy pairs, 7,424 (79.0%) carry observed DDS at CDR-L2 and
1,970 (21.0%) qualify on the allele call alone.  This script recomputes the
Table 1C features in each subset so the published percentages can be shown
to be robust (or not) to allele assignment.

Feature definitions are copied verbatim from analyze_paired_iglv321.py:
  K31E          : germline-anchored (NIGSKS in v_germline_alignment_aa_light,
                  K at +4), glutamate in the query at that position
  IGHJ3         : gene of first j_call_heavy == IGHJ3
  acidic HCDR3  : cdr3_aa_heavy contains >= 3 D/E
  VITT IGHV     : gene of first v_call_heavy in the five reported genes
  five-feature  : K31E AND VITT IGHV AND IGHJ3 AND acidic

Nothing else is computed.

Usage:
  python3 paired_dds_sensitivity.py paried-naive-human/*.csv.gz
"""

import argparse
import csv
import gzip
import sys
from collections import Counter

GERMLINE_CONTEXT = "NIGSKS"
GERMLINE_K_OFFSET = 4

VITT_IGHV = {
    "IGHV3-7",
    "IGHV5-51",
    "IGHV2-26",
    "IGHV3-30",
    "IGHV7-4-1",
}

FEATURES = [
    ("k31e", "K31E"),
    ("ighj3", "IGHJ3"),
    ("acidic_hcdr3", "acidic HCDR3 >=3 D/E"),
    ("vitt_ighv", "VITT-reported IGHV"),
    ("vitt_ighv_ighj3", "VITT IGHV + IGHJ3"),
    ("complete_five_feature", "complete five-feature"),
]


def first_call(value):
    return (value or "").split(",")[0].strip()


def gene_from_call(value):
    call = first_call(value)
    return call.split("*")[0] if call else ""


def allele_from_call(value):
    call = first_call(value)
    if "*" not in call:
        return ""
    return call.split("*", 1)[1]


def aligned_k31e(row):
    seq = row.get("v_sequence_alignment_aa_light", "") or ""
    germ = row.get("v_germline_alignment_aa_light", "") or ""
    p = germ.find(GERMLINE_CONTEXT)
    if p < 0:
        return None
    kpos = p + GERMLINE_K_OFFSET
    if kpos >= len(germ) or kpos >= len(seq):
        return None
    if germ[kpos] != "K":
        return None
    return seq[kpos] == "E"


def acidic_hcdr3(row):
    h3 = row.get("cdr3_aa_heavy", "") or ""
    return (h3.count("D") + h3.count("E")) >= 3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    args = ap.parse_args()

    # three populations: full proxy, DDS observed, allele-only
    pops = {"full": Counter(), "dds": Counter(), "allele_only": Counter()}

    for path in args.files:
        try:
            with gzip.open(path, "rt", newline="") as fh:
                next(fh)                       # OAS metadata line
                reader = csv.DictReader(fh)
                need = {
                    "v_call_light",
                    "v_sequence_alignment_aa_light",
                    "v_germline_alignment_aa_light",
                    "cdr2_aa_light",
                    "v_call_heavy",
                    "j_call_heavy",
                    "cdr3_aa_heavy",
                }
                missing = need - set(reader.fieldnames or [])
                if missing:
                    print(f"[skip] {path}: missing {sorted(missing)}",
                          file=sys.stderr)
                    continue

                for row in reader:
                    lv_call = first_call(row.get("v_call_light", ""))
                    if gene_from_call(lv_call) != "IGLV3-21":
                        continue

                    allele = allele_from_call(lv_call)
                    lcdr2 = (row.get("cdr2_aa_light", "") or "").strip()

                    by_allele = allele in {"02", "03"}
                    by_dds = lcdr2 == "DDS"
                    if not (by_allele or by_dds):
                        continue

                    k31 = aligned_k31e(row)
                    k31e = k31 is True

                    hv = gene_from_call(row.get("v_call_heavy", ""))
                    hj = gene_from_call(row.get("j_call_heavy", ""))
                    is_vitt = hv in VITT_IGHV
                    is_j3 = hj == "IGHJ3"
                    is_acid = acidic_hcdr3(row)

                    which = ["full"]
                    which.append("dds" if by_dds else "allele_only")

                    for key in which:
                        c = pops[key]
                        c["n"] += 1
                        if k31 is None:
                            c["k31_unmapped"] += 1
                        if k31e:
                            c["k31e"] += 1
                        if is_j3:
                            c["ighj3"] += 1
                        if is_acid:
                            c["acidic_hcdr3"] += 1
                        if is_vitt:
                            c["vitt_ighv"] += 1
                        if is_vitt and is_j3:
                            c["vitt_ighv_ighj3"] += 1
                        if k31e and is_vitt and is_j3 and is_acid:
                            c["complete_five_feature"] += 1
        except Exception as e:
            print(f"[skip] {path}: {e}", file=sys.stderr)

    def pct(a, b):
        return 100.0 * a / b if b else 0.0

    nf = pops["full"]["n"]
    nd = pops["dds"]["n"]
    na = pops["allele_only"]["n"]

    print("\n=== Table 1C sequence-only sensitivity analysis ===\n")
    print(f"full DDS(D) proxy            : {nf:,}")
    print(f"  observed DDS at CDR-L2     : {nd:,} ({pct(nd, nf):.2f}%)")
    print(f"  allele call (*02/*03) only : {na:,} ({pct(na, nf):.2f}%)")
    print(f"\nK31 unmapped  full {pops['full']['k31_unmapped']:,}   "
          f"DDS {pops['dds']['k31_unmapped']:,}   "
          f"allele-only {pops['allele_only']['k31_unmapped']:,}")

    print(f"\n{'feature':<24}{'full proxy':>20}{'DDS observed':>20}"
          f"{'allele-only':>20}")
    print(f"{'':<24}{f'(n={nf:,})':>20}{f'(n={nd:,})':>20}"
          f"{f'(n={na:,})':>20}")
    print("-" * 84)
    for key, label in FEATURES:
        f_ = pops["full"][key]
        d_ = pops["dds"][key]
        a_ = pops["allele_only"][key]
        print(f"{label:<24}"
              f"{f'{f_:,} ({pct(f_, nf):.2f}%)':>20}"
              f"{f'{d_:,} ({pct(d_, nd):.2f}%)':>20}"
              f"{f'{a_:,} ({pct(a_, na):.2f}%)':>20}")

    print("\nRead: if the DDS-observed column tracks the full-proxy column,")
    print("the paired result is robust to allele assignment and can be")
    print("stated as such in one clause. Material divergence between the")
    print("DDS-observed and allele-only columns would mean the published")
    print("percentages are a mixture of two populations.")


if __name__ == "__main__":
    main()

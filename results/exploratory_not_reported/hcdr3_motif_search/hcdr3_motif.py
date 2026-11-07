#!/usr/bin/env python3
"""
hcdr3_motif.py (v3)  --  VITT HCDR3 motif audit, on the Table 1C denominator.

Replicates the DDS(D) proxy exactly as analyze_paired_iglv321.py defines it:
    light-chain gene == IGLV3-21
    AND (light-chain allele in {02,03}  OR  cdr2_aa_light == "DDS")
(exact string equality on CDR-L2, applied after the gene filter).

Reports three denominators so the motif result is comparable to the paper:
    (1) all usable paired HCDR3s
    (2) IGLV3-21 pairs
    (3) DDS(D)-proxy pairs            <-- the Table 1C population

Also reports the proxy BREAKDOWN -- how many proxy pairs qualify by observed
DDS only, by *02/*03 call only, or by both -- which determines whether the
paired result can be stated as sequence-defined rather than allele-dependent.

MOTIF PROVENANCE:
  Wang JJ et al., Blood 2022;140(15):1738-42 -> published consensus "G/NLED"
    (4 residues; regex [GN]LED).
  Ivanov et al., JACS 2023;145:25203-13 -> one VITT mAb with HCDR3 GLEDAFD,
    described as resembling "(G/N)LED(A/V/T)FD".  That 7-residue string is
    Ivanov's EXTENSION, and its DAFD tail overlaps IGHJ3 germline
    (DAFDI/DAFDV), so matches to it largely restate IGHJ3 usage.

NULLS are simple amino-acid-composition expectations computed from the
observed HCDR3s of the denominator in question: product of position-class
frequencies x approximate number of placement sites x N.  HCDR3 residues are
neither independent nor positionally homogeneous, so obs/exp is DESCRIPTIVE.
No P value is computed and none should be quoted.

A second null is computed within IGHJ3-using sequences only, so motif
enrichment can be read conditional on IGHJ3 rather than alongside it.

Usage:
  python3 hcdr3_motif.py jaffe_paired/*.csv.gz --label paired
  python3 hcdr3_motif.py jaffe_paired/*.csv.gz paried-naive-human/*.csv.gz \
      --label paired_all_dirs
"""

import argparse
import re
import sys
from collections import Counter

import pandas as pd

CHUNK = 100_000
AA = set("ACDEFGHIKLMNPQRSTVWY")

TIERS = [
    ("T1_Wang2022_[GN]LED",      re.compile(r"[GN]LED"),
     [["G", "N"], ["L"], ["E"], ["D"]]),
    ("T2_Ivanov_[GN]LED[AVT]FD", re.compile(r"[GN]LED[AVT]FD"),
     [["G", "N"], ["L"], ["E"], ["D"], ["A", "V", "T"], ["F"], ["D"]]),
    ("T3_LED",                   re.compile(r"LED"),
     [["L"], ["E"], ["D"]]),
    ("T4_ED",                    re.compile(r"ED"),
     [["E"], ["D"]]),
    ("T5_acidic_ge3_DE",         None, None),
]

VITT_IGHV = {"IGHV3-7", "IGHV5-51", "IGHV2-26", "IGHV3-30", "IGHV7-4-1"}


def first_call(v):
    return (str(v) if v is not None else "").split(",")[0].strip()


def gene_of(v):
    c = first_call(v)
    return c.split("*")[0] if c else ""


def allele_of(v):
    c = first_call(v)
    return c.split("*", 1)[1] if "*" in c else ""


class Bucket:
    """Accumulates HCDR3s and motif hits for one denominator."""

    def __init__(self, name):
        self.name = name
        self.n = 0
        self.comp = Counter()
        self.len_sum = 0
        self.tier = Counter()
        self.j3 = 0
        self.j3_comp = Counter()
        self.j3_len_sum = 0
        self.tier_j3 = Counter()
        self.tier_vv = Counter()
        self.examples = {t[0]: [] for t in TIERS}

    def add(self, h3, is_j3, is_vv):
        self.n += 1
        self.comp.update(h3)
        self.len_sum += len(h3)
        if is_j3:
            self.j3 += 1
            self.j3_comp.update(h3)
            self.j3_len_sum += len(h3)
        for name, rx, _ in TIERS:
            hit = ((h3.count("D") + h3.count("E")) >= 3 if rx is None
                   else bool(rx.search(h3)))
            if not hit:
                continue
            self.tier[name] += 1
            if is_j3:
                self.tier_j3[name] += 1
            if is_vv:
                self.tier_vv[name] += 1
            if len(self.examples[name]) < 20:
                self.examples[name].append(h3)

    def _null(self, comp, len_sum, n, classes):
        tot = sum(comp.values())
        if not tot or not n or not classes:
            return None
        f = {a: v / tot for a, v in comp.items()}
        p = 1.0
        for cl in classes:
            p *= sum(f.get(a, 0) for a in cl)
        mean_len = len_sum / n
        sites = max(0.0, mean_len - len(classes) + 1)
        return p * sites * n

    def report(self):
        if not self.n:
            print(f"\n=== {self.name}: no sequences ===")
            return
        print(f"\n{'='*70}")
        print(f"=== {self.name} ===")
        print(f"usable HCDR3s     : {self.n:,}")
        print(f"mean HCDR3 length : {self.len_sum/self.n:.2f}")
        print(f"IGHJ3 usage       : {self.j3:,} ({100*self.j3/self.n:.3f}%)")

        print(f"\n  all sequences in this denominator")
        print(f"  {'tier':<26}{'obs':>10}{'obs %':>9}"
              f"{'expected':>11}{'obs/exp':>9}")
        for name, rx, cls in TIERS:
            obs = self.tier[name]
            exp = self._null(self.comp, self.len_sum, self.n, cls)
            e = f"{exp:,.0f}" if exp else "--"
            r = f"{obs/exp:.2f}" if exp else "--"
            print(f"  {name:<26}{obs:>10,}{100*obs/self.n:>8.3f}%{e:>11}{r:>9}")

        if self.j3:
            print(f"\n  conditional on IGHJ3  (n = {self.j3:,})")
            print(f"  {'tier':<26}{'obs':>10}{'obs %':>9}"
                  f"{'expected':>11}{'obs/exp':>9}")
            for name, rx, cls in TIERS:
                obs = self.tier_j3[name]
                exp = self._null(self.j3_comp, self.j3_len_sum, self.j3, cls)
                e = f"{exp:,.0f}" if exp else "--"
                r = f"{obs/exp:.2f}" if exp else "--"
                print(f"  {name:<26}{obs:>10,}{100*obs/self.j3:>8.3f}%"
                      f"{e:>11}{r:>9}")

        print(f"\n  tiers crossed with VITT-reported IGHV")
        for name, _, _ in TIERS:
            print(f"    {name:<26}{self.tier_vv[name]:>8,}")

        for name, _, _ in TIERS:
            if self.examples[name] and name.startswith(("T1", "T2")):
                print(f"\n  example HCDR3s, {name}:")
                for s in self.examples[name][:20]:
                    print(f"    {s}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--label", default="paired")
    args = ap.parse_args()

    b_all = Bucket("(1) ALL usable paired HCDR3s")
    b_lv = Bucket("(2) IGLV3-21 pairs")
    b_px = Bucket("(3) DDS(D)-proxy pairs  [Table 1C denominator]")

    n_pairs = 0
    px_break = Counter()   # 'dds_only' | 'allele_only' | 'both'
    missing_cols = set()

    for path in args.files:
        try:
            reader = pd.read_csv(path, skiprows=1, chunksize=CHUNK,
                                 low_memory=False)
        except Exception as e:
            print(f"[skip] {path}: {e}", file=sys.stderr)
            continue
        got = 0
        for chunk in reader:
            need = ["v_call_light", "cdr2_aa_light", "v_call_heavy",
                    "j_call_heavy", "cdr3_aa_heavy"]
            miss = [c for c in need if c not in chunk.columns]
            if miss:
                missing_cols.update(miss)
                break
            n_pairs += len(chunk)
            got += len(chunk)

            for lv, l2, hv, hj, h3 in zip(
                    chunk["v_call_light"], chunk["cdr2_aa_light"],
                    chunk["v_call_heavy"], chunk["j_call_heavy"],
                    chunk["cdr3_aa_heavy"]):
                h3 = str(h3)
                if len(h3) < 4 or not set(h3) <= AA:
                    continue
                is_j3 = gene_of(hj) == "IGHJ3"
                is_vv = gene_of(hv) in VITT_IGHV
                b_all.add(h3, is_j3, is_vv)

                if gene_of(lv) != "IGLV3-21":
                    continue
                b_lv.add(h3, is_j3, is_vv)

                al = allele_of(lv)
                l2s = (str(l2) if l2 is not None else "").strip()
                by_allele = al in {"02", "03"}
                by_dds = l2s == "DDS"
                if not (by_allele or by_dds):
                    continue
                b_px.add(h3, is_j3, is_vv)
                px_break["both" if (by_allele and by_dds)
                         else ("allele_only" if by_allele else "dds_only")] += 1
        print(f"[read] {path}: {got:,} pairs", file=sys.stderr)

    if missing_cols:
        print(f"\n[warn] some files lacked columns: "
              f"{sorted(missing_cols)}", file=sys.stderr)

    print(f"\n########  {args.label}  ########")
    print(f"paired records read : {n_pairs:,}")

    tot_px = sum(px_break.values())
    print(f"\n--- DDS(D) proxy breakdown (n = {tot_px:,}) ---")
    for k, lab in (("dds_only", "observed DDS at CDR-L2 only"),
                   ("allele_only", "*02/*03 allele call only"),
                   ("both", "both criteria")):
        v = px_break[k]
        print(f"  {lab:<32}{v:>9,}  ({100*v/tot_px:.2f}%)"
              if tot_px else f"  {lab:<32}{v:>9,}")
    seq_defined = px_break["dds_only"] + px_break["both"]
    if tot_px:
        print(f"\n  sequence-defined (DDS observed): {seq_defined:,} "
              f"({100*seq_defined/tot_px:.2f}%)")
        print(f"  allele-dependent only          : "
              f"{px_break['allele_only']:,} "
              f"({100*px_break['allele_only']/tot_px:.2f}%)")
        print("\n  If allele-only is a small minority, the paired result can be")
        print("  stated as sequence-defined with the disjunction as a")
        print("  sensitivity analysis.")

    for b in (b_all, b_lv, b_px):
        b.report()

    print(f"\n{'='*70}")
    print("Reading notes:")
    print("  * obs/exp is DESCRIPTIVE. The null assumes positional")
    print("    independence of HCDR3 residues, which does not hold.")
    print("    Quote no P value.")
    print("  * T1 is Wang 2022's published 4-residue consensus; a 4-mer")
    print("    occurs by chance, so read obs/exp, not obs %.")
    print("  * T2's DAFD tail overlaps IGHJ3 germline, so T2 largely")
    print("    restates IGHJ3 usage; the IGHJ3-conditional block is the")
    print("    one to read for it.")
    print("  * T5 is the manuscript's existing criterion, for comparison.")


if __name__ == "__main__":
    main()

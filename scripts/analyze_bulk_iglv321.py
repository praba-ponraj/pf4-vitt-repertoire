#!/usr/bin/env python3

"""
Analyze IGLV3-21 features in OAS bulk light-chain repertoires.

K31E definition:
Glutamate at the sequence position aligned to the IGLV3-21 germline
lysine corresponding to IMGT position 37 (K31 in Wang et al.).

This script intentionally does NOT use:
- the NIGSES substring
- a fixed CDR-L1 index
"""

import argparse
import csv
import gzip
import glob
import os
from collections import Counter


GERMLINE_CONTEXT = "NIGSKS"
GERMLINE_K_OFFSET = 4


def first_v_call(v_call):
    """Return first V assignment if multiple calls are present."""
    return (v_call or "").split(",")[0].strip()


def allele_from_v_call(v_call):
    v = first_v_call(v_call)
    if "*" not in v:
        return ""
    return v.split("*", 1)[1]


def aligned_k31e(row):
    """
    Determine K31E from V-region amino-acid alignment.

    Returns:
        True  = aligned residue is E
        False = aligned residue is not E
        None  = germline position could not be mapped
    """
    seq = row.get("v_sequence_alignment_aa", "") or ""
    germ = row.get("v_germline_alignment_aa", "") or ""

    p = germ.find(GERMLINE_CONTEXT)

    if p < 0:
        return None

    kpos = p + GERMLINE_K_OFFSET

    if kpos >= len(germ) or kpos >= len(seq):
        return None

    # Confirm that the anchor really is germline K.
    if germ[kpos] != "K":
        return None

    return seq[kpos] == "E"


def analyze_file(path):
    stats = Counter()

    with gzip.open(path, "rt", newline="") as handle:
        # OAS bulk files contain a metadata JSON record before the CSV header.
        next(handle)
        reader = csv.DictReader(handle)

        required = {
            "v_call",
            "v_sequence_alignment_aa",
            "v_germline_alignment_aa",
        }
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise RuntimeError(
                f"{path}: missing required columns: {sorted(missing)}"
            )

        for row in reader:
            stats["total_sequences"] += 1

            v = first_v_call(row.get("v_call", ""))

            if v.split("*")[0] != "IGLV3-21":
                continue

            stats["iglv321"] += 1

            allele = allele_from_v_call(v)
            if allele in {"02", "03"}:
                stats["allele_02_03"] += 1

            result = aligned_k31e(row)

            if result is None:
                stats["k31_unmapped"] += 1
            elif result:
                stats["k31e"] += 1
            else:
                stats["k31_non_e"] += 1

    return stats


def pct(n, d):
    return 100.0 * n / d if d else 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_dir",
        help="Directory containing Subject-HIP1, Subject-HIP2, Subject-HIP3",
    )
    parser.add_argument(
        "--output",
        default="bulk_iglv321_results.tsv",
        help="Output TSV file",
    )
    args = parser.parse_args()

    subjects = ["Subject-HIP1", "Subject-HIP2", "Subject-HIP3"]

    subject_stats = {}
    grand = Counter()
    rows = []

    for subject in subjects:
        pattern = os.path.join(args.input_dir, subject, "*.csv.gz")
        files = sorted(glob.glob(pattern))

        if not files:
            raise RuntimeError(f"No files found for {subject}: {pattern}")

        combined = Counter()

        for path in files:
            s = analyze_file(path)
            combined.update(s)

            rows.append({
                "level": "run",
                "subject": subject,
                "run": os.path.basename(path),
                **s,
            })

        subject_stats[subject] = combined
        grand.update(combined)

        rows.append({
            "level": "subject",
            "subject": subject,
            "run": "",
            **combined,
        })

    rows.append({
        "level": "total",
        "subject": "Total",
        "run": "",
        **grand,
    })

    fields = [
        "level",
        "subject",
        "run",
        "total_sequences",
        "iglv321",
        "iglv321_pct",
        "k31e",
        "k31e_pct",
        "k31_unmapped",
        "allele_02_03",
        "allele_02_03_pct",
    ]

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()

        for r in rows:
            d = dict(r)

            d["iglv321_pct"] = f"{pct(d.get('iglv321', 0), d.get('total_sequences', 0)):.4f}"
            d["k31e_pct"] = f"{pct(d.get('k31e', 0), d.get('iglv321', 0)):.4f}"
            d["allele_02_03_pct"] = (
                f"{pct(d.get('allele_02_03', 0), d.get('iglv321', 0)):.4f}"
            )

            writer.writerow({k: d.get(k, 0) for k in fields})

    print("\nSUBJECT RESULTS")
    print("-" * 72)

    for subject in subjects:
        s = subject_stats[subject]

        print(
            f"{subject}: "
            f"IGLV3-21={s['iglv321']:,}; "
            f"K31E={s['k31e']:,} "
            f"({pct(s['k31e'], s['iglv321']):.3f}%); "
            f"*02/*03={pct(s['allele_02_03'], s['iglv321']):.1f}%; "
            f"unmapped={s['k31_unmapped']:,}"
        )

    print("\nTOTAL")
    print("-" * 72)

    print(f"Total sequences: {grand['total_sequences']:,}")
    print(f"IGLV3-21:       {grand['iglv321']:,}")
    print(
        f"K31E:           {grand['k31e']:,} "
        f"({pct(grand['k31e'], grand['iglv321']):.3f}%)"
    )
    print(
        f"*02/*03:        {grand['allele_02_03']:,} "
        f"({pct(grand['allele_02_03'], grand['iglv321']):.1f}%)"
    )
    print(f"K31 unmapped:   {grand['k31_unmapped']:,}")

    print(f"\nResults written to: {args.output}")


if __name__ == "__main__":
    main()

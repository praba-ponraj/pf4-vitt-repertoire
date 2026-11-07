#!/usr/bin/env python3

"""
Analyze VITT-associated antibody features in paired OAS repertoires.

K31E definition:
Glutamate at the sequence position aligned to the IGLV3-21 germline
lysine in the NIGSKS context (K31 in Wang et al.). The germline alignment
is searched per sequence; no fixed query index or IMGT coordinate is used.

DDS(D) proxy:
IGLV3-21 pair assigned to *02/*03 OR containing DDS in CDR-L2.

This script intentionally does NOT use:
- NIGSES substring matching
- a fixed CDR-L1 index
- GSES/GRES exact matching
"""

import argparse
import csv
import gzip
import glob
import os
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
    """
    Determine light-chain K31E using the V-region amino-acid alignment.

    Returns:
        True  = aligned residue is E
        False = aligned residue is not E
        None  = germline K position could not be mapped
    """
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
    """Return True when heavy-chain CDR3 contains >=3 D/E residues."""
    h3 = row.get("cdr3_aa_heavy", "") or ""
    return (h3.count("D") + h3.count("E")) >= 3


def analyze_file(path):
    stats = Counter()

    with gzip.open(path, "rt", newline="") as handle:
        # OAS paired files contain metadata JSON before the CSV header.
        next(handle)
        reader = csv.DictReader(handle)

        required = {
            "v_call_light",
            "v_sequence_alignment_aa_light",
            "v_germline_alignment_aa_light",
            "cdr2_aa_light",
            "v_call_heavy",
            "j_call_heavy",
            "cdr3_aa_heavy",
        }

        missing = required - set(reader.fieldnames or [])
        if missing:
            raise RuntimeError(
                f"{path}: missing required columns: {sorted(missing)}"
            )

        for row in reader:
            stats["total_pairs"] += 1

            lv_call = first_call(row.get("v_call_light", ""))
            lv_gene = gene_from_call(lv_call)

            if lv_gene != "IGLV3-21":
                continue

            stats["iglv321"] += 1

            k31 = aligned_k31e(row)

            if k31 is None:
                stats["k31_unmapped_iglv321"] += 1
            elif k31:
                stats["k31e_all_iglv321"] += 1

            allele = allele_from_call(lv_call)
            lcdr2 = (row.get("cdr2_aa_light", "") or "").strip()

            dds_proxy = allele in {"02", "03"} or lcdr2 == "DDS"

            if not dds_proxy:
                continue

            stats["dds_proxy"] += 1

            if k31 is None:
                stats["k31_unmapped_proxy"] += 1

            k31e = k31 is True

            hv_gene = gene_from_call(row.get("v_call_heavy", ""))
            hj_gene = gene_from_call(row.get("j_call_heavy", ""))

            is_vitt_ighv = hv_gene in VITT_IGHV
            is_ighj3 = hj_gene == "IGHJ3"
            is_acidic = acidic_hcdr3(row)

            if k31e:
                stats["k31e_proxy"] += 1

            if is_ighj3:
                stats["ighj3"] += 1

            if is_acidic:
                stats["acidic_hcdr3"] += 1

            if is_vitt_ighv:
                stats["vitt_ighv"] += 1

            if is_vitt_ighv and is_ighj3:
                stats["vitt_ighv_ighj3"] += 1

            if (
                k31e
                and is_vitt_ighv
                and is_ighj3
                and is_acidic
            ):
                stats["complete_five_feature"] += 1

    return stats


def pct(n, d):
    return 100.0 * n / d if d else 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_dir",
        help="Directory containing paired OAS .csv.gz files",
    )
    parser.add_argument(
        "--output",
        default="paired_iglv321_results.tsv",
        help="Output TSV file",
    )
    args = parser.parse_args()

    files = sorted(
        glob.glob(os.path.join(args.input_dir, "*.csv.gz"))
    )

    if not files:
        raise RuntimeError(
            f"No .csv.gz files found in {args.input_dir}"
        )

    total = Counter()
    per_file = []

    for path in files:
        s = analyze_file(path)
        total.update(s)

        per_file.append({
            "file": os.path.basename(path),
            **s,
        })

    fields = [
        "file",
        "total_pairs",
        "iglv321",
        "k31e_all_iglv321",
        "k31_unmapped_iglv321",
        "dds_proxy",
        "k31e_proxy",
        "k31_unmapped_proxy",
        "ighj3",
        "acidic_hcdr3",
        "vitt_ighv",
        "vitt_ighv_ighj3",
        "complete_five_feature",
    ]

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()

        for row in per_file:
            writer.writerow(
                {k: row.get(k, 0) for k in fields}
            )

        writer.writerow({
            "file": "TOTAL",
            **{k: total.get(k, 0) for k in fields if k != "file"},
        })

    print("\nPAIRED ANALYSIS")
    print("-" * 72)

    print(f"Files analyzed:                    {len(files):,}")
    print(f"Total paired BCRs:                 {total['total_pairs']:,}")
    print(f"IGLV3-21 pairs:                    {total['iglv321']:,}")

    print(
        f"K31E among all IGLV3-21:           "
        f"{total['k31e_all_iglv321']:,} "
        f"({pct(total['k31e_all_iglv321'], total['iglv321']):.3f}%)"
    )

    print(
        f"K31 unmapped among IGLV3-21:       "
        f"{total['k31_unmapped_iglv321']:,}"
    )

    print(f"DDS(D) proxy pairs:                {total['dds_proxy']:,}")

    print(
        f"With K31E:                         "
        f"{total['k31e_proxy']:,} "
        f"({pct(total['k31e_proxy'], total['dds_proxy']):.2f}%)"
    )

    print(
        f"K31 unmapped among proxy pairs:    "
        f"{total['k31_unmapped_proxy']:,}"
    )

    print(
        f"With IGHJ3:                        "
        f"{total['ighj3']:,} "
        f"({pct(total['ighj3'], total['dds_proxy']):.1f}%)"
    )

    print(
        f"With acidic HCDR3 >=3 D/E:         "
        f"{total['acidic_hcdr3']:,} "
        f"({pct(total['acidic_hcdr3'], total['dds_proxy']):.1f}%)"
    )

    print(
        f"With VITT-reported IGHV:           "
        f"{total['vitt_ighv']:,} "
        f"({pct(total['vitt_ighv'], total['dds_proxy']):.1f}%)"
    )

    print(
        f"With VITT IGHV + IGHJ3:            "
        f"{total['vitt_ighv_ighj3']:,} "
        f"({pct(total['vitt_ighv_ighj3'], total['dds_proxy']):.1f}%)"
    )

    print(
        f"Complete five-feature combination: "
        f"{total['complete_five_feature']:,}"
    )

    print(f"\nResults written to: {args.output}")


if __name__ == "__main__":
    main()

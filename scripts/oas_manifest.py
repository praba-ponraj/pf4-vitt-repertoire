#!/usr/bin/env python3
"""Create a provenance manifest from bulk or paired OAS CSV archives."""

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path


FIELDS = [
    "Run", "Author", "Link", "Species", "Subject", "Age", "BSource",
    "BType", "Disease", "Vaccine", "Longitudinal", "Isotype", "Chain",
    "Unique sequences",
]


def parse_metadata_line(line, path):
    """Parse the JSON metadata record that precedes an OAS CSV header."""
    text = line.strip()
    try:
        metadata = json.loads(text)
    except json.JSONDecodeError:
        # Some OAS downloads represent the JSON record as one CSV-quoted field.
        try:
            fields = next(csv.reader([text]))
            if len(fields) != 1:
                raise ValueError("metadata record has multiple CSV fields")
            metadata = json.loads(fields[0])
        except (csv.Error, json.JSONDecodeError, StopIteration, ValueError) as exc:
            raise ValueError(f"invalid OAS metadata record: {exc}") from exc

    if not isinstance(metadata, dict):
        raise ValueError(f"{path}: OAS metadata record is not a JSON object")
    return metadata


def discover_files(input_dirs):
    """Return unique OAS archives found recursively below the input directories."""
    files = {}
    for directory in input_dirs:
        if not directory.is_dir():
            print(f"ERROR: input directory does not exist: {directory}", file=sys.stderr)
            continue
        for path in directory.rglob("*.csv.gz"):
            if path.is_file():
                files[path.resolve()] = path
    return sorted(files.values(), key=lambda path: str(path))


def read_rows(paths):
    """Read rows and report unavailable metadata fields to stderr."""
    rows = []
    errors = 0
    seen_names = set()
    for path in paths:
        try:
            with gzip.open(path, "rt", newline="") as handle:
                first_line = next(handle)
            metadata = parse_metadata_line(first_line, path)
        except (OSError, EOFError, StopIteration, UnicodeError, ValueError) as exc:
            print(f"ERROR: {path}: {exc}", file=sys.stderr)
            errors += 1
            continue

        missing = [field for field in FIELDS if metadata.get(field) in (None, "")]
        if missing:
            print(
                f"WARNING: {path.name}: blank or unavailable metadata fields: "
                + ", ".join(missing),
                file=sys.stderr,
            )
        if path.name in seen_names:
            print(
                f"WARNING: duplicate filename in manifest inputs: {path.name}",
                file=sys.stderr,
            )
        seen_names.add(path.name)
        rows.append([path.name] + [metadata.get(field, "") for field in FIELDS])
    return rows, errors


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Write a TSV provenance manifest for bulk or paired OAS .csv.gz files. "
            "TSV output is written to stdout; diagnostics are written to stderr."
        )
    )
    parser.add_argument(
        "input_dirs", metavar="DIR", nargs="+", type=Path,
        help="input directory to search recursively (one or more)",
    )
    args = parser.parse_args()

    paths = discover_files(args.input_dirs)
    if not paths:
        print("ERROR: no .csv.gz files found in the supplied directories", file=sys.stderr)
        return 1

    rows, errors = read_rows(paths)
    if errors:
        print(
            f"ERROR: manifest not written because {errors} file(s) could not be read",
            file=sys.stderr,
        )
        return 1

    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(["file"] + FIELDS)
    writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

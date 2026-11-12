"""Rank one cotton GWAS locus from a curated evidence CSV."""

import argparse
import csv
import math
from pathlib import Path


DEFAULT = dict(lof=10, missense=10, function=10, expression=10, homology=5,
               literature=5, mapping=30, regulation=20)
OPTIONAL = {"mapping", "regulation"}
OUTPUT = ["rank", "priority_score", "observed_evidence_score",
          "evidence_completeness", "strong_evidence_classes", "evidence_tier"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", "-o", type=Path, required=True)
    parser.add_argument("--weights", help="Comma-separated overrides, e.g. lof=30,missense=25")
    args = parser.parse_args()
    weights = DEFAULT.copy()
    try:
        overrides = set()
        for item in args.weights.split(",") if args.weights else []:
            name, raw = item.split("=", 1)
            name = name.strip()
            value = float(raw)
            if name not in weights or name in overrides or not math.isfinite(value) or value < 0:
                raise ValueError(f"Invalid weight: {item}")
            weights[name] = value
            overrides.add(name)
        total = sum(weights.values())
        if not math.isfinite(total) or total <= 0:
            raise ValueError("Total weight must be finite and positive.")
        weights = {name: value / total for name, value in weights.items()}
        rows = []
        seen = set()
        with args.input.open(encoding="utf-8-sig", newline="") as source:
            reader = csv.DictReader(source, strict=True)
            fields = reader.fieldnames or []
            required = {"gene_id", *(f"{name}_score" for name in weights.keys() - OPTIONAL)}
            if required - set(fields) or len(set(fields)) != len(fields) or set(fields) & set(OUTPUT):
                raise ValueError("Input needs gene_id and six unique score columns without reserved output fields.")
            for row in reader:
                gene = (row.get("gene_id") or "").strip()
                if None in row or None in row.values() or not gene or gene in seen:
                    raise ValueError(f"Malformed row or missing/duplicate gene ID on line {reader.line_num}.")
                seen.add(gene)
                points, coverage, strong = [], [], 0
                for name, weight in weights.items():
                    raw = row.get(f"{name}_score", "").strip()
                    if not raw:
                        continue
                    score = float(raw)
                    if not 0 <= score <= 1:
                        raise ValueError(f"Line {reader.line_num}: {name}_score must be in [0, 1].")
                    points.append(weight * score)
                    coverage.append(weight)
                    strong += weight > 0 and score >= .75
                rows.append((math.fsum(points), math.fsum(coverage), gene, strong, row))
        if not rows:
            raise ValueError("Input has no candidate rows.")
        if args.input.resolve() == args.output.resolve() or (
            args.output.exists() and args.input.samefile(args.output)
        ):
            raise ValueError("Input and output must be different files.")
        rows.sort(key=lambda record: (-record[0], -record[1], record[2]))
        with args.output.open("w", encoding="utf-8", newline="") as destination:
            writer = csv.DictWriter(destination, fieldnames=OUTPUT + fields)
            writer.writeheader()
            for rank, (points, coverage, _, strong, row) in enumerate(rows, 1):
                if coverage < .50:
                    tier = "Insufficient evidence"
                elif points >= .75 and coverage >= .75 and strong >= 2:
                    tier = "Tier 1"
                elif points >= .50:
                    tier = "Tier 2"
                else:
                    tier = "Tier 3"
                writer.writerow(dict(
                    row, rank=rank, priority_score=f"{points * 100:.2f}",
                    observed_evidence_score=f"{points / coverage * 100:.2f}" if coverage else "",
                    evidence_completeness=f"{coverage:.2f}",
                    strong_evidence_classes=strong, evidence_tier=tier,
                ))
    except (OSError, ValueError, csv.Error) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
repeat_counter.py

Genome-wide repeat read counter for BAM files.

What it does:
- scans BAM records directly
- includes unmapped reads if they have sequence
- removes PCR duplicates by default
- use -d / --include-duplicates to keep duplicates
- accepts multiple motifs in one run - specify these in your command line 
- change length of tract detected by altering min-bp
- automatically pools each motif with its reverse complement
- writes separate output files per motif
- keeps a simple ChIP vs Input comparison per motif

Example of what to put in command line to detect a selection of trinucleotides: 

  python3 repeat_counter.py \
    --chip chip.bam \
    --input input.bam \
    --motif AAT AAC AAG CAG CGG AGG \
    --min-bp 18 \
    --outdir /mnt/d/your/folder/here \
    --prefix Cell_line_condition_

"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import pysam


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def sanitize_token(text: str) -> str:
    text = text.upper().strip()
    text = re.sub(r"[^A-Z0-9_+-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "MOTIF"


def count_repeat_tracts(seq: str, motif: str, min_units: int) -> Tuple[bool, int]:
    """
    Return:
      - whether seq contains at least one tract of motif repeated >= min_units
      - total repeat units across all qualifying tracts in that orientation

    Designed for speed on read-length strings.
    """
    seq = seq.upper()
    motif = motif.upper()
    n = len(motif)

    pos = 0
    positive = False
    total_units = 0

    while True:
        idx = seq.find(motif, pos)
        if idx < 0:
            break

        # Skip internal shifted hits inside an already-counted tract.
        if idx >= n and seq[idx - n : idx] == motif:
            pos = idx + 1
            continue

        j = idx + n
        units = 1
        while seq.startswith(motif, j):
            units += 1
            j += n

        if units >= min_units:
            positive = True
            total_units += units

        pos = idx + 1

    return positive, total_units


def make_chrom_stats() -> dict:
    return {
        "raw_reads_seen": 0,
        "eligible_reads": 0,
        "mapped_reads": 0,
        "unmapped_reads": 0,
        "repeat_positive_reads": 0,
        "motif_positive_reads": 0,
        "rc_positive_reads": 0,
        "repeat_units": 0,
        "no_sequence_reads": 0,
    }


def build_motif_specs(motifs: List[str]) -> List[dict]:
    specs = []
    seen = set()
    for motif in motifs:
        m = sanitize_token(motif)
        rc = sanitize_token(revcomp(motif))
        # Keep user order, but avoid exact duplicate motif entries.
        key = m
        if key in seen:
            continue
        seen.add(key)
        specs.append(
            {
                "motif": m,
                "rc": rc,
                "same_as_rc": m == rc,
            }
        )
    return specs


def process_bam(
    bam_path: str,
    motif_specs: List[dict],
    min_bp: int,
    include_duplicates: bool,
    threads: int = 4,
    progress_every: int = 5_000_000,
) -> dict:
    # Precompute minimum repeat units per motif.
    for spec in motif_specs:
        spec["min_units"] = math.ceil(min_bp / len(spec["motif"]))

    by_motif = {
        spec["motif"]: defaultdict(make_chrom_stats)
        for spec in motif_specs
    }

    raw_reads_seen = 0
    duplicate_reads_seen = 0
    duplicate_reads_skipped = 0
    secondary_reads_skipped = 0
    supplementary_reads_skipped = 0

    with pysam.AlignmentFile(bam_path, "rb", threads=threads) as bam:
        for read in bam.fetch(until_eof=True):
            raw_reads_seen += 1

            if progress_every and raw_reads_seen % progress_every == 0:
                print(
                    f"[{Path(bam_path).name}] processed {raw_reads_seen:,} BAM records...",
                    file=sys.stderr,
                    flush=True,
                )

            if read.is_duplicate:
                duplicate_reads_seen += 1
                if not include_duplicates:
                    duplicate_reads_skipped += 1
                    continue

            if read.is_secondary:
                secondary_reads_skipped += 1
                continue

            if read.is_supplementary:
                supplementary_reads_skipped += 1
                continue

            seq = read.query_sequence
            if not seq:
                chrom = read.reference_name if read.reference_name else "unmapped"
                for spec in motif_specs:
                    stats = by_motif[spec["motif"]][chrom]
                    stats["raw_reads_seen"] += 1
                    stats["eligible_reads"] += 1
                    stats["no_sequence_reads"] += 1
                continue

            seq = seq.upper()
            chrom = "unmapped" if read.is_unmapped or not read.reference_name else read.reference_name

            for spec in motif_specs:
                motif = spec["motif"]
                rc = spec["rc"]
                min_units = spec["min_units"]
                same_as_rc = spec["same_as_rc"]

                stats = by_motif[motif][chrom]
                stats["raw_reads_seen"] += 1
                stats["eligible_reads"] += 1

                if read.is_unmapped:
                    stats["unmapped_reads"] += 1
                else:
                    stats["mapped_reads"] += 1

                motif_pos, motif_units = count_repeat_tracts(seq, motif, min_units)
                if same_as_rc:
                    rc_pos, rc_units = False, 0
                else:
                    rc_pos, rc_units = count_repeat_tracts(seq, rc, min_units)

                if motif_pos or rc_pos:
                    stats["repeat_positive_reads"] += 1
                    stats["motif_positive_reads"] += 1 if motif_pos else 0
                    stats["rc_positive_reads"] += 1 if rc_pos else 0
                    stats["repeat_units"] += motif_units + rc_units

    result = {
        "bam_path": bam_path,
        "min_bp": min_bp,
        "include_duplicates": include_duplicates,
        "raw_reads_seen": raw_reads_seen,
        "duplicate_reads_seen": duplicate_reads_seen,
        "duplicate_reads_skipped": duplicate_reads_skipped,
        "secondary_reads_skipped": secondary_reads_skipped,
        "supplementary_reads_skipped": supplementary_reads_skipped,
        "by_motif": {},
    }

    for spec in motif_specs:
        motif = spec["motif"]
        motif_stats = by_motif[motif]
        eligible = sum(d["eligible_reads"] for d in motif_stats.values()) or 1
        repeat_positive_reads = sum(d["repeat_positive_reads"] for d in motif_stats.values())
        repeat_units = sum(d["repeat_units"] for d in motif_stats.values())
        result["by_motif"][motif] = {
            "motif": motif,
            "reverse_complement": spec["rc"],
            "same_as_rc": spec["same_as_rc"],
            "min_units": spec["min_units"],
            "eligible_reads": eligible,
            "repeat_positive_reads": repeat_positive_reads,
            "repeat_units": repeat_units,
            "by_chrom": motif_stats,
        }

    return result


def chrom_sort_key(chrom: str):
    c = chrom
    if c.lower().startswith("chr"):
        c = c[3:]
    if c.isdigit():
        return (0, int(c))
    order = {"x": 23, "y": 24, "m": 25, "mt": 25, "unmapped": 1000}
    return (1, order.get(c.lower(), 2000), chrom)


def write_summary(out_path: Path, result: dict):
    eligible = result["eligible_reads"] or 1
    repeat_reads = result["repeat_positive_reads"]
    repeat_units = result["repeat_units"]

    with out_path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric", "value"])
        rows = [
            ("bam_path", result["bam_path"]),
            ("motif", result["motif"]),
            ("reverse_complement", result["reverse_complement"]),
            ("same_as_rc", result["same_as_rc"]),
            ("min_bp", result["min_bp"]),
            ("min_repeat_units", result["min_units"]),
            ("include_duplicates", result["include_duplicates"]),
            ("raw_reads_seen", result["raw_reads_seen"]),
            ("eligible_reads", result["eligible_reads"]),
            ("duplicate_reads_seen", result["duplicate_reads_seen"]),
            ("duplicate_reads_skipped", result["duplicate_reads_skipped"]),
            ("secondary_reads_skipped", result["secondary_reads_skipped"]),
            ("supplementary_reads_skipped", result["supplementary_reads_skipped"]),
            ("repeat_positive_reads", repeat_reads),
            ("repeat_positive_rpm", f"{repeat_reads / eligible * 1_000_000:.6f}"),
            ("repeat_units", repeat_units),
            ("repeat_units_per_million_reads", f"{repeat_units / eligible * 1_000_000:.6f}"),
            ("positive_fraction", f"{repeat_reads / eligible:.8f}"),
        ]
        for row in rows:
            w.writerow(row)


def write_by_chrom(out_path: Path, result: dict):
    by_chrom = result["by_chrom"]
    eligible = result["eligible_reads"] or 1

    with out_path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "chromosome",
            "raw_reads_seen",
            "eligible_reads",
            "mapped_reads",
            "unmapped_reads",
            "repeat_positive_reads",
            "repeat_positive_rpm",
            "repeat_units",
            "repeat_units_per_million_reads",
            "motif_positive_reads",
            "rc_positive_reads",
            "no_sequence_reads",
        ])
        for chrom in sorted(by_chrom.keys(), key=chrom_sort_key):
            d = by_chrom[chrom]
            w.writerow([
                chrom,
                d["raw_reads_seen"],
                d["eligible_reads"],
                d["mapped_reads"],
                d["unmapped_reads"],
                d["repeat_positive_reads"],
                f"{d['repeat_positive_reads'] / eligible * 1_000_000:.6f}",
                d["repeat_units"],
                f"{d['repeat_units'] / eligible * 1_000_000:.6f}",
                d["motif_positive_reads"],
                d["rc_positive_reads"],
                d["no_sequence_reads"],
            ])


def write_comparison(out_path: Path, chip: dict, inp: dict):
    chip_eligible = chip["eligible_reads"] or 1
    inp_eligible = inp["eligible_reads"] or 1

    chip_rpm = chip["repeat_positive_reads"] / chip_eligible * 1_000_000
    inp_rpm = inp["repeat_positive_reads"] / inp_eligible * 1_000_000
    log2fc = math.log2((chip_rpm + 1e-9) / (inp_rpm + 1e-9))

    with out_path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric", "chip", "input", "value"])
        w.writerow(["repeat_positive_reads", chip["repeat_positive_reads"], inp["repeat_positive_reads"], ""])
        w.writerow(["repeat_positive_rpm", f"{chip_rpm:.6f}", f"{inp_rpm:.6f}", ""])
        w.writerow(["log2_fold_change_rpm", "", "", f"{log2fc:.6f}"])


def make_output_paths(outdir: Path, prefix: str, motif: str) -> dict:
    base = outdir / f"{prefix}_{motif}"
    return {
        "chip_summary": base.with_name(base.name + "_chip_summary.tsv"),
        "chip_by_chrom": base.with_name(base.name + "_chip_by_chrom.tsv"),
        "input_summary": base.with_name(base.name + "_input_summary.tsv"),
        "input_by_chrom": base.with_name(base.name + "_input_by_chrom.tsv"),
        "comparison": base.with_name(base.name + "_comparison.tsv"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chip", required=True, help="ChIP BAM")
    ap.add_argument("--input", required=True, help="Input BAM")
    ap.add_argument("--motif", nargs="+", required=True,
                    help="One or more repeat motifs, e.g. AAT AAC AAG CAG CGG AGG")
    ap.add_argument("--min-bp", type=int, default=18, help="Minimum consecutive repeat bp")
    ap.add_argument("--outdir", default=".", help="Output directory")
    ap.add_argument("--prefix", required=True, help="Output file prefix")
    ap.add_argument("-d", "--include-duplicates", action="store_true",
                    help="Include PCR duplicates instead of removing them")
    ap.add_argument("--threads", type=int, default=4, help="HTSlib decompression threads")
    ap.add_argument("--progress-every", type=int, default=5_000_000,
                    help="Print progress every N BAM records (0 disables)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    motifs = [sanitize_token(m) for m in args.motif]
    motif_specs = build_motif_specs(motifs)

    chip_all = process_bam(
        bam_path=args.chip,
        motif_specs=motif_specs,
        min_bp=args.min_bp,
        include_duplicates=args.include_duplicates,
        threads=args.threads,
        progress_every=args.progress_every,
    )

    inp_all = process_bam(
        bam_path=args.input,
        motif_specs=motif_specs,
        min_bp=args.min_bp,
        include_duplicates=args.include_duplicates,
        threads=args.threads,
        progress_every=args.progress_every,
    )

    for spec in motif_specs:
        motif = spec["motif"]
        paths = make_output_paths(outdir, args.prefix, motif)

        chip = chip_all["by_motif"][motif].copy()
        chip.update({
            "bam_path": chip_all["bam_path"],
            "min_bp": chip_all["min_bp"],
            "include_duplicates": chip_all["include_duplicates"],
            "raw_reads_seen": chip_all["raw_reads_seen"],
            "duplicate_reads_seen": chip_all["duplicate_reads_seen"],
            "duplicate_reads_skipped": chip_all["duplicate_reads_skipped"],
            "secondary_reads_skipped": chip_all["secondary_reads_skipped"],
            "supplementary_reads_skipped": chip_all["supplementary_reads_skipped"],
        })

        inp = inp_all["by_motif"][motif].copy()
        inp.update({
            "bam_path": inp_all["bam_path"],
            "min_bp": inp_all["min_bp"],
            "include_duplicates": inp_all["include_duplicates"],
            "raw_reads_seen": inp_all["raw_reads_seen"],
            "duplicate_reads_seen": inp_all["duplicate_reads_seen"],
            "duplicate_reads_skipped": inp_all["duplicate_reads_skipped"],
            "secondary_reads_skipped": inp_all["secondary_reads_skipped"],
            "supplementary_reads_skipped": inp_all["supplementary_reads_skipped"],
        })

        write_summary(paths["chip_summary"], chip)
        write_by_chrom(paths["chip_by_chrom"], chip)

        write_summary(paths["input_summary"], inp)
        write_by_chrom(paths["input_by_chrom"], inp)

        write_comparison(paths["comparison"], chip, inp)

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
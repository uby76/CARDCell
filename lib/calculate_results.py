#!/usr/bin/env python3
"""Calculate KMA-Depth/CARD ARG copies per 16S-derived bacterial cell."""
import argparse
import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

csv.field_size_limit(sys.maxsize)
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
REFERENCE_16S_LENGTH = 1432.0


def open_fastq(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def fastq_stats(path):
    reads = bases = 0
    with open_fastq(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().rstrip("\r\n")
            plus = handle.readline()
            qual = handle.readline().rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+") or len(seq) != len(qual):
                raise ValueError(f"Invalid FASTQ structure near read {reads + 1}: {path}")
            reads += 1
            bases += len(seq)
    if reads == 0:
        raise ValueError(f"Empty FASTQ: {path}")
    return reads, bases, bases / reads


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_lineages(paths):
    lineages, mate_counts = [], []
    for path in paths:
        count = 0
        with open(path) as handle:
            for line in handle:
                fields = line.rstrip().split("\t")
                if len(fields) < 2:
                    continue
                taxa = [x.strip() for x in fields[1].split(";") if x.strip()]
                if taxa and taxa[0].casefold() == "bacteria":
                    lineages.append(taxa[:7]); count += 1
        mate_counts.append(count)
    return lineages, mate_counts


def rrndb_index(path):
    result = {}
    for row in read_tsv(path):
        if row.get("mean"):
            result[(row["rank"].lower(), row["name"].strip().casefold())] = row
    return result


def group_name(lineage):
    if len(lineage) >= 6:
        return "genus", lineage[5]
    idx = len(lineage) - 1
    return RANKS[idx], f"Unclassified_{lineage[idx]}"


def match_rrndb(lineage, db):
    for idx in range(min(5, len(lineage) - 1), -1, -1):
        row = db.get((RANKS[idx], lineage[idx].casefold()))
        if row:
            return lineage[idx], RANKS[idx], float(row["mean"])
    return "", "", None


def labels(text):
    return [x.strip() for x in text.split(";") if x.strip()] or ["Unspecified"]


def qc_flags(mapped_reads, percent_coverage, low_reads, low_coverage):
    flags = []
    if mapped_reads < low_reads:
        flags.append(f"LOW_MAPPED_READS_lt_{low_reads}")
    if percent_coverage < low_coverage:
        flags.append(f"LOW_PERCENT_COVERAGE_lt_{low_coverage:g}")
    return ";".join(flags) if flags else "PASS"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", required=True); ap.add_argument("--r2", required=True)
    ap.add_argument("--sample", required=True); ap.add_argument("--sample-dir", required=True)
    ap.add_argument("--rrndb", required=True)
    ap.add_argument("--qc-low-mapped-reads", type=int, default=3)
    ap.add_argument("--qc-low-percent-coverage", type=float, default=10.0)
    ap.add_argument("--qc-low-16s-reads", type=int, default=100)
    ap.add_argument("--comparison-ko30", type=float)
    ap.add_argument("--comparison-previous-16s", type=float)
    ap.add_argument("--comparison-previous-16s-count", type=int)
    args = ap.parse_args()
    out = Path(args.sample_dir).resolve(); out.mkdir(parents=True, exist_ok=True)

    r1_reads, r1_bases, r1_mean = fastq_stats(args.r1)
    r2_reads, r2_bases, r2_mean = fastq_stats(args.r2)
    if r1_reads != r2_reads:
        raise ValueError(f"Paired read-count mismatch: R1={r1_reads}, R2={r2_reads}")
    total_reads = r1_reads + r2_reads
    combined_mean = (r1_bases + r2_bases) / total_reads
    write_tsv(out / "01_read_statistics.tsv", ["sample", "R1_reads", "R2_reads", "total_reads", "R1_mean_length", "R2_mean_length", "combined_mean_length"], [{
        "sample": args.sample, "R1_reads": r1_reads, "R2_reads": r2_reads, "total_reads": total_reads,
        "R1_mean_length": f"{r1_mean:.12f}", "R2_mean_length": f"{r2_mean:.12f}", "combined_mean_length": f"{combined_mean:.12f}"}])

    prefix = out / "rgi_bwt" / args.sample
    card_rows = read_tsv(Path(str(prefix) + ".allele_mapping_data.txt"))
    depth_rows, qc_depth_rows = [], []
    for row in card_rows:
        raw = {"sample": args.sample, "Reference_Sequence": row["Reference Sequence"], "ARO_Term": row["ARO Term"],
            "ARO_Accession": row["ARO Accession"], "Mapped_Reads": row["All Mapped Reads"],
            "Reference_Model_Type": row["Reference Model Type"], "Reference_DB": row["Reference DB"],
            "Reference_Allele_Source": row["Reference Allele Source"],
            "Completely_Mapped_Reads": row["Completely Mapped Reads"],
            "Mapped_Reads_with_Flanking_Sequence": row["Mapped Reads with Flanking Sequence"],
            "Percent_Coverage": row["Percent Coverage"], "Length_Coverage_bp": row["Length Coverage (bp)"],
            "Reference_Length": row["Reference Length"], "KMA_Depth": f"{float(row['Depth']):.12f}",
            "AMR_Gene_Family": row["AMR Gene Family"], "Drug_Class": row["Drug Class"],
            "Resistance_Mechanism": row["Resistance Mechanism"]}
        depth_rows.append(raw)
        flagged = dict(raw)
        flagged["QC_Flag"] = qc_flags(int(float(row["All Mapped Reads"])), float(row["Percent Coverage"]), args.qc_low_mapped_reads, args.qc_low_percent_coverage)
        flagged["QC_Action"] = "retain_review_recommended" if flagged["QC_Flag"] != "PASS" else "retain"
        qc_depth_rows.append(flagged)
    depth_fields = ["sample", "Reference_Sequence", "ARO_Term", "ARO_Accession", "Reference_Model_Type",
        "Reference_DB", "Reference_Allele_Source", "Completely_Mapped_Reads",
        "Mapped_Reads_with_Flanking_Sequence", "Mapped_Reads", "Percent_Coverage",
        "Length_Coverage_bp", "Reference_Length", "KMA_Depth", "AMR_Gene_Family", "Drug_Class",
        "Resistance_Mechanism"]
    write_tsv(out / "02_CARD_KMA_depth.tsv", depth_fields, depth_rows)
    write_tsv(out / "02_CARD_KMA_depth_QC.tsv", depth_fields + ["QC_Flag", "QC_Action"], qc_depth_rows)

    tax_files = [out / "metaxa2" / f"{args.sample}_R1.taxonomy.txt", out / "metaxa2" / f"{args.sample}_R2.taxonomy.txt"]
    lineages, mate_16s = read_lineages(tax_files)
    n16s = len(lineages)
    if n16s == 0:
        raise ValueError("No bacterial 16S reads detected; cell denominator cannot be calculated")
    coverage16s = n16s * combined_mean / REFERENCE_16S_LENGTH
    write_tsv(out / "03_Metaxa2_16S.tsv", ["sample", "N_16S_R1", "N_16S_R2", "N_16S_total", "combined_mean_read_length", "16S_reference_length", "16S_coverage"], [{
        "sample": args.sample, "N_16S_R1": mate_16s[0], "N_16S_R2": mate_16s[1], "N_16S_total": n16s,
        "combined_mean_read_length": f"{combined_mean:.12f}", "16S_reference_length": int(REFERENCE_16S_LENGTH), "16S_coverage": f"{coverage16s:.12f}"}])

    db = rrndb_index(args.rrndb)
    groups = {}
    for lineage in lineages:
        key = group_name(lineage)
        if key not in groups: groups[key] = {"count": 0, "lineages": Counter()}
        groups[key]["count"] += 1; groups[key]["lineages"][tuple(lineage)] += 1
    rrn_rows, matched_reads, weighted_raw = [], 0, 0.0
    for (rank, taxon), info in groups.items():
        rel = info["count"] / n16s
        lineage = list(info["lineages"].most_common(1)[0][0])
        match_name, match_rank, copy_number = match_rrndb(lineage, db)
        matched = copy_number is not None
        contribution = rel * copy_number if matched else 0.0
        if matched: matched_reads += info["count"]; weighted_raw += contribution
        rrn_rows.append({"Taxon": taxon, "Taxonomic_Rank": rank, "Relative_Abundance": f"{rel:.12f}",
            "rrnDB_Matched_Name": match_name, "rrnDB_Matched_Rank": match_rank,
            "rrn_Copy_Number": "" if copy_number is None else f"{copy_number:.6f}",
            "Fallback_Level": "none" if rank == "genus" and match_rank == "genus" else (match_rank or "unmatched"),
            "Weighted_Contribution": f"{contribution:.12f}", "Match_Status": "matched" if matched else "unmatched"})
    rrn_rows.sort(key=lambda x: (-float(x["Relative_Abundance"]), x["Taxon"]))
    rrn_fields = ["Taxon", "Taxonomic_Rank", "Relative_Abundance", "rrnDB_Matched_Name", "rrnDB_Matched_Rank", "rrn_Copy_Number", "Fallback_Level", "Weighted_Contribution", "Match_Status"]
    write_tsv(out / "04_rrnDB_copy_number.tsv", rrn_fields, rrn_rows)
    matched_fraction = matched_reads / n16s
    if matched_fraction == 0: raise ValueError("No Metaxa2 taxon matched rrnDB")
    weighted_copy = weighted_raw / matched_fraction
    cell_coverage = coverage16s / weighted_copy
    write_tsv(out / "05_cell_coverage.tsv", ["sample", "N_16S_reads", "mean_read_length", "16S_reference_length", "16S_coverage", "weighted_mean_16S_copy_number", "cell_equivalent_coverage", "matched_taxonomic_abundance_fraction"], [{
        "sample": args.sample, "N_16S_reads": n16s, "mean_read_length": f"{combined_mean:.12f}", "16S_reference_length": int(REFERENCE_16S_LENGTH),
        "16S_coverage": f"{coverage16s:.12f}", "weighted_mean_16S_copy_number": f"{weighted_copy:.12f}",
        "cell_equivalent_coverage": f"{cell_coverage:.12f}", "matched_taxonomic_abundance_fraction": f"{matched_fraction:.12f}"}])

    copies_rows, qc_copies_rows = [], []
    for raw, qc in zip(depth_rows, qc_depth_rows):
        depth = float(raw["KMA_Depth"])
        result = {"sample": args.sample, "Reference_Sequence": raw["Reference_Sequence"], "ARO_Term": raw["ARO_Term"],
            "ARO_Accession": raw["ARO_Accession"], "KMA_Depth": raw["KMA_Depth"], "cell_equivalent_coverage": f"{cell_coverage:.12f}",
            "ARG_copies_per_cell": f"{depth/cell_coverage:.12f}", "Mapped_Reads": raw["Mapped_Reads"],
            "Percent_Coverage": raw["Percent_Coverage"], "AMR_Gene_Family": raw["AMR_Gene_Family"],
            "Drug_Class": raw["Drug_Class"], "Resistance_Mechanism": raw["Resistance_Mechanism"]}
        copies_rows.append(result)
        flagged = dict(result); flagged["QC_Flag"] = qc["QC_Flag"]; flagged["QC_Action"] = qc["QC_Action"]
        qc_copies_rows.append(flagged)
    copies_fields = ["sample", "Reference_Sequence", "ARO_Term", "ARO_Accession", "KMA_Depth", "cell_equivalent_coverage", "ARG_copies_per_cell", "Mapped_Reads", "Percent_Coverage", "AMR_Gene_Family", "Drug_Class", "Resistance_Mechanism"]
    write_tsv(out / "06_ARG_copies_per_cell.tsv", copies_fields, copies_rows)
    write_tsv(out / "06_ARG_copies_per_cell_QC.tsv", copies_fields + ["QC_Flag", "QC_Action"], qc_copies_rows)

    total_depth = sum(float(x["KMA_Depth"]) for x in copies_rows)
    total_copies = sum(float(x["ARG_copies_per_cell"]) for x in copies_rows)
    aro_count = len({x["ARO_Accession"] for x in copies_rows})
    sample_flag = f"LOW_16S_READS_lt_{args.qc_low_16s_reads}" if n16s < args.qc_low_16s_reads else "PASS"
    write_tsv(out / "07_sample_summary.tsv", ["sample", "total_reads", "mean_read_length", "N_16S_reads", "16S_coverage", "weighted_mean_16S_copy_number", "cell_equivalent_coverage", "number_CARD_references_detected", "number_ARO_terms_detected", "total_KMA_ARG_Depth", "total_ARG_copies_per_cell", "sample_16S_QC_flag"], [{
        "sample": args.sample, "total_reads": total_reads, "mean_read_length": f"{combined_mean:.12f}", "N_16S_reads": n16s,
        "16S_coverage": f"{coverage16s:.12f}", "weighted_mean_16S_copy_number": f"{weighted_copy:.12f}",
        "cell_equivalent_coverage": f"{cell_coverage:.12f}", "number_CARD_references_detected": len(copies_rows),
        "number_ARO_terms_detected": aro_count, "total_KMA_ARG_Depth": f"{total_depth:.12f}",
        "total_ARG_copies_per_cell": f"{total_copies:.12f}", "sample_16S_QC_flag": sample_flag}])

    summaries = []
    for level, field in [("ARO Term", "ARO_Term"), ("AMR Gene Family", "AMR_Gene_Family"), ("Drug Class", "Drug_Class"), ("Resistance Mechanism", "Resistance_Mechanism")]:
        agg = defaultdict(float)
        for row in copies_rows:
            item_labels = [row[field]] if level == "ARO Term" else labels(row[field])
            share = float(row["ARG_copies_per_cell"]) / len(item_labels)
            for label in item_labels: agg[label] += share
        for label, value in sorted(agg.items(), key=lambda x: (-x[1], x[0])):
            summaries.append({"Summary_Level": level, "Label": label, "ARG_copies_per_cell": f"{value:.12f}",
                "Aggregation_Rule": "sum reference values by ARO term" if level == "ARO Term" else "fractional allocation across multiple labels"})
    write_tsv(out / "07_abundance_by_annotation.tsv", ["Summary_Level", "Label", "ARG_copies_per_cell", "Aggregation_Rule"], summaries)

    validation = [
        {"Metric": "new_16S_derived_cell_equivalent_coverage", "Value": f"{cell_coverage:.12f}", "Use": "primary denominator"},
        {"Metric": "matched_abundance_fraction", "Value": f"{matched_fraction:.12f}", "Use": "rrnDB QC"},
        {"Metric": "unmatched_abundance_fraction", "Value": f"{1-matched_fraction:.12f}", "Use": "rrnDB QC"},
        {"Metric": "matched_abundances_renormalized", "Value": "yes" if matched_fraction < 1 else "no_not_needed", "Use": "method record"}]
    if args.comparison_previous_16s is not None:
        validation += [{"Metric": "previous_16S_derived_cell_equivalent_coverage", "Value": f"{args.comparison_previous_16s:.12f}", "Use": "comparison only"},
            {"Metric": "new_minus_previous_16S_difference", "Value": f"{cell_coverage-args.comparison_previous_16s:.12f}", "Use": "comparison only"}]
    if args.comparison_previous_16s_count is not None:
        validation.append({"Metric": "agreement_with_previous_16S_read_count", "Value": "yes" if n16s == args.comparison_previous_16s_count else "no", "Use": "comparison only"})
    if args.comparison_ko30 is not None:
        validation += [{"Metric": "ARGsOAP_KO30_nCell", "Value": f"{args.comparison_ko30:.15f}", "Use": "comparison only; never denominator"},
            {"Metric": "16S_to_KO30_ratio", "Value": f"{cell_coverage/args.comparison_ko30:.12f}", "Use": "comparison only"}]
    write_tsv(out / "validation_comparison.tsv", ["Metric", "Value", "Use"], validation)

    print(f"Sample: {args.sample}"); print(f"Total reads: {total_reads}"); print(f"Combined mean read length: {combined_mean:.12f}")
    print(f"Bacterial 16S reads: {n16s}"); print(f"16S coverage: {coverage16s:.12f}")
    print(f"Weighted mean 16S copies/cell: {weighted_copy:.12f}"); print(f"Cell-equivalent coverage: {cell_coverage:.12f}")
    print(f"CARD references detected: {len(copies_rows)}"); print(f"ARO terms detected: {aro_count}")
    print(f"Total KMA ARG Depth: {total_depth:.12f}"); print(f"Total ARG copies/cell: {total_copies:.12f}")
    print("ARG copies/cell formula: KMA Depth / 16S-derived cell-equivalent coverage")
    if args.comparison_previous_16s is not None:
        print(f"Previous 16S-derived cell-equivalent coverage: {args.comparison_previous_16s:.12f}")
        print(f"Difference: {cell_coverage-args.comparison_previous_16s:.12f}")
    if args.comparison_ko30 is not None:
        print(f"Comparison-only ARGs-OAP KO30 nCell: {args.comparison_ko30:.15f}")
        print(f"16S/KO30 ratio: {cell_coverage/args.comparison_ko30:.12f}")


if __name__ == "__main__": main()

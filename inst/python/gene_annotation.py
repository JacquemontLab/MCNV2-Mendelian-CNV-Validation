#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import polars as pl
import subprocess
import tempfile
import shutil
import sys
import os

"""
CNV–Gene Intersection DB (Polars, robust)
- CNV input: first 5 columns are positionally mapped to Chr, Start, Stop, Type, SampleID (renamed exactly like this).
- All remaining CNV columns are kept as-is and appear after those 5 columns.
- gene_resource: 4 cols (chrom, start, end, attributes with GTF-like key="value"; ...).
- prob_regions: TSV/BED (with/without header); optional GenomeVersion filter; no Region filter.
- Output column order:
  [Chr, Start, Stop, Type, SampleID] + [other CNV columns] + [
     t_Start, t_End, transcript_length,
     gene_id, gene_name, gene_type, transcript_id, transcript_name,
     bp_overlap, cnv_problematic_region_overlap
  ]
"""

GENE_REGION_COLUMNS = [
    "t_Start", "t_End", "transcript_length",
    "gene_id", "gene_name", "gene_type", "transcript_id", "transcript_name",
    "bp_overlap", "cnv_problematic_region_overlap",
]

GENE_REGION_DTYPES = {
    "t_Start": pl.Int64,
    "t_End": pl.Int64,
    "transcript_length": pl.Int64,
    "gene_id": pl.Utf8,
    "gene_name": pl.Utf8,
    "gene_type": pl.Utf8,
    "transcript_id": pl.Utf8,
    "transcript_name": pl.Utf8,
    "bp_overlap": pl.Int64,
    "cnv_problematic_region_overlap": pl.Float64,
}

def _normalize_gv_arg(gv: str) -> str:
    s = (gv or "").strip().lower()
    if "grch38" in s or "hg38" in s or s == "38":
        return "38"
    if "grch37" in s or "hg19" in s or s in {"19", "37"}:
        return "37"
    digits = "".join(ch for ch in s if ch.isdigit())
    if digits == "38":
        return "38"
    if digits in {"37", "19"}:
        return "37"
    return "38"

def _pattern_for_gv(norm: str) -> str:
    return r"(grch38|hg38|38)" if norm == "38" else r"(grch37|hg19|37|19)"

def _load_cnv_positional(path: str) -> pl.LazyFrame:
    """
    Read CNV table (TSV). Regardless of header/names, the FIRST FIVE columns
    are mapped to Chr, Start, Stop, Type, SampleID. Remaining columns are kept.
    """
    try:
        lf = pl.scan_csv(path, separator="\t", has_header=True)
        cols = lf.collect_schema().names()
        if len(cols) < 5:
            raise ValueError
        c1, c2, c3, c4, c5 = cols[0], cols[1], cols[2], cols[3], cols[4]
    except Exception:
        lf = pl.scan_csv(path, separator="\t", has_header=False)
        cols = lf.collect_schema().names()
        if len(cols) < 5:
            raise ValueError("CNV file must have at least 5 columns (Chr, Start, Stop, Type, SampleID).")
        c1, c2, c3, c4, c5 = cols[0], cols[1], cols[2], cols[3], cols[4]

    other_cols = cols[5:]

    lf = (
        lf.select(
            [
                pl.col(c1).alias("Chr"),
                pl.col(c2).cast(pl.Int64, strict=False).alias("Start"),
                pl.col(c3).cast(pl.Int64, strict=False).alias("Stop"),
                pl.col(c4).cast(pl.Utf8,  strict=False).alias("Type"),
                pl.col(c5).cast(pl.Utf8,  strict=False).alias("SampleID"),
            ] + [pl.col(c) for c in other_cols]
        )
        .drop_nulls(["Start", "Stop"])
    )
    return lf

def _load_prob_regions_flex(path: str, genome_version: str) -> pl.LazyFrame:
    """
    Read prob_regions (TSV/BED) with/without header.
    - If header: try to find Chr/Start/End|Stop; if GenomeVersion exists, filter by alias.
    - If no header: first 3 columns are Chr, Start, Stop.
    """
    norm = _normalize_gv_arg(genome_version)
    pattern = _pattern_for_gv(norm)

    try:
        lf = pl.scan_csv(path, separator="\t", has_header=True)
        cols = lf.collect_schema().names()
        if len(cols) < 3:
            raise ValueError

        def pick(name_opts, fallback):
            for opt in name_opts:
                for c in cols:
                    if c.lower() == opt:
                        return c
            return fallback

        chr_col   = pick(["chr", "chrom", "chromosome", "contig"], cols[0])
        start_col = pick(["start", "pos", "begin", "position"],    cols[1])
        end_col   = pick(["stop", "end", "stoppos", "endpos"],     cols[2])

        gv_col = next((c for c in cols if c.lower() == "genomeversion"), None)
        if gv_col:
            lf = lf.filter(pl.col(gv_col).cast(pl.Utf8).str.to_lowercase().str.contains(pattern))
        else:
            print(f"[WARN] 'GenomeVersion' not found in {path}. No genome-version filtering.", file=sys.stderr)

        lf = (
            lf.select([
                pl.col(chr_col).alias("Chr"),
                pl.col(start_col).cast(pl.Int64, strict=False).alias("Start"),
                pl.col(end_col).cast(pl.Int64, strict=False).alias("Stop"),
            ])
            .drop_nulls(["Start", "Stop"])
        )
        return lf

    except Exception:
        lf = pl.scan_csv(path, separator="\t", has_header=False)
        return (
            lf.select([
                pl.col("column_1").alias("Chr"),
                pl.col("column_2").cast(pl.Int64, strict=False).alias("Start"),
                pl.col("column_3").cast(pl.Int64, strict=False).alias("Stop"),
            ])
            .drop_nulls(["Start", "Stop"])
        )

def _finalize_order(joined: pl.LazyFrame | pl.DataFrame, cnv_cols_after_norm: list[str]) -> pl.DataFrame:
    """
    Output order:
      1) Chr, Start, Stop, Type, SampleID
      2) other CNV columns (original CNV order), excluding the 5 fixed ones
      3) hard-coded gene/region columns (created if missing)
    """
    front = ["Chr", "Start", "Stop", "Type", "SampleID"]
    lf = joined.lazy() if isinstance(joined, pl.DataFrame) else joined
    existing = set(lf.collect_schema().names())

    cnv_extras = [c for c in cnv_cols_after_norm if c not in set(front)]

    add_exprs = []
    for col in GENE_REGION_COLUMNS:
        if col not in existing:
            add_exprs.append(pl.lit(None, dtype=GENE_REGION_DTYPES.get(col, pl.Utf8)).alias(col))
    if add_exprs:
        lf = lf.with_columns(add_exprs)

    selects = [pl.col(c) for c in front]
    selects += [pl.col(c) for c in cnv_extras if c in existing]
    selects += [
        pl.col(col).cast(GENE_REGION_DTYPES.get(col, pl.Utf8), strict=False).alias(col)
        for col in GENE_REGION_COLUMNS
    ]
    return lf.select(selects).collect()

def main():
    parser = argparse.ArgumentParser(
        description="Build CNV–gene intersection DB with fixed front columns (Chr, Start, Stop, Type, SampleID) and hard-coded gene columns."
    )
    parser.add_argument("--cnv", required=True, help="CNV TSV. First five columns are Chr, Start, Stop, Type, SampleID (positional).")
    parser.add_argument("--gene_resource", required=True, help="Gene resource: 4 columns [chrom, start, end, attributes] (GTF-like).")
    parser.add_argument("--prob_regions", required=True, help="Problematic regions TSV/BED (with/without header).")
    parser.add_argument("--genome_version", default="GRCh38", help="Genome alias among GRCh38, GRCh37, hg19, hg38, 19, 38.")
    parser.add_argument("--bedtools_path", default="bedtools", help="Absolute path to bedtools executable")
    parser.add_argument("--out", default="data/cnv_geneDB.tsv", help="Output TSV (default: data/cnv_geneDB.tsv)")
    args = parser.parse_args()

    # CNV
    cnvLF = _load_cnv_positional(args.cnv)  # Chr, Start, Stop, Type, SampleID + extras
    cnv_cols_after_norm = cnvLF.collect_schema().names()

    # prob_regions
    pRegions = _load_prob_regions_flex(args.prob_regions, args.genome_version)

    tmpdir = tempfile.mkdtemp()
    try:
        # Write 3-col BEDs for bedtools
        cnv_bed = os.path.join(tmpdir, "tmp_cnvs.bed")
        cnvLF.select(["Chr", "Start", "Stop"]).unique().sink_csv(cnv_bed, separator="\t", include_header=False)

        pRegions_bed = os.path.join(tmpdir, "tmp_pRegions.bed")
        pRegions.select(["Chr", "Start", "Stop"]).sink_csv(pRegions_bed, separator="\t", include_header=False)

        # bedtools: CNV × Genes
        inter_bed = os.path.join(tmpdir, "tmp_intersect.bed")
        with open(inter_bed, "w") as fout:
            subprocess.run(
                ["bedtools", "intersect", "-a", cnv_bed, "-b", args.gene_resource, "-F", "0.1", "-wao"],
                check=True, stdout=fout
            )

        # bedtools: CNV × prob_regions
        preg_inter_bed = os.path.join(tmpdir, "tmp_pRegion_intersect.bed")
        with open(preg_inter_bed, "w") as fout:
            subprocess.run(
                ["bedtools", "intersect", "-a", cnv_bed, "-b", pRegions_bed, "-wao"],
                check=True, stdout=fout
            )

        # Parse CNV × Genes (B has 4 columns: chrom, start, end, attributes)
        lf_inter = pl.scan_csv(inter_bed, separator="\t", has_header=False)
        intersected_genes = (
            lf_inter
              .select([
                  pl.col("column_1").alias("Chr"),
                  pl.col("column_2").cast(pl.Int64, strict=False).alias("Start"),
                  pl.col("column_3").cast(pl.Int64, strict=False).alias("Stop"),
                  pl.col("column_5").cast(pl.Int64, strict=False).alias("t_Start"),
                  pl.col("column_6").cast(pl.Int64, strict=False).alias("t_End"),
                  pl.col("column_7").cast(pl.Utf8,   strict=False).alias("attrs"),
                  pl.col("column_8").cast(pl.Int64,  strict=False).alias("bp_overlap"),
              ])
              # extract directly from 'attrs' without creating 'attrs_norm'
              .with_columns([
                  pl.col("attrs").fill_null("").cast(pl.Utf8).alias("attrs"),
              ])
              .with_columns([
                  pl.col("attrs").str.replace_all('""', '"').alias("attrs_clean"),
              ])
              .with_columns([
                  pl.col("attrs_clean").str.strip_chars('"').alias("attrs_clean"),
                  pl.col("attrs_clean").str.extract(r'gene_id\s+"([^"]+)"').alias("gene_id"),
                  pl.col("attrs_clean").str.extract(r'transcript_id\s+"([^"]+)"').alias("transcript_id"),
                  pl.col("attrs_clean").str.extract(r'gene_type\s+"([^"]+)"').alias("gene_type"),
                  pl.col("attrs_clean").str.extract(r'gene_name\s+"([^"]+)"').alias("gene_name"),
                  pl.col("attrs_clean").str.extract(r'transcript_name\s+"([^"]+)"').alias("transcript_name"),
                  (pl.col("t_End") - pl.col("t_Start")).alias("transcript_length"),
              ])
              # don't drop aggressively; we'll select final columns later
        )

        # Materialize to avoid lazy planning quirks
        intersected_genes_df = intersected_genes.collect()

        # Parse CNV × prob_regions; compute max fraction per CNV
        lf_preg = pl.scan_csv(preg_inter_bed, separator="\t", has_header=False)
        cnvs_with_pRegion = (
            lf_preg
              .select([
                  pl.col("column_1").alias("Chr"),
                  pl.col("column_2").cast(pl.Int64, strict=False).alias("Start"),
                  pl.col("column_3").cast(pl.Int64, strict=False).alias("Stop"),
                  pl.col("column_7").cast(pl.Int64, strict=False).alias("bp_overlap"),
              ])
              .with_columns([
                  (pl.col("bp_overlap") / (pl.col("Stop") - pl.col("Start") + 1.0))
                  .alias("cnv_problematic_region_overlap")
              ])
              .group_by(["Chr", "Start", "Stop"])
              .agg(pl.max("cnv_problematic_region_overlap").alias("cnv_problematic_region_overlap"))
        )
        cnvs_with_pRegion_df = cnvs_with_pRegion.collect()

        # Join & finalize order (use the materialized DataFrames converted back to Lazy)
        joined = (
            cnvLF
              .join(intersected_genes_df.lazy(), on=["Chr", "Start", "Stop"], how="left")
              .join(cnvs_with_pRegion_df.lazy(),  on=["Chr", "Start", "Stop"], how="left")
        )

        final_df = _finalize_order(joined, cnv_cols_after_norm)
        final_df.write_csv(args.out, separator="\t")

    finally:
        shutil.rmtree(tmpdir)

if __name__ == "__main__":
    main()

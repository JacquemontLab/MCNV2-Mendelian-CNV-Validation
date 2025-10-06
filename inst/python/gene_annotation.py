#!/usr/bin/env python3
import argparse
import polars as pl
import subprocess
import tempfile
import shutil
import sys

"""
CNV-Gene Intersection Database Builder (STOP au lieu de End, prob_regions partout, sans filtre Region)

- CNV : on prend toujours les 3 premières colonnes comme (Chr, Start, STOP), avec/sans en-tête.
- Toutes les autres colonnes CNV sont conservées telles quelles dans la sortie.
- prob_regions : on lit le fichier passé via --prob_regions.
  * en-tête ou non : on détecte Chr/Start/(End|STOP) ; sinon on prend les 3 premières colonnes.
  * si 'GenomeVersion' existe (casse ignorée), on filtre en acceptant les alias:
      --genome_version ∈ {GRCh38, GRCh37, hg19, hg38, 19, 38} (casse ignorée).
    Sinon: WARNING et pas de filtre (on utilise toutes les lignes).
  * on ignore toute colonne 'Region' (aucun filtre dessus).
- bedtools : intersections sur les BED 3 colonnes (Chr, Start, STOP) du run en tmpdir.
- pRegions : on calcule la fraction d'overlap, puis on prend le MAX par (Chr, Start, STOP).
- Join final sur (Chr, Start, STOP).
"""

def _ci_lookup(cols, wanted):
    w = wanted.lower()
    for c in cols:
        if c.lower() == w:
            return c
    return None

def _normalize_gv_arg(gv: str) -> str:
    """
    Normalise l'argument genome_version en un canon court:
    - '38' pour {GRCh38, hg38, 38}
    - '37' pour {GRCh37, hg19, 19}
    """
    s = (gv or "").strip().lower()
    if "grch38" in s or "hg38" in s or s == "38":
        return "38"
    if "grch37" in s or "hg19" in s or s == "19" or s == "37":
        return "37"
    digits = "".join(ch for ch in s if ch.isdigit())
    if digits == "38":
        return "38"
    if digits in {"37", "19"}:
        return "37"
    # défaut raisonnable
    return "38"

def _pattern_for_gv(norm: str) -> str:
    """Regex (lowercase) pour matcher les valeurs de GenomeVersion."""
    if norm == "38":
        return r"(grch38|hg38|38)"
    return r"(grch37|hg19|37|19)"

def _load_cnv_any_header(path: str) -> pl.LazyFrame:
    """
    Charge un TSV CNV (avec/sans en-tête). On prend les 3 premières colonnes
    comme Chr, Start, STOP (cast Start/STOP en int). On conserve toutes les colonnes d'origine.
    """
    # essai avec header
    lf = pl.scan_csv(path, separator="\t", has_header=True)
    cols = lf.collect_schema().names()
    header_looks_fake = any(c.startswith("column_") for c in cols[:min(3, len(cols))])

    if header_looks_fake:
        # relire sans header
        lf = pl.scan_csv(path, separator="\t", has_header=False)
        cols = lf.collect_schema().names()

    if len(cols) < 3:
        raise ValueError("Le fichier CNV doit contenir au moins 3 colonnes (Chr, Start, STOP).")

    c1, c2, c3 = cols[0], cols[1], cols[2]

    lf = (
        lf.with_columns([
            pl.col(c1).alias("Chr"),
            pl.col(c2).cast(pl.Int64, strict=False).alias("Start"),
            pl.col(c3).cast(pl.Int64, strict=False).alias("STOP"),
        ])
        .drop_nulls(["Start", "STOP"])
    )
    return lf

def _load_prob_regions_flex(path: str, genome_version: str) -> pl.LazyFrame:
    """
    Lit --prob_regions (TSV/BED) avec ou sans en-tête.
    - Si en-tête : détecte Chr/Start/(End|STOP).
      * Si 'GenomeVersion' (casse ignorée) existe -> filtre par alias de --genome_version.
      * Sinon -> WARNING et PAS de filtre (on garde toutes les lignes).
    - Si pas d'en-tête : colonnes 1..3 = Chr, Start, STOP.
    Retourne un LazyFrame avec colonnes normalisées Chr, Start, STOP (int).
    """
    norm = _normalize_gv_arg(genome_version)
    pattern = _pattern_for_gv(norm)

    try:
        lf = pl.scan_csv(path, separator="\t", has_header=True)
        cols = lf.collect_schema().names()
        if len(cols) < 3:
            raise ValueError

        chr_col   = _ci_lookup(cols, "Chr")   or cols[0]
        start_col = _ci_lookup(cols, "Start") or cols[1]
        # accepter End ou STOP
        stop_col  = _ci_lookup(cols, "STOP") or _ci_lookup(cols, "End") or cols[2]

        gv_col = _ci_lookup(cols, "GenomeVersion")
        if gv_col:
            # filtre par alias (regex, insensible à la casse via to_lowercase)
            lf = lf.filter(pl.col(gv_col).cast(pl.Utf8).str.to_lowercase().str.contains(pattern))
        else:
            print(
                f"[WARN] 'GenomeVersion' absente dans {path}. "
                f"Aucun filtrage appliqué : toutes les lignes de regions seront utilisées.",
                file=sys.stderr
            )

        lf = (
            lf.select([
                pl.col(chr_col).alias("Chr"),
                pl.col(start_col).cast(pl.Int64, strict=False).alias("Start"),
                pl.col(stop_col).cast(pl.Int64, strict=False).alias("STOP"),
            ])
            .drop_nulls(["Start", "STOP"])
        )
        return lf

    except Exception:
        # pas d'en-tête fiable -> 3 premières colonnes
        lf = pl.scan_csv(path, separator="\t", has_header=False)
        return (
            lf.select([
                pl.col("column_1").alias("Chr"),
                pl.col("column_2").cast(pl.Int64, strict=False).alias("Start"),
                pl.col("column_3").cast(pl.Int64, strict=False).alias("STOP"),
            ])
            .drop_nulls(["Start", "STOP"])
        )

def main():
    parser = argparse.ArgumentParser(description="Build CNV-gene intersection database (STOP + prob_regions partout, sans filtre Region)")
    parser.add_argument("--cnv", help="Path to CNV input TSV", required=True)
    parser.add_argument("--gene_resource", help="Path to gene resource file (BED/TSV)", required=True)
    parser.add_argument("--out", default="data/cnv_geneDB.tsv",
                        help="Output path for CNV-gene database (default: data/cnv_geneDB.tsv)")
    parser.add_argument("--genome_version", default="GRCh38",
                        help="Genome version alias: GRCh38, GRCh37, hg19, hg38, 19, 38 (case-insensitive)")
    parser.add_argument("--prob_regions", required=True, help="Path to problematic_regions TSV/BED")
    parser.add_argument("--bedtools_path", default="bedtools",
                        help="Absolute path to bedtools executable")
    args = parser.parse_args()

    # CNV (flex) + pRegions (flex)
    cnvDB = _load_cnv_any_header(args.cnv)
    pRegions = _load_prob_regions_flex(args.prob_regions, args.genome_version)

    tmpdir = tempfile.mkdtemp()
    try:
        # CNVs -> BED (3 colonnes)
        cnv_bed = f"{tmpdir}/tmp_cnvs.bed"
        cnvDB.select(["Chr", "Start", "STOP"]).unique().sink_csv(
            cnv_bed, separator="\t", include_header=False
        )

        # pRegions -> BED (3 colonnes)
        pRegions_bed = f"{tmpdir}/tmp_pRegions.bed"
        pRegions.select(["Chr", "Start", "STOP"]).sink_csv(
            pRegions_bed, separator="\t", include_header=False
        )

        # Intersect CNV x Genes
        inter_bed = f"{tmpdir}/tmp_intersect.bed"
        subprocess.run(
            f"{args.bedtools_path} intersect -a {cnv_bed} -b {args.gene_resource} -F 0.1 -wao > {inter_bed}",
            shell=True, check=True
        )

        # Intersect CNV x pRegions
        preg_inter_bed = f"{tmpdir}/tmp_pRegion_intersect.bed"
        subprocess.run(
            f"{args.bedtools_path} intersect -a {cnv_bed} -b {pRegions_bed} -wao > {preg_inter_bed}",
            shell=True, check=True
        )

        # Colonnes de l'intersect gènes (identiques à l’original, STOP au lieu de End)
        input_columns = [
            pl.col("column_1").alias("Chr"),
            pl.col("column_2").alias("Start"),
            pl.col("column_3").alias("STOP"),
            pl.col("column_5").alias("t_Start"),
            pl.col("column_6").alias("t_End"),
            pl.col("column_8").alias("gene_type"),
            pl.col("column_9").alias("transcript"),
            pl.col("column_10").alias("gene_name"),
            pl.col("column_12").alias("LOEUF"),
            pl.col("column_13").alias("bp_overlap"),
        ]

        intersected_genes = (
            pl.scan_csv(inter_bed, separator="\t", has_header=False)
              .select(input_columns)
              .with_columns(
                  (pl.col(c).replace([".", "NA"], None) for c in ["gene_type","gene_name","transcript","LOEUF"])
              )
              .with_columns(
                  (pl.col("t_End") - pl.col("t_Start")).alias("transcript_length"),
                  pl.col("LOEUF").cast(pl.Float64).alias("LOEUF"),
              )
        )

        # pRegions overlap fraction, puis MAX par CNV (Chr, Start, STOP)
        p_region_columns = [
            pl.col("column_1").alias("Chr"),
            pl.col("column_2").alias("Start"),
            pl.col("column_3").alias("STOP"),
            pl.col("column_7").alias("bp_overlap"),
        ]
        cnvs_with_pRegion = (
            pl.scan_csv(preg_inter_bed, separator="\t", has_header=False)
              .select(p_region_columns)
              .with_columns(
                  (pl.col("bp_overlap") / (pl.col("STOP") - pl.col("Start") + 1))
                  .alias("cnv_problematic_region_overlap")
              )
              .group_by(["Chr", "Start", "STOP"])
              .agg(pl.max("cnv_problematic_region_overlap").alias("cnv_problematic_region_overlap"))
        )

        # Join final sur (Chr, Start, STOP)
        (
            cnvDB.join(intersected_genes, on=["Chr", "Start", "STOP"], how="left")
                 .join(cnvs_with_pRegion, on=["Chr", "Start", "STOP"], how="left")
                 .sink_csv(args.out, separator="\t")
        )

    finally:
        shutil.rmtree(tmpdir)

if __name__ == "__main__":
    main()

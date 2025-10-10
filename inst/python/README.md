CNV Inheritance Annotation (child vs. parents)

Annotate whether each child CNV is inherited (True) or de novo (False) using reciprocal overlap with parental CNVs, and mark gene-level transmission as True/False or intergenic.

The script expects a gene-annotated CNV table and a pedigree. It outputs CHILD rows only with two added columns:

transmitted_cnv (bool): True = inherited, False = de novo

transmitted_gene (bool / "intergenic"): gene-level inheritance for each child row

Overlap is computed with bioframe.overlap.

#####################################################################################

Requirements

Python ≥ 3.8

Packages:

pandas

numpy

bioframe

Install
python3 -m pip install --upgrade pip
python3 -m pip install pandas numpy bioframe

###############################################################################
##############################################################################

Inputs
1) Gene-annotated CNV table (--cnv_geneAnnot)

Tab-delimited (TSV); .tsv or .tsv.gz both supported.

Expected columns (exact or auto-detected variants in parentheses):

SampleID (sampleid, sample, id, iid)

Chr (chr, chrom, chromosome, contig)

Start (start, pos, begin, position)

End (end, stop, STOP, END)

Type (type, svtype) — values normalized to DEL/DUP (also accepts DELETION/DUPLICATION)

gene_name (optional; gene, symbol)

CNV type must match between child and parent (DEL↔DEL, DUP↔DUP) for inheritance.

2) Pedigree (--pedigree)

Tab-delimited (TSV).

Auto-detects the following (case-insensitive, prefix-match):

Child: ChildID (or SampleID, child, proband, subject, …)

Father: FatherID (or father, biofather, …)

Mother: MotherID (or mother, biomother, …)

Only complete trios are used (ChildID, FatherID, MotherID all present).

A unique TrioKey = ChildID_FatherID_MotherID is created.

#######################################################################
#######################################################################
Output

A TSV at --output containing CHILD rows only, plus:

transmitted_cnv (bool)

transmitted_gene (bool / "intergenic" / NaN)

All original CNV columns are preserved.

######################################################################

Method

Trios: Build TrioKey and long table (TrioKey, SampleID, family_statue ∈ {child,father,mother}).

Join CNVs with trio roles by SampleID. Many-to-many is allowed (parents can appear in multiple trio keys).

Normalize:

Type → DEL/DUP (accepts DELETION/DUPLICATION)

Coordinates to numeric; filter End > Start and Type ∈ {DEL,DUP}.

CNV-level inheritance:

For each unique child CNV (TrioKey, Chr, Start, End, Type), compute overlaps with parent CNVs of the same (TrioKey, Type).

Reciprocal overlap threshold: ov_len/child_len ≥ T and ov_len/parent_len ≥ T where T = --overlap (default 0.5).

Mark transmitted_cnv = True for those child CNVs; else False.

Gene-level inheritance per row:

If gene_name present: for a child row with a gene, set transmitted_gene = True if the same (TrioKey, gene_name, Type) exists in any parent; else False.

If the child row has no gene, set transmitted_gene = "intergenic".

Parent rows are not output.

Overlap definition follows bioframe (0-based, half-open).

##############################################################
#############################################################

Usage

```
python cnv_inheritance.py \
  --cnv_geneAnnot "/path/to/SPARK_CNV_annotated_By_Genes.tsv" \
  --pedigree "/path/to/trios_WGS_SPARK.tsv" \
  --output "annotated_child_cnv.tsv" \
  --overlap 0.5
```

--overlap is the single reciprocal fraction applied to both sides (child and parent).

Use quotes around paths with spaces.

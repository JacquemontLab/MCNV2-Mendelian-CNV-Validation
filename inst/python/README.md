

### Bin Contents


#### CNV Annotation:

annotate_CNV.py 

requires bedtools on system Path
```bash
    python annotate_CNV.py \
        --cnv data/cnv_input.tsv \
        --gene_resource data/gene_annotations.bed \
        --prob_regions data/problematic_regions.tsv \
        --out data/cnv_geneDB.tsv \
        --genome_version GRCh38
```
Will produce a header like this:
```
┌─────────┬───────────┬───────────┬─────────┬───────────┬───────────┬───────────┬────────────────┬─────────────────┬───────────┬────────┬────────────┬───────────────────┬────────────────────────────────┐
│   Chr   │   Start   │    End    │  TYPE   │ SampleID  │  t_Start  │   t_End   │   gene_type    │   transcript    │ gene_name │ LOEUF  │ bp_overlap │ transcript_length │ cnv_problematic_region_overlap │
│ varchar │   int64   │   int64   │ varchar │  varchar  │   int64   │   int64   │    varchar     │     varchar     │  varchar  │ double │   int64    │       int64       │             double             │
├─────────┼───────────┼───────────┼─────────┼───────────┼───────────┼───────────┼────────────────┼─────────────────┼───────────┼────────┼────────────┼───────────────────┼────────────────────────────────┤
│ chr3    │ 195778253 │ 195791547 │ DEL     │ SP0042077 │ 195746771 │ 195811929 │ protein_coding │ ENST00000463781 │ MUC4      │  0.953 │      13294 │             65158 │             0.9888679954870252 │
│ chr22   │  21208539 │  21216717 │ DUP     │ SP0042109 │        -1 │        -1 │                │                 │           │        │          0 │                 0 │             0.9998777356645067 │
└─────────┴───────────┴───────────┴─────────┴───────────┴───────────┴───────────┴────────────────┴─────────────────┴───────────┴────────┴────────────┴───────────────────┴────────────────────────────────┘
```
#### Annotate CNV-Inheritance 

For using just genomic coordinates
cnv_trio_inheritence.py

```bash
python cnv_inheritance.py \
    --cnv data/merged_cnvs.tsv.gz \
    --pedigree data/filtered_pedigree.tsv \
    --output results/annotated_child_cnv.tsv \
    --type_col TYPE \
    --overlap 0.5,0.1
```
Will produce a header like this with one column per requested overlap condition:
```
┌───────────┬─────────┬───────────┬───────────┬─────────┬────────────────────────┐
│ SampleID  │   Chr   │   Start   │    End    │  TYPE   │ Observed_in_Parent_0.5 │
│  varchar  │ varchar │   int64   │   int64   │ varchar │        boolean         │
├───────────┼─────────┼───────────┼───────────┼─────────┼────────────────────────┤
│ SP0042124 │ chr1    │ 152600682 │ 152614148 │ DUP     │ false                  │
│ SP0042124 │ chr1    │ 152600682 │ 152614148 │ DUP     │ false                  │
└───────────┴─────────┴───────────┴───────────┴─────────┴────────────────────────┘
```

For a gene-level inheritance assessment:


```bash
    python get_transmitted_genes.py \
        --pedigree data/filtered_pedigree.tsv \
        --anno_cnv_file data/annotated_cnvs.tsv \  #from annotate_CNV.py
        --output cnv_anno_with_transmitted_gene.tsv
```

Which makes this header:
```
┌─────────┬──────────┬──────────┬─────────┬───────────┬──────────┬──────────┬────────────────┬─────────────────┬───────────┬────────┬────────────┬───────────────────┬────────────────────────────────┬─────────────────┐
│   Chr   │  Start   │   End    │  TYPE   │ SampleID  │ t_Start  │  t_End   │   gene_type    │   transcript    │ gene_name │ LOEUF  │ bp_overlap │ transcript_length │ cnv_problematic_region_overlap │ transmittedGene │
│ varchar │  int64   │  int64   │ varchar │  varchar  │  int64   │  int64   │    varchar     │     varchar     │  varchar  │ double │   int64    │       int64       │             double             │     boolean     │
├─────────┼──────────┼──────────┼─────────┼───────────┼──────────┼──────────┼────────────────┼─────────────────┼───────────┼────────┼────────────┼───────────────────┼────────────────────────────────┼─────────────────┤
│ chr20   │  1603760 │  1619994 │ DEL     │ SP0042125 │  1561385 │  1620009 │ protein_coding │ ENST00000381605 │ SIRPB1    │  1.189 │      16234 │             58624 │             0.7433323067446874 │ true            │
│ chr1    │ 25246394 │ 25307901 │ DEL     │ SP0042125 │ 25242249 │ 25247454 │ protein_coding │ ENST00000243189 │ RSRP1     │  1.423 │       1060 │              5205 │              0.796953241854718 │ true            │
└─────────┴──────────┴──────────┴─────────┴───────────┴──────────┴──────────┴────────────────┴─────────────────┴───────────┴────────┴────────────┴───────────────────┴────────────────────────────────┴─────────────────┘
```
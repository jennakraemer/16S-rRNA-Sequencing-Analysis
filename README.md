# 16S rRNA Sequencing Re-analysis

Re-analysis of Illumina 16S rRNA amplicon sequencing data from Kinzelman *et al.* (2020), focused on microbial community differences across sampling locations at a Lake Michigan beach (Racine, WI).

## Project goals

-   Reproduce and compare key community patterns across locations using an updated ASV-based pipeline.
-   Replace OTU-based processing (QIIME 1 / Greengenes) with DADA2 + SILVA taxonomy, and run modern differential abundance (MaAsLin3).

## Data

-   **NCBI BioProject:** PRJNA640667
-   **Samples:** SAMN15327323–SAMN15327376
-   **Important note:** SRA FASTQs were uploaded as merged reads even though the original run was 2×250 paired-end (MiSeq), so this workflow treats the data as single-end reads.

## Methods overview

### Original study

-   QIIME 1.9.1 OTU clustering, Greengenes 13_5, rarefaction; alpha/beta diversity and LEfSe.

### This re-analysis

-   **ASV inference:** DADA2
-   **Taxonomy:** SILVA v138.2
-   **Alpha diversity:** Shannon (non-parametric permutation t-tests; 999 resamples)
-   **Beta diversity:** ordination (PCoA) and PERMANOVA (RRPP)
-   **Differential abundance:** MaAsLin3

## Key results

-   Location was significant in PERMANOVA (p = 0.001), consistent with the original report.
-   Differential abundance results supported enrichment of **Arcobacter** in stormwater (IEB Overflow) relative to several other sites; patterns for **Flavobacterium** differed from the paper in this re-analysis.

## Repository contents (typical)

-   `*.Rmd` / `*.md`: analysis writeup(s)
-   `figures/`: plots and paper figure reproductions referenced by the report

## Reproducibility

### Requirements

-   R/RStudio
-   R packages used in the analysis include: `dada2`, `ShortRead`, `Biostrings`, `tidyverse`, `data.table`, `vegan`, `RRPP`, `coin`, `multcompView`, `phyloseq`, and `microbiome`.
-   Download FASTQs for BioProject PRJNA640667 from SRA BioProject and place them in local fastq/ directory.
-   Update file paths inside the `.Rmd` to point to your FASTQ directory and your SILVA training set.

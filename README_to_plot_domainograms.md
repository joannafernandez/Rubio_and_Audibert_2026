## Domainograms:

#### Processing raw FASTQs:  
Refer to https://github.com/joannafernandez/cut-run_for_soutoglou_laboratory for information on how the Soutoglou Laboratory processes CUT&RUN sequencing data

Using nf-core cutandrun v 3.2.2  
https://github.com/nf-core/cutandrun

> Tamara L. Hodgetts, Charlotte West, Nicholas M. Luscombe, Jernej Ule, James Briscoe& Chris Cheshire. (2024). nf-core/cutandrun: A Nextflow pipeline for the analysis of CUT&RUN, CUT&Tag and TIP-seq datasets (Version 3.2.2) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.10606804

> The nf-core framework for community-curated bioinformatics pipelines.
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.

The nf-core cutandrun pipeline was run with the following command:

```bash
nextflow run nf-core/cutandrun -r 3.2.2 -profile docker \
  --input ~/samplesheet.csv \
  --genome hg38 \
  --normalisation_mode CPM \
  --normalisation_binsize 1 \
  --outdir /media/bafj7/jo4/laminnospike \
  --work-dir ~/nf_work \
  --save_reference \
  --save_trimmed \
  --use_control true \
  --max_memory '24.GB' \
  --max_cpus 8 \
  --cores 30 \
  -with-report report.html \
  -with-trace trace.txt \
  -with-timeline timeline.html \
  -with-dag dag.png

```

---
#### Downstream 
1. igg_adjust_LB1.sh (run on CPM normalised bigwigs from nfcore cutandrun pipeline)
2. table generated with multibigwigsummary:

```bash
  multiBigwigSummary bins \
  --bwfiles \
  /Volumes/jo4/laminnospike/03_peak_calling/03_bed_to_bigwig/derived_bigwigs/LB1_DMSO_CTRL_R1.adj_log2.bw \
  /Volumes/jo4/laminnospike/03_peak_calling/03_bed_to_bigwig/derived_bigwigs/LB1_DMSO_LAD_R1.adj_log2.bw \
  /Volumes/jo4/laminnospike/03_peak_calling/03_bed_to_bigwig/derived_bigwigs/LB1_DMSO_CTRL_R2.adj_log2.bw \
  /Volumes/jo4/laminnospike/03_peak_calling/03_bed_to_bigwig/derived_bigwigs/LB1_DMSO_LAD_R2.adj_log2.bw \
  /Volumes/jo4/laminnospike/03_peak_calling/03_bed_to_bigwig/derived_bigwigs/LB1_DMSO_CTRL_R3.adj_log2.bw \
  /Volumes/jo4/laminnospike/03_peak_calling/03_bed_to_bigwig/derived_bigwigs/LB1_DMSO_LAD_R3.adj_log2.bw \
  --binSize 5000 \
  --outRawCounts /Volumes/jo4/laminnospike/03_peak_calling/03_bed_to_bigwig/derived_bigwigs/DMSO_5kb_rawcounts.tsv \
  -o /Volumes/jo4/laminnospike/03_peak_calling/03_bed_to_bigwig/derived_bigwigs/DMSO_5kb_summary.npz
```

3. Plotted with python (see Figure3c.py)

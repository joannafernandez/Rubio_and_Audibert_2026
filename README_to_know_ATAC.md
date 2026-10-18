## README_to_know_ATAC.md

   > read me for direction on ATAC-seq pipeline used and plotting
   > ATAC_differentialsdeltadelta.py used in Fig 4a, Ext Fig 6 i,j,k for delta plots
   > Figure4a_ExtFig6ijk_heatmaps.py is script used to plot heatmaps only

   > nb disregard ATAC_heatmap.py and barebones.py, outdated


---
#### Data processing 

- **nf-core pipeline introduction** —
  https://nf-co.re/atacseq/2.1.2/

- **nf-core ATAC-seq 2.1.2 usage** —
  https://nf-co.re/atacseq/2.1.2/docs/usage/

- **nf-core ATAC-seq 2.1.2 parameters** —
  https://nf-co.re/atacseq/2.1.2/parameters/

References: 
> The nf-core framework for community-curated bioinformatics pipelines.
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.

>Harshil Patel, Jose Espinosa-Carrasco, Björn Langer, Phil Ewels, nf-core bot, Maxime U Garcia, Robert Syme, Alexander Peltzer, Adam Talbot, Drew Behrens, Gisela Gabernet, Mingda Jin, Matthias Hörtenhuber, Jonatan Gonzalez Rodriguez, Kevin Menden& Ömer An. (2023). nf-core/atacseq: [2.1.2] - 2022-08-07 (Version 2.1.2) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.8222875


---
The nf-core atac-seq pipeline was run with the following command:

```bash

nextflow run nf-core/atacseq -r 2.1.2 -profile docker \
  --input ~/atac_samplesheet.csv \
  --genome hg38 \
  --read_length 150 \
  --outdir /media/bafj7/jo5/atac/atacseq_results \
  --work-dir ~/nf_work \
  --save_reference \
  --save_trimmed \
  --max_memory '24.GB' \
  --max_cpus 8 \
  --cores 30 \
  -with-report atacseq_report.html \
  -with-trace atacseq_trace.txt \
  -with-timeline atacseq_timeline.html \
  -with-dag atacseq_dag.png

```

# Rubio_and_Audibert_compile
Compilation of python and shell scripts used for analysis and plotting  in the article "DNA damage at Lamina-Associated Domains triggers nuclear envelope reorganization and chromatin detachment to prevent nuclear envelope blebbing and genome instability." 



---

### Directory:

1. README_to_plot_domainograms.md
   > read me to plot domainograms
   >
   > requires use of igg_adjust_LB1.sh
   >
   > requires use of Figure3c.py for plotting

2. Fig2g_means.py
   > Averaged bigwigs generated using bw_averaging.sh
   >
   > Plotting averaged bigwigs (not spike normalised) of LaminB1 and 53BP1 signal in DMSO condition (with and without DNA damage)
   >
   
3. TIDEKO_verify_withrelativepercentage.py
   > plotting TIDE data (relative summation of significant indels only) of knock out pools with various single guide RNAs
   >
  

4. ExtendedFig5f_pATMyH2AXcoloc.py
   > analysis and plotting of per nuclei and foci colocalisation rates between yH2AX and pATM in WT and 53BP1-/- cells
   >

5. NEW_ExtendedFig7.py (disregard ExtendedFigt_TIDE.py, outdated)
   > Processing raw TIDE data and plotting for sgLAD and sgiLAD experiments

6. README_to_know_ATAC.md
   > read me for direction on ATAC-seq pipeline used and plotting
   >
   > ATAC_differentialsdeltadelta.py used in Fig 4a, Ext Fig 6 i,j,k for delta plots
   >
   > Figure4a_ExtFig6ijk_heatmaps.py is script used to plot heatmaps only
   >
   

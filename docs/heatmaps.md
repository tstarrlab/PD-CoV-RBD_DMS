---
layout: heatmaps
permalink: /RBD-heatmaps/
---

---

*NOTE data are preliminary*

### Overview

You can use this tool to explore the experimentally determined impacts of amino acid mutations on APN binding in the PD-CoV receptor-binding domain (RBD). 

#### Instruction

To use this tool, select the metric that you wish to display in the heatmap (either change in galline APN binding affinity (-log10 $$K_D$$), or change in mean fluorescence intensity (MFI) for human and porcine APN as determined at 1uM ligand incubation, by selecting that metric in the corresponding drop down menu. Hover over individual mutations to see exact numerical details.

#### Technical Details

The impact on APN receptor-binding of every single amino-acid mutation in PD-CoV RBDs, as determined by high-throughput FACS-seq assays. Wildtype amino acids are indicated by an 'x', and gray squares indicate missing mutations from each library. The number of internally replicated barcodes with which a mutation was measured is visible as `Barcode Count` in the tooltips, where higher numbers indicate higher-confidence measurements.

Galline APN binding was determined from full titration curves and thus is analogous to a Kd measurement. Human and porcine APN binding was more sparse in the library, so we did not collect entire titrations, but rather just determined mutant MFI at a 1uM ligand incubation.

### Data

Raw data  can be found [here](https://github.com/tstarrlab/PD-CoV-RBD_DMS/blob/main/results/final_variant_scores/final_variant_scores.csv). The code used to make these plots can be found [here](https://github.com/tstarrlab/PD-CoV-RBD_DMS/blob/main/RBD-Heatmaps-Interactive-Visualization.ipynb).

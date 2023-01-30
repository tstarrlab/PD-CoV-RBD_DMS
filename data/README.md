# Input data
This directory contains input data for the analysis.

## Basic information about sequencing mapping and parsing specs

These files are used for the basic processing of the deep sequencing data to call variants by barcode and count barcodes:

   - `PacBio_amplicon_PDCoV.gb`: the amplicon being sequenced by PacBio.

   - `feature_parse_specs_PDCoV.yaml`: how to parse the amplicon when handling the PacBio data.

   - [PacBio_runs.csv](PacBio_runs.csv): list of the PacBio runs used to call the variants.

   - [barcode_runs.csv](barcode_runs.csv): list of the Illumina runs used to count the barcodes for different samples.

## Structures

   - [6m0j.pdb](6m0j.pdb): ACE2-bound SARS-CoV-2 RBD structure
   
   - [./VOC_structures](./VOC_structures): structurally aligned RBDs of ACE2-bound RBD structures for the Wuhan-Hu-1, alpha, beta, and omicron structures. See the description in the [structural_shifts](../structural_shifts.Rmd) analysis of how these files were created.

## Reference sequences

   - Genbank plasmid maps `171` and `171lib` showing the parental PD-CoV yeast display expression construct and the library version that contains the architecture with the inserted N16 barcode fragment
   
   - [dcov_RBD_aligned.fasta](./dcov_RBD_aligned.fasta): our alignment of porcine and avian deltacoronavirus sequences
   
   - [RBD_sites.csv](./RBD_sites.csv): table giving mapping of alignment, RBD, and spike indexed numbering for our PD-CoV construct

   
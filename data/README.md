# Input data
This directory contains input data for the analysis.

## Basic information about sequences and alignments

These files are used for the basic processing of the deep sequencing data to call variants by barcode and count barcodes:

   - [PacBio_amplicons.gb](PacBio_amplicons.gb): the amplicons being sequenced by PacBio.
     Note that there are a variety of possible amplicons because in addition to the SARS-CoV-2 RBD there are RBDs from a variety of other viral strains.

   - [feature_parse_specs.yaml](feature_parse_specs.yaml): how to parse the amplicon when handling the PacBio data.

   - [PacBio_runs.csv](PacBio_runs.csv): list of the PacBio runs used to call the variants.

   - [barcode_runs.csv](barcode_runs.csv): list of the Illumina runs used to count the barcodes for different samples. This file is built from [barcode_runs_orig-names.csv](barcode_runs_orig-names.csv) by the Jupyter notebook [build_barcode_runs.ipynb](build_barcode_runs.ipynb).

   - [RBD_sites.csv](RBD_sites.csv): gives site and residue information for SARS-CoV-2, including alignment of the RBD integer numbering with the Spike numbering for SARS-CoV-2 RBD, alignment to SARS-CoV, and structural annotations as detailed below.

## Alignments of different Spikes / RBDs

## Input data from previous studies

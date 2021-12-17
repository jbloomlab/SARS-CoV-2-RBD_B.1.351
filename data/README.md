# Input data
This directory contains input data for the analysis.

## Basic information about sequences and alignments

These files are used for the basic processing of the deep sequencing data to call variants by barcode and count barcodes:
  
   - [wildtype_sequence.fasta](wildtype_sequence.fasta): The sequence of the unmutated B.1.351 RBD.

   - [wuhan1_wildtype_sequence.fasta](wuhan1_wildtype_sequence.fasta): The sequence of the unmutated Wuhan-Hu-1 RBD.

   - [RBD_sites.csv](RBD_sites.csv): gives site and residue information for SARS-CoV-2, including alignment of the RBD integer numbering with the Spike numbering for SARS-CoV-2 RBD, alignment to SARS-CoV, and structural annotations as detailed below.

## Files related to processing NGS sequencing data

  - [PacBio_amplicons.gb](PacBio_amplicons.gb): the amplicons being sequenced by PacBio.
    Note that there are a variety of possible amplicons because in addition to the SARS-CoV-2 RBD there are RBDs from a variety of other viral strains.

  - [feature_parse_specs.yaml](feature_parse_specs.yaml): how to parse the amplicon when handling the PacBio data.

  - [PacBio_runs.csv](PacBio_runs.csv): list of the PacBio runs used to call the variants.

  - [barcode_runs.csv](barcode_runs.csv): list of the Illumina runs used to count the barcodes for different samples. This file is built from [barcode_runs_orig-names.csv](barcode_runs_orig-names.csv) by the Jupyter notebook [build_barcode_runs.ipynb](build_barcode_runs.ipynb).

  - [barcode_runs_helper_files](barcode_runs_helper_files): subdirectory containing auxillary notebooks and files to help generate the `barcode_runs.csv` file that maps sample names to file paths on the Hutch cluster.

## Plasmid maps

[plasmid_maps](plasmid_maps/): This subdirectory contains the full Genbank maps for key plasmids used in the study:

      - [2800_HDM_Spikedelta21_D614G.gb](plasmid_maps/2800_HDM_Spikedelta21_D614G.gb): the wildtype D614G and

      - [2957_HDM_Spikedelta21_B.1.351.gb](plasmid_maps/2957_HDM_Spikedelta21_B.1.351.gb): B.1.351 spike expression plasmids used to make pseudotyped lentiviral particles

      - [3021_pETcon-SARS-CoV-2-RBD_K417N_E484K_N501Y.gb](plasmid_maps/3021_pETcon-SARS-CoV-2-RBD_K417N_E484K_N501Y.gb): the B.1.351 wildtype RBD yeast-display plasmid, and

      - [pETcon-SARS-CoV-2-RBD-B1351_lib-assembled.gb](plasmid_maps/pETcon-SARS-CoV-2-RBD-B1351_lib-assembled.gb): the fully assembled plasmid, including the Illumina adaptors and Nx16 barcode, that was used as the template for designing the Twist site-saturation variant libraries.

## For visualizing serum-escape data:

These files are used for visualizing the antibody- or serum-escape data:

  - [site_color_schemes.csv](site_color_schemes.csv): Schemes for how to color sites (can be used in escape profiles). Here are details on these schemes.

  - [escape_profiles_config.yaml](escape_profiles_config.yaml): Information on how to plot the escape profiles; manually edit this to alter their plotting.

  - [early2020_escape_profiles_config.yaml](early2020_escape_profiles_config.yaml): Same as above, but for early 2020 plasmas mapped against the Wuhan-Hu-1 RBD library.

  - [lineplots_config.yaml](lineplots_config.yaml): Config file for making line plots to compare serum-escape scores between cohorts.

  - [dms-view_metadata.md](dms-view_metadata.md): Used for rendering the dms-view page to visualize data.

## Input data from previous studies

  - [final_variant_scores.csv](final_variant_scores.csv): Comprehensive measurement of RBD expression and ACE2 binding as determined by deep mutational scanning experiments in a site-saturation single-mutant library in the B.1.351 RBD background similar to those described in [Starr, et al. (2020)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7418704/). These scores are used to filter variants that are highly deleterious for ACE2 binding or RBD expression before calculation of antibody-escape scores in the [counts_to_scores.ipynb](counts_to_scores.ipynb) notebook.

## Data from GISAID surveillance

  - [210801_mutation_counts.csv](210801_mutation_counts.csv): The counts of all RBD mutations in the SARS-CoV-2 spike alignment downloaded from [GISAID](https://www.gisaid.org/) on Aug. 1, 2021. These data are used to help set filters on the ACE2 binding and RBD expression scores, to ideally retain mutations that are observed an appreciable number of times in GISAID sequences.

## PDB files in [pdbs](pdbs/) subdirectory

  - [6M0J](pdbs/6M0J.pdb): Wuhan-Hu-1 RBD bound to huACE2.

  - [7LYQ](pdbs/7LYQ.pdb): B.1.351 spike trimer

  - [7LYQ_RBD_ACE2](pdbs/7LYQ_RBD_ACE2.pdb): RBD from 7LYQ aligned with 6M0J and shown bound to ACE2

  - [7LYQ_RBD](pdbs/7LYQ_RBD.pdb): RBD from 7LYQ

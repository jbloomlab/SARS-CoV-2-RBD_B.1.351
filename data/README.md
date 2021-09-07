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
   
   - [plasmid_maps](plasmid_maps/): This subdirectory contains the full Genbank maps for key plasmids used in the study: 
    - the wildtype D614G and B.1.351 spike expression plasmids used to make pseudotyped lentiviral particles, 
    - the B.1.351 wildtype RBD yeast-display plasmid, and 
    - the fully assembled plasmid, including the Illumina adaptors and Nx16 barcode, that was used as the template for designing the Twist site-saturation variant libraries.

## For visualizing serum-escape data:

These files are used for visualizing the antibody- or serum-escape data:

  - [site_color_schemes.csv](site_color_schemes.csv): Schemes for how to color sites (can be used in escape profiles). Here are details on these schemes.

    - The `subdomain` scheme colors sites orange if they directly contact ACE2 (within 4 angstroms in PDB 6m0j) in the SARS-CoV-2 structure (residues 417, 446, 449, 453, 455, 456, 475, 486, 487, 489, 493, 496, 498, 500, 501, 502, 505), blue if they are in the receptor binding motif (RBM, residue 437 to 508, inclusive), and green if they are in the core RBD domain (all other sites). These definitions match those used in [Starr et al (2020)](https://www.cell.com/cell/fulltext/S0092-8674(20)31003-5).
  
  - [escape_profiles_config.yaml](escape_profiles_config.yaml): Information on how to plot the escape profiles; manually edit this to alter their plotting.

## Alignments of different Spikes / RBDs

## Input data from previous studies

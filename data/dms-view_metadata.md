## A SARS-CoV-2 variant elicits an antibody response with shifted immunodominance hierarchy

For experimental background, see our paper **here (update with link)**.

### What data are shown here?
We are showing mutations to the SARS-CoV-2 RBD that escape binding by polyclonal antibodies from individuals infected with early 2020 or B.1.351 viral isolates, measured using mutational antigenic profiling. Raw data are available raw data are available [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_B.1.351/blob/main/results/supp_data/B1351_raw_data.csv).
The drop-down menus can be used to select the escape-mutation maps for each antibody or plasma.

When you click on sites, they will be highlighted on the protein structure of the ACE2-bound Wuhan-Hu-1 RBD structure ([PDB 6M0J](https://www.rcsb.org/structure/6M0J)) or to the RBD of the B.1.351 variant ([PDB 7LYQ](https://www.rcsb.org/structure/7LYQ)).

At the site level you can visualize one of two quantities:

 - *total escape* is the sum of the escape from all mutations at a site.

 - *max escape* is the magnitude of the largest-effect escape mutation at each site.

At the mutation level, the height of each letter is proportional to the extent to which that amino-acid mutation escapes antibody binding.
You can color the logo plot letters in four ways:

 - *escape color ACE2 bind* means color letters according to how that mutation affects ACE2 binding, with yellow meaning highly deleterious, and brown meaning neutral or beneficial for ACE2 binding.

 - *escape color RBD expr* means color letters according to how that mutation affects RBD expression.

 - *escape color gray* means color all letters gray.

 - *escape color func group* means color letters by functional group.

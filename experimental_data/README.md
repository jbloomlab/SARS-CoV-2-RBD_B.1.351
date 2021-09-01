# Experimental (i.e., non-DMS) data for B.1.351 plasma-escape study
This directory contains the data and code for analyzing ELISA and neutralization assays to accompany the B.1.351 plasma-escape mapping by DMS.

Note that all experiments with B.1.351 infection-elicited plasmas are done with B.1.351 RBD and spike proteins, and B.1.351 spike-pseudotyped lentiviral particles.

Similarly, all experiments with "early 2020" plasmas are done with Wuhan-Hu-1 RBD and spike proteins, and D614G spike-pseudotyped lentiviral particles.

**Please note** that all early 2020 experiments were previously performed and described in [Greaney, et al. _Cell Host & Microbe_ (2021)](https://www.sciencedirect.com/science/article/pii/S1931312821000822), with the exception of the point-mutant neutralization assays analyzed in [analyze_neut_data.ipynb](analyze_neut_data.ipynb), which were newly performed in for study.  

### Input data
- The [./data/entry_titers/](data/entry_titers) subdirectory contains the files for calculating the entry titers of pseudotyped lentiviral particles used in neutralization assays. This is primarily used for determining how much lentiviral supernatant to use in each experiment, to target an appropriate amount of relative luciferase units per well.

- The [./data/rbd_depletion_elisas/](data/rbd_depletion_elisas) subdirectory contains files with OD450 readings before and after depletion of RBD-binding antibodies.

- [./data/haarvi_rbd_depletion_foldchange_ic50.csv](data/haarvi_rbd_depletion_foldchange_ic50.csv) is a supplementary table from [Greaney, et al. _Cell Host & Microbe_ (2021)](https://www.sciencedirect.com/science/article/pii/S1931312821000822) that has the IC50 values and fold-change IC50s pre- and post-depletion of RBD-binding antibodies for the early 2020 plasmas.

- The [./neut_data/](neut_data) subdirectory contains the raw Excel files and plate maps for each day of assays.

    ### Formatting data for neutralization assays
    Each neutralization assay should have its data stored in an Excel file in a subdirectory of [./neut_data/](neut_data) named with the data in the format `2020-10-02`.
    The subdirectories should also contain a `sample_map.csv` file that maps the Excel file data to the samples in a format that is readable by the Python script [excel_to_fractinfect.py](excel_to_fractinfect.py) (see [here](https://github.com/jbloomlab/exceltofractinfect) for more on this script, which was written by Kate Crawford).
    The plate layouts referred to by the sample maps are in [./PlateLayouts/](PlateLayouts).

### Running the code
Now you can run the entire analysis.
The analysis consists primarily of a series of Jupyter notebooks along with some additional code in [Snakefile](Snakefile).
You can run the analysis by using [Snakemake](https://snakemake.readthedocs.io) to run [Snakefile](Snakefile):

    snakemake --j1

## Configuring the analysis
The configuration for the analysis is specified in [config.yaml](config.yaml).

## Notebooks that perform the analysis
- [calculate_entry_titer.ipynb](calculate_entry_titer.ipynb): calculates PV entry titers in RLU.

- [rbd_depletion_elisas.ipynb](rbd_depletion_elisas.ipynb): analyzes ELISAs that measure binding of B.1.351 plasma to B.1.351 RBD and spike before and after depletion of B.1.351 RBD-binding antibodies.

- [rbd_depletion_neuts.ipynb](rbd_depletion_neuts.ipynb): analyzes assays that measure neutralization of B.1.351 spike-pseudotyped lentiviral particles before and after depletion of B.1.351 RBD-binding antibodies.

- [analyze_neut_data.ipynb](analyze_neut_data.ipynb): Analyzes point-mutant neuts.



## Results
Results are placed in the [./results/](results) subdirectory.

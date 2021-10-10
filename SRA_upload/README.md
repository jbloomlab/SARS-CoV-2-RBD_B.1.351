# Uploading FASTQ files to the SRA

Details of how the raw sequencing files (FASTQ files) were uploaded to the NIH [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra).
The submitted files are in BioProject [xxx](https://www.ncbi.nlm.nih.gov/bioproject/xxx).

The Python Jupyter notebook [upload_to_SRA.ipynb](upload_to_SRA.ipynb) has instructions and does the uploading.

Because the FTP upload takes a while, you may want to run the Jupyter notebook using `slurm` so there is no timeout with::

    sbatch --wrap="jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 upload_to_SRA.ipynb" --time 2-0


## Manual entries into the SRA submission portal

### Step 3: Project info
**Project title:** Deep mutational scanning of the RBD from SARS-CoV-2 variant lineages using yeast display of barcoded libraries

**Public description:** We performed deep mutational scanning of the SARS-CoV-2 spike receptor-binding domain (RBD) from variant viral lineages using yeast display to measure the effects of all amino-acid mutations on various phenotypes, such as on binding of antibodies or polyclonal plasmas.

### Step 4: Biosample type
**NCBI packages:** Microbe

### Step 5: Biosample attributes
Note: BioSamples must have unique attributes, excluding title, sample_name, description, and accession. So I had to create a new column, that I called `sequencing platform` to distinguish the two BioSamples I wanted to upload.

**PacBio samples**
* sample_name: PacBio_CCSs
* sample_title: PacBio CCSs linking variants to barcodes for SARS-CoV-2 B.1.351 RBD deep mutational scanning
* organism: Severe acute respiratory syndrome coronavirus 2
* strain: B.1.351
* isolation_source: plasmid
* collection_date: 2021
* geo_loc_name: USA
* sample_type: plasmid
* description: long-read PacBio CCSs linking variants to barcodes for SARS-CoV-2 B.1.351 RBD deep mutational scanning
* sequencing platform: PACBIO_SMRT

**Illumina barcode sequencing samples**
* sample_name: B1351_plasma_escape_barcodes		
* sample_title: Illumina barcode sequencing for B.1.351 mutational antigenic profiling
* organism: Severe acute respiratory syndrome coronavirus 2
* strain: B.1.351
* isolation_source: plasmid
* collection_date: 2021
* geo_loc_name: USA
* sample_type: plasmid
* description: Illumina barcode sequencing for B.1.351 mutational antigenic profiling
* sequencing platform: ILLUMINA

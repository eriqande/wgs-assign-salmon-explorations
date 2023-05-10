wgs-assign-salmon-explorations
================
Eric C. Anderson

This repo/directory/Rstudio-project contains the code used to assess
WGSassign on existing Chinook salmon data.

The process is entirely specified in the Snakefile.

The file `001-doing-it.Rmd` holds some of my notes on how the input
files were created.

The BAM files needed as input can be created from the fastq files from
Thompson et al. (2020),
<https://www.science.org/doi/10.1126/science.aba9059>, which are
available from the short read archive:
<https://www.ncbi.nlm.nih.gov/sra/PRJNA667732>. The VCF file used as
input to create reference/{mprun}/reference-beagle-gl.gz was created by
running GATK for variant calling upon the BAMs. These files are too
large to store in this repository.

Should you want access to them without running an entire bioinformatic
pipeline, please contact `eric.anderson@noaa.gov`.

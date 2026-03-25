# Local Differential Editing Index (LoDEI)

## General Notes

LoDEI - the local differential editing index - offers a collection of programs to detect and analyze differentially edited A-to-I regions in two sets of RNA-seq samples.

```
lodei -h             # get a list of all commands 
lodei subcommand -h  # get info of a subcommand
```

The subcommand to detect differential A-to-I editing is `lodei find`.

Analyzing RNA-seq data typcially requires the mapping of NGS reads in fastq format to a reference genome. 
The primary input for `lodei find` are sorted BAM files as produced by NGS-read mappers like STAR.


LoDEI is free software and licensed under GPLv3.

If you use LoDEI, please cite https://doi.org/10.1038/s41467-024-53298-y

### System Requirements

* Linux operating system (our systems run on Ubuntu 22.04)
* conda/mamba *or* Podman/Docker

### Installation

Prior to installation, we recommend to follow the instructions of the provided [test data](#test-data) to be able to verify a proper installation.

Install LoDEI by using one of the following ways:

1. use the [conda/mamba](#installation-and-usage-via-conda) package manager to install LoDEI.
3. [build](#build-the-image-via-the-containerfile) a Podman/Docker image locally by using the provided Containerfile.


## Publication

Torkler, P., Sauer, M., Schwartz, U. et al. LoDEI: a robust and sensitive tool to detect transcriptome-wide differential A-to-I editing in RNA-seq data. Nat Commun 15, 9121 (2024). https://doi.org/10.1038/s41467-024-53298-y

## Test Data

We provide a small test dataset (~15MB, https://zenodo.org/doi/10.5281/zenodo.10907019) that contains all required input files to run `lodei find` to demonstrate the proper usage for detecting differentially edited A-to-I regions.
The test dataset contains sorted BAM files belonging to two different conditions that are thought to be compared against each other, genomic annotations, and the nucleotide sequences for three genes of the human genome.

Let's create an example directory and download the testdata:

```
cd ~            # change to your home directory
mkdir example   # create a new directory to store example data
cd example      # switch to the example directory
# download and unpack test data:
wget https://zenodo.org/records/10907020/files/test_data.tar.gz
tar -xzf test_data.tar.gz
```

After unpacking, the directory `data_testrun` should be in your `example` directory (see below).
The subdirectory `data_testrun/annotation` contains genomic sequences in the fasta format and genomic annotations in the GFF3 format.

The `data_testrun/bam` subdirectory contains BAM files for 10 samples whereas samples 01-05 belong to set 1 and samples 06-10 belong to set 2.

*DO NOT CHANGE ANYTHING IN THE `data_testrun` DIRECTORY*

```
user@linux ~/example $ tree
.
|-- data_testrun
|   |-- annotation
|   |   |-- example_plot_regions.txt
|   |   |-- genome.fa
|   |   |-- genome.fa.fai
|   |   `-- test_anno.gff3
|   `-- bam
|       |-- s01.bam
|       |-- s01.bam.bai
|       |-- s02.bam
|       |-- s02.bam.bai
|       |-- s03.bam
|       |-- s03.bam.bai
|       |-- s04.bam
|       |-- s04.bam.bai
|       |-- s05.bam
|       |-- s05.bam.bai
|       |-- s06.bam
|       |-- s06.bam.bai
|       |-- s07.bam
|       |-- s07.bam.bai
|       |-- s08.bam
|       |-- s08.bam.bai
|       |-- s09.bam
|       |-- s09.bam.bai
|       |-- s10.bam
|       `-- s10.bam.bai
`-- test_data.tar.gz
```

## Installation and Usage via conda

Generate a new environment and install LoDEI:

```
conda create -c conda-forge -c bioconda --name lodei
conda activate lodei
conda install lodei
```

### Run LoDEI Using Conda

To verify a proper installation we run LoDEI on the provided [test dataset](#test-data).

Make sure to get back into your example directory where your unpacked the [test dataset](#test-data). 
Create a new output directory at `~/example` where LoDEI can save the results and finally move into the example data directory:

```
cd ~/example
mkdir output_conda
cd data_testrun
```

Run LoDEI on the testdata:

```
lodei find \
--group1 bam/s01.bam bam/s02.bam bam/s03.bam bam/s04.bam bam/s05.bam \
--group2 bam/s06.bam bam/s07.bam bam/s08.bam bam/s09.bam bam/s10.bam \
-f annotation/genome.fa \
-g annotation/test_anno.gff3 \
-o ../output_conda \
-c 3 --library SR --min_coverage 5
```

Detailed explanation of parameters and arguments:

* `--group1 ...` 
  * provide the list of sorted BAM files (separated by space) of samples belonging to group 1. Note, for each input`.bam` file a corresponding `.bai` file is required to be present in the same directory. 
* `--group2 ...`
  * provide the list of sorted BAM files (separated by space) of samples belonging to group 2. Note, for each input`.bam` file a corresponding `.bai` file is required to be present in the same directory.
* `-f annotation/genome.fa`
  * provide the reference genome used to generate the provided BAM files.
* `annotation/test_anno.gff3`
  * LoDEI caculates differential editing for all entries of the provided [annotation](#what-kind-of-annotation-file-do-i-need-to-use) file.
* `-o ../output_conda`
  * define the output directory. LoDEI generates many automatically named files.
* `-c 3`
  * Number of used CPU cores.
* `--library SR`
  * provide the strandedness of your BAM files. Strandedness is defined as in `salmon` see: https://salmon.readthedocs.io/en/latest/library_type.html  
    Currently, the following types are supported:  
    SR = reverse stranded,  
    SF = forward stranded,  
    U = unstranded  
    ISR = paired-end, reverse stranded  
    ISF = paired-end, forward stranded,  
    If you are unsure what kind of library type (strandedness) your data is, have a look at the [FAQ](#how-do-i-infer-the-library-type-(strandedness)-of-my-data) and 
    https://github.com/rna-editing1/getlibtype  
* `--min_coverage 5`
  * only consider single positions that have a coverage >= min_coverage in all samples.
* [Should I use](#should-i-use---rm_snps) `--rm_snps`?

Wait until LoDEI finishes the calculation (~1-2min) and have a look at the [output](#output).

## Installation and Usage via Podman

Common Linux distributions are typically shipped with Podman. Podman is a tool to create, run and maintain containers.
For a detailed introduction of Podman we refer the reader to the primary documentation at https://podman.io.

### Build the image via the Containerfile

Let's build the image locally:

```
cd ~  # enter your home directory 
git clone https://github.com/rna-editing1/lodei.git # get the repository
cd lodei
podman build -f Containerfile -t lodei
```

### Usage via Podman

Verify that podman is able to start LoDEI by trying to run the new container:

```
podman run -it --rm localhost/lodei:latest lodei find -h
```

If your container runs successfully, 
you should see the help page of LoDEI.


### Run LoDEI Using Podman

To verify a proper installation we run LoDEI on the provided [test dataset](#test-data).

#### Mount a volume/directory into the container

The LoDEI container needs access to the provided files (annotations and bam files) as well as a directory where it can save results to.
The `-v` option is needed to make directories of your host file system available in the container.
In a nutshell, `-v` mounts directories of your file system into the container.
The general syntax is

```
-v /path/on/host/system:/path/in/container:option
```

Option can be `ro` for _read only_ and `rw` for _read and write_.

Note, the directory in the container does not need to exist there. You can specify any directory.

#### Run LoDEI

If you've followed the steps of the [test dataset](#test-data) the directory `~/example` exists.
Switch to the `~/example` directory and create a new directory where LoDEI shall save all output into:

```
cd ~/example
mkdir output_test
```

Next, we will apply LoDEI on the test dataset via calling

```
podman run --rm \
-v ~/example/data_testrun/annotation:/annotation:ro \
-v ~/example/data_testrun/bam:/bam:ro \
-v ~/example/output_test:/output:rw \
localhost/lodei:latest lodei find \
--group1 /bam/s01.bam /bam/s02.bam /bam/s03.bam /bam/s04.bam /bam/s05.bam \
--group2 /bam/s06.bam /bam/s07.bam /bam/s08.bam /bam/s09.bam /bam/s10.bam \
-f /annotation/genome.fa \
-g /annotation/test_anno.gff3 \
-o /output \
-c 3 --library SR --min_coverage 5
```

Detailed explanation of parameters and arguments:

* `-v ~/example/data_testrun/annotation:/annotation:ro`
    * mount the host directory `~/example/data_testrun/annotation` to the directoy `/annotation` in the container with read only permission.
* `-v ~/example/data_testrun/bam:/bam:ro`
    * mount the host directory `~/example/data_testrun/bam` to the directoy `/bam` in the container with read only permission.
* `-v ~/example/output_test:/output:rw`
    * mount the host directory `~/example/output_test` to the directoy `/output` in the container with read and write permissions.
* `localhost/lodei_0.0.1:latest lodei find`
    * `localhost/lodei_0.0.1:latest` is the name of the image from which a new container shall be started followed by the command line call to start `lodei find`.
* `--group1 ...` 
    * provide the list of sorted BAM files (separated by space) of samples belonging to group 1. Note, for each input`.bam` file a corresponding `.bai` file is required to be present in the same directory.
* `--group2 ...`
    * provide the list of sorted BAM files (separated by space) of samples belonging to group 2. Note, for each input`.bam` file a corresponding `.bai` file is required to be present in the same directory.
* `-f /annotation/genome.fa`
    * provide the reference genome used to generate the provided BAM files.
* `/annotation/test_anno.gff3`
    * LoDEI caculates differential editing for all entries of the provided [annotation](#what-kind-of-annotation-file-do-i-need-to-use) file.
* `-o /output`
    * define the output directory. LoDEI generates many automatically named files.
* `-c 3`
    * Number of used CPU cores.
* `--library SR`
    * provide the strandedness of your BAM files. Strandedness is defined as in `salmon` see: https://salmon.readthedocs.io/en/latest/library_type.html  
    Currently, the following types are supported:  
    SR = reverse stranded,  
    SF = forward stranded,  
    U = unstranded  
    ISR = paired-end, reverse stranded  
    ISF = paired-end, forward stranded,  
    If you are unsure what kind of library type (strandedness) your data is, have a look at the [FAQ](#how-do-i-infer-the-library-type-(strandedness)-of-my-data) and 
    https://github.com/rna-editing1/getlibtype
* `--min_coverage 5`
    * only consider single positions that have a coverage >= min_coverage in all samples.
* [Should I use](#should-i-use---rm_snps) `--rm_snps`?

Wait until LoDEI finishes the calculation (~1-2min) and have a look at the [output](#output).

## Output

The primary outputs are BED-format-like plaintext files containing the genomic coordinates, their differential editing signals and q-values of all windows. The first line is the header. Each subsequent line corresponds to a single window.

| Column name | Description |
| --- | ----------- |
| chrom | The name of the chromosome (e.g. chr2, 2) where the window was detected (string) |
| wstart | The starting position of the window (int) |
| wend | The stopping position of the window (int) |
| name | Contains the gene name where the window was detected or empty (string) |
| wEI | The calculated differential signal (see eq. 4 in the publication) (float) |
| strand | Defines the strand where the differential signals was detected. Either "+" or "-" (char) |
| q_value | Calculated q value of the detected wEI signal (float)|

LoDEI computes differential signals for all possible mismatch pairs. As a consequence, for each nucleotide mismatch X and Y an output file is generated according to the following scheme `/windows/windows_XY.txt` where X and Y are the nucleotide mismatches. Consequently, the file `/windows/windows_AG.txt` should be examined in case of A-to-I editing. Note, that the nucleotides mismatches X and Y refer to the 5'-3' orientation. If you are interested analyzing A-to-I editing you only need to look at the `_AG.txt` files. LoDEI properly handles the mismatch detection with respect to the used sequencing library and strand orientation of RNAs internally.

The results of all mismatches are located in the sub-directoy `/windows` in the output directory:

```
$ cd ~/example/output_test

$ tree windows/
windows
├── windows_AC.txt
├── windows_AG.txt
├── windows_AT.txt
├── windows_CA.txt
├── windows_CG.txt
├── windows_CT.txt
├── windows_GA.txt
├── windows_GC.txt
├── windows_GT.txt
├── windows_TA.txt
├── windows_TC.txt
└── windows_TG.txt

0 directories, 12 files
```

If windows achieve a q value < 0.1, LoDEI creates additional output files for each mismatch pair for windows with a q value < 0.1 according to the naming scheme `windows_qfiltered_XY.txt`, where X and Y are the nucleotide mismatches.


## Getting Started

Since LoDEI requires sorted BAM files as input the following steps/programs are typically performed/run prior to running LoDEI:

* `fastqc` / `multiqc` for quality control of the data
* `cutadapt` for quality filtering of the reads. The `-q` parameter specifies the quality filtering. A value of at least 20 (`-q 20`) is recommended.
* `STAR` for aligning RNA-seq data to the reference
* `samtools` for sorting and indexing the BAM files obtained from `STAR`.
* `lodei` for differential RNA editing analysis.

Keep in mind to set the `--library` parameter of `lodei` properly. 
If you are unsure what kind of library type (strandedness) your data is, have a look at the [FAQ](#how-do-i-infer-the-library-type-(strandedness)-of-my-data) and https://github.com/rna-editing1/getlibtype


## FAQ

### What kind of annotation file do I need to use?

The annotation file provided in GFF format that LoDEI takes as input should contain the genomic regions of interest for your analysis question. 
LoDEI uses a sliding window approach. 
For a given genomic region (that's an entry/line in your GFF), LoDEI calculates the differential editing for all windows that fit into that given region. 
Thus, your GFF file should fulfill the following requirements: 

* Make sure your GFF annotation file does not contain redundant entries. 
The genomic locations in your file should be unique and not overlapping with each other. 
Common genomic annotation files like basic gene annotation files obtained from gencodegenes.org should not be used without prior filtering since annotation files typically contain many redundant genomic locations since they contain genes, transcripts, and exons. A starting point for standard RNA-seq might be the set of protein-coding genes (but it depends on your experiment):
```
grep -P  "gene\t.*gene_type=protein_coding" gencode.v47.basic.annotation.gff3 > gencode.v47.basic.annotation_genes_protein_coding.gff3
```

* Ensure that the annotation file you provide to LoDEI covers a large set of genomic locations to ensure that LoDEI gets enough data to calculate q values.

### Should I use --rm_snps?

Short answer: if you are unsure, yes. 
Long answer: If you compare datasets from the same cell line you typically don't need that option. If the sets that you compare against each other contain sequencing data from different cells/samples/patients/etc. you should use this option. 

### How do I infer the library type (strandedness) of my data?

To run LoDEI it is required to specify the library type for your sequencing data via the `--library` parameter. 
We provide the additional small program `getlibtype` here https://github.com/rna-editing1/getlibtype to help you identifying your library type.
LoDEI uses the same library type specification as Salmon (https://salmon.readthedocs.io/en/latest/library_type.html). 
`getlibtype` is a small wrapper for `salmon` that utilizes `salmon` only for the library type detection. 

### Why is the library type so important?

Detecting RNA editing is based on scanning for mismatches between the sequencing data and the reference genome. 
In case of A-to-I editing, publications typically refer to the analysis of A/G mismatches between the reads and the reference. This description is correct from the perspective of 5'-3' transcript orientation and transcripts that originate from genes from the forward strand.
Unfortunately, transcripts can be located on the forward or reverse strand of the genome and the used sequencing chemistry has an impact on the type of mismatches with respect to the relative orientation of transcripts. 
In other words, the type of mismatch to look at is dependent on the location of the gene (forward or reverse) and the underlying sequencing chemistry. 
To ease the analysis and not getting down into this rabbit hole, LoDEI takes care of all of this internally.






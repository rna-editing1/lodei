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
* Podman/Docker *or* conda/mamba

To run LoDEI we recommend to use the provided [Podman/Docker](#installation-and-usage-via-podman) image or to manually install LoDEI into a [conda/mamba](#installation-and-usage-via-conda-and-pip) environment.

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

## Installation

Install LoDEI by using one of the following ways:

1. use the provided [Podman/Docker](#installation-and-usage-via-podman) image.
2. use the [conda/mamba](#installation-and-usage-via-conda-and-pip) package manager to install dependencies and LoDEI manually (LoDEI cannot be installed via `conda install` yet).

Independent the way you choose, we recommend to follow the instructions to run LoDEI on the provided [test data](#test-data) to familiarize with the usage and verify a proper installation.

### Installation and Usage via Podman

Common Linux distributions are typically shipped with Podman. Podman is a tool to create, run and maintain containers.
For a detailed introduction of Podman we refer the reader to the primary documentation at https://podman.io.

Get the container image from DockerHub (https://hub.docker.com/r/lodei/lodei) via 

```
podman pull docker.io/lodei/lodei:latest
```

Verify that the new image is part of your container image storage. You should find an entry similar to the example shown below:

```
podman images
REPOSITORY                     TAG         IMAGE ID      CREATED         SIZE
docker.io/lodei/lodei          latest      ea6601c991f9  42 minutes ago  2.3 GB
```

Verify that podman is able to start LoDEI by trying to run the new container:

```
podman run --rm docker.io/lodei/lodei:latest lodei find -h
```

You should see the help page of LoDEI:

```
usage: lodei find [-h] -a GROUP1 [GROUP1 ...] -b GROUP2 [GROUP2 ...] -f FASTA -g GFF -o OUTPUT [-c CORES] [-i ID] [-w WINDOW_SIZE] [-s STEP_SIZE] [-l LIBRARY] [-m MIN_COVERAGE] [-p] [--self] [--subprocessid SUBPROCESSID] [-v]

'lodei find' detects differentially edited regions between two groups of samples for a given set of genomic locations.

optional arguments:
  -h, --help            show this help message and exit
  -a GROUP1 [GROUP1 ...], --group1 GROUP1 [GROUP1 ...]
                        Input BAM files for group 1 separated by white space.
  -b GROUP2 [GROUP2 ...], --group2 GROUP2 [GROUP2 ...]
                        Input BAM files for group 2 separated by white space.
  -f FASTA, --fasta FASTA
                        Fasta file. Must be the same fasta file used for producing BAM 
                        files of group1 and group2.
  -g GFF, --gff GFF     Annotation in GFF file format. The program tries to find local editing
                        regions within the given GFF annotations.
  -o OUTPUT, --output OUTPUT
                        Output directory
  -c CORES, --cores CORES
                        Number of cores used. [Default: 1]
  -i ID, --id ID        Optional run ID. [Default: 1]
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        Half-window size for editing index calculation. The complete 
                        window size is: 2*w+1 [Default: 50]
  -s STEP_SIZE, --step_size STEP_SIZE
                        The editing index windows will be shifted by -s positions. [Default: 15]
  -l LIBRARY, --library LIBRARY
                        Specifies the sequencing library type. [Default: SR]
  -m MIN_COVERAGE, --min_coverage MIN_COVERAGE
                        Minimum required read coverage at regions. Only positions 
                        with a coverage >= m in all samples of all groups are 
                        used for editing index calculation. [Default: 5]
  -p, --rm_snps         Simple heuristic to remove possible SNPs [Default: off]
  --self                Only used if the the program calls itself 
                        for multi processing reasons. [Default: False]
  --subprocessid SUBPROCESSID
                        For creating a separate log file when using -c > 1. [Default: 0]
  -v, --verbose         Verbose mode. [Default: off]
```

#### Run LoDEI Using Podman

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

##### Run LoDEI

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
docker.io/lodei/lodei:latest lodei find \
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
    * provide the list of sorted BAM files (separated by space) of samples belonging to group 1.
* `--group2 ...`
    * provide the list of sorted BAM files (separated by space) of samples belonging to group 2.
* `-f /annotation/genome.fa`
    * provide the reference genome used to generate the provided BAM files.
* `/annotation/test_anno.gff3`
    * LoDEI caculates differential editing for all entries of the provided annotation file.
* `-o /output`
    * define the output directory. LoDEI generates many automatically named files.
* `-c 3`
    * Number of used CPU cores.
* `--library SR`
    * provide the strandedness of your BAM files. SR = reverse stranded, SF = forward stranded, U = unstranded. Note, currently LoDEI only handles stranded single-end or unstranded RNA-seq data. To run stranded paired-end samples you need to generate BAM files for `_1` and `_2` reads seperately and run LoDEI for the resulting BAM files individually and merge the final results.
* `--min_coverage 5`
    * only consider single positions that have a coverage >= min_coverage in all samples.

Wait until LoDEI finishes the calculation (~1-2min) and have a look at the [output](#output).

### Installation and Usage via conda and pip

Currently, LoDEI is not directly available via `conda`.
Here, a new `conda` environment is created and requirements are installed.
Next, LoDEI is downloaded and installed into the newly created environment.

#### Create Environment and Install Dependencies

Generate a new environment and install required packages:

```
conda create -c conda-forge -c bioconda --name lodei
conda activate lodei
conda install -y python=3.8 pysamstats=1.1.2 pandas=2.0.3 matplotlib=3.7.1
```

Download LoDEI and install it via `pip`.

```
cd ~
git clone https://github.com/rna-editing1/lodei.git
cd lodei
# Make sure that the pyproject.toml of LoDEI is in your current directory
# and that the correct conda environment is active
python -m pip install . --no-deps --ignore-installed
```

#### Run LoDEI Using Conda

To verify a proper installation we run LoDEI on the provided [test dataset](#test-data).

Make sure to get back into your example directory where your unpacked the testdata, create a new output directory at `~/example` where LoDEI can save the results and finally move into the example data directory:

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
  * provide the list of sorted BAM files (separated by space) of samples belonging to group 1.
* `--group2 ...`
  * provide the list of sorted BAM files (separated by space) of samples belonging to group 2.
* `-f annotation/genome.fa`
  * provide the reference genome used to generate the provided BAM files.
* `annotation/test_anno.gff3`
  * LoDEI caculates differential editing for all entries of the provided annotation file.
* `-o ../output_conda`
  * define the output directory. LoDEI generates many automatically named files.
* `-c 3`
  * Number of used CPU cores.
* `--library SR`
  * provide the strandedness of your BAM files. SR = reverse stranded, SF = forward stranded, U = unstranded. Note, currently LoDEI only handles stranded single-end or unstranded RNA-seq data. To run stranded paired-end samples you need to generate BAM files for `_1` and `_2` reads seperately and run LoDEI for the resulting BAM files individually and merge the final results.
* `--min_coverage 5`
  * only consider single positions that have a coverage >= min_coverage in all samples.

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

LoDEI computes differential signals for all possible mismatch pairs. As a consequence, for each nucleotide mismatch X and Y an output file is generated according to the following scheme `/windows/windows_XY.txt` where X and Y are the nucleotide mismatches. Consequently, the file `/windows/windows_AG.txt` should be examined in case of A-to-I editing.

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

The output for the provided test dataset is available at Zenodo: https://zenodo.org/doi/10.5281/zenodo.10907019




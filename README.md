# Differential Peak Enrichment Analysis

The repository contains utility scripts to find Differentially Enriched Regions (DER) of histone modification peaks, and a DockerFile with the directives on how to produce a container, which is also available on [DockerHub](https://hub.docker.com/r/alessandro201/differential-peak-analysis).

- `diffbind_analysis.R` uses [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) to perform the differential analysis. It then saves relevant data and statistics, and produces plots.
- `homer2gtf.py` joins the peaks annotated with [HOMER](http://homer.ucsd.edu/homer/ngs/annotation.html) with the metadata produced by DiffBind (like FDR, pvalue, and log fold change) in a TSV file, and conveniently prepares a GTF file of the significant peaks (loadable as a track in [IGV](https://igv.org/), for example).

## Table Of Contents
- [Differential Peak Enrichment Analysis](#differential-peak-enrichment-analysis)
  * [Workflow](#workflow)
  * [Samplesheet preparation](#samplesheet-preparation)
  * [Parameters](#parameters)
  * [Installation](#installation)
  * [Detailed Workflow](#detailed-workflow)
- [Manual installation](#manual-installation)
  * [Container building](#container-building)
  * [Conda environment](#conda-environment)
  * [Dependencies](#dependencies)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>



## Workflow

1. Perform differential analysis of the samples in `samplesheet.csv` using the TISSUE `control` as the base condition and save the significant peaks and the plots in `diffbind_results/`.

   ```bash
   diffbind_analysis.R samplesheet.csv -o DER -t control
   ```

2. Annotate the peaks with HOMER using the GRCh38 gencode_v47 annotated scaffold as reference (which can be downloaded [here](https://www.gencodegenes.org/human/)). You can choose any GTF file as reference. For more information on the parameters refer to the official documentation of [HOMER](http://homer.ucsd.edu/homer/ngs/annotation.html).

   ```bash
   annotatePeaks.pl DER/DER_control.tsv hg38 -gtf gencode.v47.chr_patch_hapl_scaff.annotation.gtf.gz > DER_control_annotated.tsv
   ```

3. Join the two outputs

   ```bash
   homer2gtf.py -d diffbind_results/DER_control.tsv DER_control_annotated.tsv -o DER_annotated/
   ```

## Samplesheet preparation

DiffBind requires a `samplesheet.csv` containing the `.bam` alignments and the peaks called of each sample.
The example below uses the output from the [nf-core/cutandrun](https://nf-co.re/cutandrun/3.2.2/) pipeline:

```csv
SampleID  ,Tissue  ,Factor   ,Condition ,Treatment  ,Replicate ,bamReads                                                                                  ,Peaks                                                                                    ,PeakCaller
treated_1 ,tissue1 ,H3K4_me1 ,treated   ,Full-Media ,1         ,cutandrun/results/02_alignment/bowtie2/target/markdup/treated_1.target.markdup.sorted.bam ,cutandrun/results/03_peak_calling/04_called_peaks/macs2/treated_1.macs2.peaks.cut.bed.gz ,bed
treated_2 ,tissue1 ,H3K4_me1 ,treated   ,Full-Media ,2         ,cutandrun/results/02_alignment/bowtie2/target/markdup/treated_2.target.markdup.sorted.bam ,cutandrun/results/03_peak_calling/04_called_peaks/macs2/treated_2.macs2.peaks.cut.bed.gz ,bed
treated_3 ,tissue1 ,H3K4_me1 ,treated   ,Full-Media ,3         ,cutandrun/results/02_alignment/bowtie2/target/markdup/treated_3.target.markdup.sorted.bam ,cutandrun/results/03_peak_calling/04_called_peaks/macs2/treated_3.macs2.peaks.cut.bed.gz ,bed
control_1 ,tissue2 ,H3K4_me1 ,control   ,Full-Media ,1         ,cutandrun/results/02_alignment/bowtie2/target/markdup/control_1.target.markdup.sorted.bam ,cutandrun/results/03_peak_calling/04_called_peaks/macs2/control_1.macs2.peaks.cut.bed.gz ,bed
control_2 ,tissue2 ,H3K4_me1 ,control   ,Full-Media ,2         ,cutandrun/results/02_alignment/bowtie2/target/markdup/control_2.target.markdup.sorted.bam ,cutandrun/results/03_peak_calling/04_called_peaks/macs2/control_2.macs2.peaks.cut.bed.gz ,bed
control_3 ,tissue2 ,H3K4_me1 ,control   ,Full-Media ,3         ,cutandrun/results/02_alignment/bowtie2/target/markdup/control_3.target.markdup.sorted.bam ,cutandrun/results/03_peak_calling/04_called_peaks/macs2/control_3.macs2.peaks.cut.bed.gz ,bed
```

If you wish to use consensus samples, they need to be specified in a different samplesheet that follows the same format as the above. Be sure to use the same SampleID for the samples and the corresponding consensus.

## Parameters

These are the parameters of `diffbind_analysis.R`:

- `-o` is the output directory in which the significant peaks and the plots will be saved
- `-t` defines the `TISSUE` that will be used as base condition, eg. the denominator in the fold change. The default is `BULK`
- `-c` specify a samplesheet containing the consensus peakset. If given, diffbind will avoid computing its own consensus peakset.
- `-b` specify the blacklist to use:
   - `TRUE`: DEFAULT - Automatically inferr the blacklist to use
   - `DBA_BLACKLIST_HG19`: Homo sapiens 19 (chromosomes have 'chr')
   - `DBA_BLACKLIST_HG38`: Homo sapiens 38 (chromosomes have 'chr')
   - `DBA_BLACKLIST_GRCH37`: Homo sapiens 37 (chromosomes are numbers)
   - `DBA_BLACKLIST_GRCH38`: Homo sapiens 38 (chromosomes are numbers)
   - `DBA_BLACKLIST_MM9`: Mus musculus 9
   - `DBA_BLACKLIST_MM10`: Mus musculus 10
   - `DBA_BLACKLIST_CE10`: C. elegans 10
   - `DBA_BLACKLIST_CE11`: C. elegans 11
   - `DBA_BLACKLIST_DM3`: Drosophila melanogaster 3
   - `DBA_BLACKLIST_DM6`: Drosophila melanogaster 6
- `-m` Defines the minimum number of replicates a peak must be in to be considered consensus. Used by diffbind to find a consensus peakset if not given
- `-s` defines how the summit will be computed. DiffBind by default shortens the peaks around ±200bp from the summit, hence 401 total. Depending on the histon marks, it may be appropiate or not. Choices:
  - `-s INTEGER` choose a specific distance from the summit
  - `-s false` disable the centering of the peaks
  - `-s media` use the median peak length. **_BEWARE_** that the median is computed on all samples, thus if you have multiple histon marks and you want to use the median peak length, you should use a samplesheet per histon mark.

```text
$ diffbind_analysis.R -h
usage: diffbind_analysis.R [--] [--help] [--opts OPTS] [--outdir
       OUTDIR] [--tissue-contrast TISSUE-CONTRAST] [--summit SUMMIT]
       [--blacklist BLACKLIST] [--minOverlap MINOVERLAP] [--consensus
       CONSENSUS] samplesheet

Launch DiffBind on the data and perform the plots

positional arguments:
  samplesheet            Samplesheet in csv format

flags:
  -h, --help             show this help message and exit

optional arguments:
  -x, --opts             RDS file containing argument values
  -o, --outdir           Output directory [default: .]
  -t, --tissue-contrast  Set the baseline condition, eg. the
                         denominator in the fold change. You should
                         have two distinct tissues, and you should
                         specify the one you want as baseline
                         [default: BULK]
  -s, --summit           Re-center each peak interval around its point
                         of highest pileup.  '--summit 200' will select
                         -200/+200 bp around the point of highest
                         pileup giving peaks of 401bp.  '--summit
                         false' will disable the re-centering.
                         '--summit median' will create peaks about the
                         median peak size of ALL SAMPLES. THERE WILL BE
                         NO DISTINCTION BETWEEN CONDITION, TISSUE OR
                         ANY OTHER METRIC.  Choices: ['false',
                         'median', INTEGER] [default: 200]
  -b, --blacklist        Apply the ENCODE blacklist to the peaks.
                         Choices: 'TRUE': automatically infer the
                         genome and use the corresponding ENCODE
                         blacklist 'FALSE': does not apply a blacklist
                         'DBA_BLACKLIST_HG19': Homo sapiens 19
                         (chromosomes have 'chr') 'DBA_BLACKLIST_HG38':
                         Homo sapiens 38 (chromosomes have 'chr')
                         'DBA_BLACKLIST_GRCH37': Homo sapiens 37
                         (chromosomes are numbers)
                         'DBA_BLACKLIST_GRCH38': Homo sapiens 38
                         (chromosomes are numbers) 'DBA_BLACKLIST_MM9':
                         Mus musculus 9 'DBA_BLACKLIST_MM10': Mus
                         musculus 10 'DBA_BLACKLIST_CE10': C. elegans
                         10 'DBA_BLACKLIST_CE11': C. elegans 11
                         'DBA_BLACKLIST_DM3': Drosophila melanogaster 3
                         'DBA_BLACKLIST_DM6': Drosophila melanogaster 6
                         [default: TRUE]
  -m, --minOverlap       Set the minimum number of replicates a peak
                         must be in to be considered consensus.  If
                         it's between 0 and 1, peaks will be included
                         if they are present in at least this
                         proportion of replicates.  [default: 2]
  -c, --consensus        Specify a samplesheet containing the consensus
                         peakset for each of the Tissue in the
                         Samplesheet. It will be used instead of the
                         one generated by diffbind. Incompatible with
                         --minOverlap.
```

```text
$ homer2gtf.py -h
usage: homer2gtf [-h] [-d [DIFFBIND]] [-f FOLD] [-r FDR] [-p PVALUE] [-t TSS] [--type [TYPE ...]] [--no-type [NO_TYPE ...]] [-o [OUTPUT_PREFIX]] FILE

Transform the output of HOMER annotation to a bed file readable by IGV

positional arguments:
  FILE                  Output from HOMER annotation as .tsv file

options:
  -h, --help            show this help message and exit
  -d [DIFFBIND], --diffbind [DIFFBIND]
                        Table of differentially significant peaks given to HOMER. This will be used to get the differentially bound strength
  -o [OUTPUT_PREFIX], --output-prefix [OUTPUT_PREFIX]
                        Prefix of the output. Default: [FILE]_joined

filters:
  -f FOLD, --fold FOLD  Keep only peaks with abs(log10 fold change) >= FOLD
  -r FDR, --fdr FDR     Keep only peaks with -log10(fdr) >= FDR
  -p PVALUE, --pvalue PVALUE
                        Keep only peaks with -log10(pvalue) >= PVALUE
  -t TSS, --tss TSS     Keep only peaks with abs(distance from tss) <= TSS
  --type [TYPE ...]     Keep only peaks that have type TYPE ['exon', 'intron', 'intergenic', ...]. To add multiple types divide them with commas
  --no-type [NO_TYPE ...]
                        Filter out peaks that have type TYPE ['exon', 'intron', 'intergenic', ...]. To add multiple types divide them with commas
```

## Installation

You can find a docker container with the software already installed on [DockerHub](https://hub.docker.com/r/alessandro201/differential-peak-analysis). To download it and run an interactive shell do:

```bash
docker pull alessandro201/differential-peak-analysis
docker run --rm -it alessandro201/differential-peak-analysis bash

# OR
singularity pull docker://alessandro201/differential-peak-analysis
singularity run differential-peak-analysis_latest.sif bash
```

Singularity will download the image from docker and convert it into a `.sif ` file in the current directory.

Clone the repo to have the latest version of the script available:
```
git clone https://github.com/Alessandro201/differential-peak-analysis.git
```

## Detailed Workflow
To perform the analysis on the current directory, using the samplesheet from the example above, do:

```bash
# CURRENT_DIRECTORY/
# ├── differential-peak-analysis               # [THIS_REPO]
# │       ├── homer2gtf.py
# │       ├── diffbind_analysis.R
# │       └── ...
# ├── cutandrun/
# │   └── results/
# │       ├── 02_alignment/
# │       │   └── ...
# │       └── 03_peak_calling/
# │           └── ...
# ├── samplesheet.csv
# ├── consensus_samplesheet.csv
# └── gencode.v4y.chr_patch_hapl_scaff.annotation.gtf.gz

docker run --rm -it -v .:/mnt -w /mnt alessandro201/differential-peak-analysis ./differential-peak-analysis/diffbind_analysis.R samplesheet.csv -o DER/ -t control -c consensus_samplesheet.csv --summit median

mkdir DER_CONTROL_HOMER/
docker run --rm -it -v .:/mnt -w /mnt alessandro201/differential-peak-analysis bash -c 'annotatePeaks.pl DER/DER_control.bed hg38 -gtf gencode.v47.chr_patch_hapl_scaff.annotation.gtf.gz -annStats DER_CONTROL_HOMER/annotations_stats.txt -go DER_CONTROL_HOMER/go_annotations -genomeOntology DER_CONTROL_HOMER/genome_ontology -cpu 4 > DER_CONTROL_HOMER/DER_control_annotated.tsv'

mkdir DER_CONTROL_FINAL/
docker run --rm -it -v .:/mnt -w /mnt alessandro201/differential-peak-analysis ./differential-peak-analysis/homer2gtf.py DER_CONTROL_HOMER/DER_control_annotated.tsv --diffbind DER/DER_control.bed --output DER_CONTROL_FINAL/ --tss 100000 --fdr 0.001
```

Or with singularity
```bash
# CURRENT_DIRECTORY/
# ├── differential-peak-analysis_latest.sif    # [singularity image]
# ├── differential-peak-analysis               # [THIS_REPO]
# │       ├── homer2gtf.py
# │       ├── diffbind_analysis.R
# │       └── ...
# ├── cutandrun/
# │   └── results/
# │       ├── 02_alignment/
# │       │   └── ...
# │       └── 03_peak_calling/
# │           └── ...
# ├── samplesheet.csv
# ├── consensus_samplesheet.csv
# └── gencode.v47.chr_patch_hapl_scaff.annotation.gtf.gz

singularity run -B .:/mnt -W /mnt ./differential-peak-analysis_latest.sif ./differential-peak-analysis/diffbind_analysis.R samplesheet.csv -o DER/ -t control -c consensus_samplesheet.csv --summit median

mkdir DER_CONTROL_HOMER/
singularity run -B .:/mnt -W /mnt ./differential-peak-analysis_latest.sif -c 'annotatePeaks.pl DER/DER_control.bed hg38 -gtf gencode.v47.chr_patch_hapl_scaff.annotation.gtf.gz -annStats DER_CONTROL_HOMER/annotations_stats.txt -go DER_CONTROL_HOMER/go_annotations -genomeOntology DER_CONTROL_HOMER/genome_ontology -cpu 4 > DER_CONTROL_HOMER/DER_control_annotated.tsv'

mkdir DER_CONTROL_FINAL/
singularity run -B .:/mnt -W /mnt ./differential-peak-analysis_latest.sif ./differential-peak-analysis/homer2gtf.py DER_CONTROL_HOMER/DER_control_annotated.tsv --diffbind DER/DER_control.bed --output DER_CONTROL_FINAL/ --tss 100000 --fdr 0.001
```

> [!NOTE]
> `annotatePeaks.pl` outputs the file to STDOUT, hence we need to redirect it to a file using `>`. However, while all the commands are passed to the container, the redirection sign and what comea after is excluded. Thus we need to pass the whole command as a string to bash so that it can be executed fully from within the container as expected.


# Manual installation
## Container building

To build the container clone this repository and build it with:

```bash
git clone https://github.com/Alessandro201/differential-peak-analysis.git
cd differential-peak-analysis
docker build -t IIT/differential-peak-analysis . --platform=linux/amd64
```

## Conda environment
If you don't want to create or use a docker container you can install the environment with conda (or mamba or micromamba):
```bash
git clone https://github.com/Alessandro201/differential-peak-analysis.git
cd differential-peak-analysis
conda env create -f env.yml
conda activate differential-analysis

# Install other necessary R dependencies:
Rscript install_dependencies.R

# Activate the environment with conda, mamba or micromamba
micromamba activate differential-analysis

# Install the reference genome for HOMER:
$CONDA_PREFIX/share/homer/configureHomer.pl -install hg38
```

> [!NOTE]
> Depending on the package manager or how you installed HOMER, `configureHomer.pl` could be in different locations. 

## Dependencies

These are the dependencies needed by the scripts and the programs, that will be installed following the previous steps.

The R script `differential_analysis.R` depends on:

- [argparser](https://cran.r-project.org/web/packages/argparser/index.html)
- [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
- [Bioconductor](https://bioconductor.org/install/)
- [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html)
- [profileplyr](https://bioconductor.org/packages/release/bioc/html/profileplyr.html)

The Python script `homer2gtf.py` depends on:

- [numpy](https://numpy.org/install/)
- [pandas](https://pandas.pydata.org/)

To use HOMER you need to:

- Install [HOMER](http://homer.ucsd.edu/homer/introduction/install.html)

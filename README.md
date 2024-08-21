# Differential Peak Analysis

This repository contains utility scripts to perform differential peak analysis.

- `diffbind_analysis.R` uses [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) to perform the differential analysis. It then saves relevant data and statistics, and produces plots.
- `homer2igv.py` joins the peaks annotated with [HOMER](http://homer.ucsd.edu/homer/ngs/annotation.html) with the metadata produced by DiffBind (like FDR, pvalue, and log fold change) in a TSV file, and conveniently prepares a GTF file of the significant peaks (loadable as a track in [IGV](https://igv.org/), for example).

## Workflow

1. Perform differential analysis of the samples in `samplesheet.csv` using `control` as the base condition and save the significant peaks and the plots in `diffbind_results/`.

   ```bash
   diffbind_analysis.R samplesheet.csv -o diffbind_results/ -c control
   ```

2. Annotate the peaks with HOMER using the GRCh38 gencode_v46 annotated scaffold as reference. You can choose any GTF file as reference. For more information on the parameters refer to the official documentation of [HOMER](http://homer.ucsd.edu/homer/ngs/annotation.html).

   ```bash
   annotatePeaks.pl diffbind_results/differentially_bound_sites.tsv hg38 -gtf gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz > differentially_bound_sites_annotated.tsv
   ```

3. Join the two outputs

   ```bash
   homer2igv.py -d diffbind_results/differentially_bound_sites.tsv differentially_bound_sites_annotated.tsv
   ```

## Samplesheet preparation

DiffBind requires a `samplesheet.csv` containing the `.bam` alignments and the peaks called of each sample.
The example below uses the output from the [nf-core/cutandrun](https://nf-co.re/cutandrun/3.2.2/) pipeline:

```csv
SampleID  ,Tissue  ,Factor   ,Condition ,Treatment  ,Replicate ,bamReads                                                                                  ,Peaks                                                                                    ,PeakCaller
treated_1 ,tissue1 ,H3k4_me1 ,treated   ,Full-Media ,1         ,cutandrun/results/02_alignment/bowtie2/target/markdup/treated_1.target.markdup.sorted.bam ,cutandrun/results/03_peak_calling/04_called_peaks/macs2/treated_1.macs2.peaks.cut.bed.gz ,bed
treated_2 ,tissue1 ,H3k4_me1 ,treated   ,Full-Media ,2         ,cutandrun/results/02_alignment/bowtie2/target/markdup/treated_2.target.markdup.sorted.bam ,cutandrun/results/03_peak_calling/04_called_peaks/macs2/treated_2.macs2.peaks.cut.bed.gz ,bed
treated_3 ,tissue1 ,H3k4_me1 ,treated   ,Full-Media ,3         ,cutandrun/results/02_alignment/bowtie2/target/markdup/treated_3.target.markdup.sorted.bam ,cutandrun/results/03_peak_calling/04_called_peaks/macs2/treated_3.macs2.peaks.cut.bed.gz ,bed
control_1 ,tissue2 ,H3k4_me1 ,control   ,Full-Media ,1         ,cutandrun/results/02_alignment/bowtie2/target/markdup/control_1.target.markdup.sorted.bam ,cutandrun/results/03_peak_calling/04_called_peaks/macs2/control_1.macs2.peaks.cut.bed.gz ,bed
control_2 ,tissue2 ,H3k4_me1 ,control   ,Full-Media ,2         ,cutandrun/results/02_alignment/bowtie2/target/markdup/control_2.target.markdup.sorted.bam ,cutandrun/results/03_peak_calling/04_called_peaks/macs2/control_2.macs2.peaks.cut.bed.gz ,bed
control_3 ,tissue2 ,H3k4_me1 ,control   ,Full-Media ,3         ,cutandrun/results/02_alignment/bowtie2/target/markdup/control_3.target.markdup.sorted.bam ,cutandrun/results/03_peak_calling/04_called_peaks/macs2/control_3.macs2.peaks.cut.bed.gz ,bed
```

## Parameters

These are the parameters of `diffbind_analysis.R`:

- `-o` is the output directory in which the significant peaks and the plots will be saved
- `-c` defined the base condition, hence the denominator in the fold change. The defualt is `control`
- `-s` defines how the summit will be computed. DiffBind by default centers the peaks keeping -200bp/+200bp from the summit, hence each peak will have 401bp total. Depending on the histon marks, it may be appropriate or not:
  - `-s INTEGER` choose a specific distance from the summit
  - `-s false` disable the centering of the peaks
  - `-s median` use the median peak length. **_BEWARE_** that the median is computed on all samples, thus if you have multiple histon marks and you want to use the median peak length, you should use a samplesheet per histon mark.

```text
$ diffbind_analysis.R -h
usage: diffbind_analysis.R [--] [--help] [--opts OPTS] [--outdir
       OUTDIR] [--contrast CONTRAST] samplesheet

Launch DiffBind on the data and perform the plots

positional arguments:
  samplesheet     Samplesheet in csv format

flags:
  -h, --help      show this help message and exit

optional arguments:
  -x, --opts      RDS file containing argument values
  -o, --outdir    Output directory [default: .]
  -c, --contrast  Set the baseline condition, eg. the denominator in
                  the fold change [default: control]
  -s, --summit    Re-center each peak interval around its point of
                  highest pileup.  '--summit 200' will select -200/+200
                  bp around the point of highest pileup giving peaks of
                  401bp.  '--summit false' will disable the
                  re-centering.  '--summit median' will create peaks
                  about the median peak size of all samples.  Choices:
                  ['false', 'median', INTEGER] [default: 200]
```

```text
$ homer2igv.py -h
usage: homer2igv [-h] [-d [DIFFBIND]] [-f FOLD] [-r FDR] [-p PVALUE] [-t TSS] [--type [TYPE ...]] [--no-type [NO_TYPE ...]] [-o [OUTPUT_PREFIX]] FILE

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

## Software installation

You can find a docker container with the software already installed at [CONTAINER POSITION]. To download it and run an interactive shell do:

```bash
docker pull [CONTAINER POSITION]
docker run --rm -it [CONTAINER NAME] bash
```

To perform the analysis on the current directory, using the samplesheet from the example above, do:

```text
CURRENT_DIRECTORY/
├── cutandrun/
│   └── results/
│       ├── 02_alignment/
│       │   └── ...
│       └── 03_peak_calling/
│           └── ...
├── samplesheet.csv
└── gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz
```

```bash
docker run --rm -it -v .:/mnt -w /mnt [CONTAINER NAME] diffbind_analysis.R samplesheet.csv -o diffbind_results/ -c control
docker run --rm -it -v .:/mnt -w /mnt [CONTAINER NAME] bash -c 'annotatePeaks.pl diffbind_results/differentially_bound_sites.tsv hg38 -gtf gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz > differentially_bound_sites_annotated.tsv'
docker run --rm -it -v .:/mnt -w /mnt [CONTAINER NAME] homer2igv.py -d diffbind_results/differentially_bound_sites.tsv differentially_bound_sites_annotated.tsv
```

> [!NOTE] > `annotatePeaks.pl` outputs the file to STDOUT, hence we need to redirect it to a file using `>`. However, while all the commands are passed to the container, the redirection sign is excluded and puts all the output of the command executed by docker (in this case `annotatePeaks.pl`) in the specified file. Given that in the example above we mounted the current directory on the container because we would like to keep all the outputs, it doesn't change the final output. However if we were to mount a directory which is not the current one, the output would be saved in the current directory and not in the mounted one. This is why in the example above the command is passed as a string to bash, to be executed inside the container as expected.

### Container building

To build the container clone this repository and build it with:

```bash
git clone https://github.com/Alessandro201/differential-peak-analysis.git
cd differential-peak-analysis
docker build -t IIT/differential-peak-analysis . --platform=linux/amd64
```

### PIXI (Conda) environment

If you don't want to build a docker container, the recommended way to install all the dependencies is to use the [PIXI](https://github.com/prefix-dev/pixi) package manager, which is built on the conda ecosystem but is faster and better than conda.

```bash
git clone https://github.com/Alessandro201/differential-peak-analysis.git
cd differential-peak-analysis
pixi shell
```

Alternatively with conda:

```bash
git clone https://github.com/Alessandro201/differential-peak-analysis.git
cd differential-peak-analysis
conda env create -f env.yml
conda activate differential-analysis
```

Install other necessary R dependencies:
```bash
Rscript install_dependencies.R
```

Install the reference genome for HOMER:

```bash
.pixi/envs/default/share/homer/configureHomer.pl -install hg38
```

> [!NOTE]
> Depending on the package manager or how you installed HOMER, `configureHomer.pl` could be in different locations.

### Dependencies

These are the dependencies needed by the scripts and the programs, that will be installed following the previous steps.

The R script `differential_analysis.R` depends on:

- [argparser](https://cran.r-project.org/web/packages/argparser/index.html)
- [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
- [Bioconductor](https://bioconductor.org/install/)
- [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html)
- [profileplyr](https://bioconductor.org/packages/release/bioc/html/profileplyr.html)

The Python script `homer2igv.py` depends on:

- [numpy](https://numpy.org/install/)
- [pandas](https://pandas.pydata.org/)

To use HOMER you need to:

- Install [HOMER](http://homer.ucsd.edu/homer/introduction/install.html)
- Install the reference genome for HOMER with:
  ```bash
  configureHomer.pl -install hg38
  ```

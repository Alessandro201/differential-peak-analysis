# Differential Peak Analysis

This repository contains utility scripts to perform differential peak analysis.

- `diffbind_analysis.R` uses [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) to perform the differential analysis. It then saves relevant data and statistics, and produces plots.
- `homer2igv.py` joins the peaks annotated with [HOMER](http://homer.ucsd.edu/homer/ngs/annotation.html) with the metadata produced by DiffBind (like FDR, pvalue, and log fold change) in a TSV file, and conveniently prepares a GTF file of the significant peaks (loadable as a track in [IGV](https://igv.org/), for example).


## Workflow

1) Perform differential analysis using `control` as the base condition.
   ```fish
   ./diffbind_analysis.R samplesheet.csv -o diffbind_results/ -c control
   ```

1) Annotate the peaks with HOMER using the GRCh38 gencodev46 annotated scaffold as reference.
   ```fish
   annotatePeaks.pl diffbind_results/differentially_bound_sites.tsv hg38 -gtf gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz > differentially_bound_sites_annotated.tsv
   ```

2) Join the two outputs
   ```fish
   ./homer2igv.py -d diffbind_results/differentially_bound_sites.tsv differentially_bound_sites_annotated.tsv
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

```
$ ./diffbind_analysis.R -h
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
```

```
$ ./homer2igv.py -h
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

The recommended way to install all the dependencies is using [pixi](https://github.com/prefix-dev/pixi) or CONDA.

```bash
git clone https://github.com/Alessandro201/differential-peak-analysis.git
cd differential-peak-analysis
pixi shell
```
> [!IMPORTANT]
> Remember to call `pixi shell` from the project folder to activate the virtual environment and let the script access the dependencies

To install the reference genome for HOMER you need to run:

```bash
.pixi/envs/default/share/homer/configureHomer.pl -install hg38
```
> [!NOTE]
> Depending on how you installed HOMER, `configureHomer.pl` could be in different locations.



### Dependencies

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
  ```fish
  configureHomer.pl -install hg38
  ```

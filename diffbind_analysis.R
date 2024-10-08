#!/usr/bin/env Rscript

suppressMessages(library(argparser, quietly = TRUE))
suppressMessages(library(stringr, quietly = TRUE))

# Create a parser
p <- arg_parser("Launch DiffBind on the data and perform the plots")

# Add command line arguments
p <- add_argument(p, "samplesheet", help = "Samplesheet in csv format", type = "character", nargs = 1)
p <- add_argument(p, "--outdir", help = "Output directory", default = ".", type = "character", nargs = 1)
p <- add_argument(p, "--contrast", help = "Set the baseline condition, eg. the denominator in the fold change", default = "control", type = "character", nargs = 1)
p <- add_argument(p, "--summit", help = "Re-center each peak interval around its point of highest pileup.
    '--summit 200' will select -200/+200 bp around the point of highest pileup giving peaks of 401bp.
    '--summit false' will disable the re-centering.
    '--summit median' will create peaks about the median peak size of all samples.
    Choices: ['false', 'median', INTEGER]", default = "200", type = "character", nargs = 1)
p <- add_argument(p, "--blacklist",
    help = "Apply the ENCODE blacklist to the peaks. Choices:
    'TRUE': automatically infer the genome and use the corresponding ENCODE blacklist
    'FALSE': does not apply a blacklist
    'DBA_BLACKLIST_HG19': Homo sapiens 19 (chromosomes have 'chr')
    'DBA_BLACKLIST_HG38': Homo sapiens 38 (chromosomes have 'chr')
    'DBA_BLACKLIST_GRCH37': Homo sapiens 37 (chromosomes are numbers)
    'DBA_BLACKLIST_GRCH38': Homo sapiens 38 (chromosomes are numbers)
    'DBA_BLACKLIST_MM9': Mus musculus 9
    'DBA_BLACKLIST_MM10': Mus musculus 10
    'DBA_BLACKLIST_CE10': C. elegans 10
    'DBA_BLACKLIST_CE11': C. elegans 11
    'DBA_BLACKLIST_DM3': Drosophila melanogaster 3
    'DBA_BLACKLIST_DM6': Drosophila melanogaster 6 ",
    default = "TRUE", type = "character", nargs = 1
)
p <- add_argument(p, "--minOverlap",
    help = "Set the minimum number of replicates a peak must be in to be considered consensus.
    If it's between 0 and 1, peaks will be included if they are present in at least this proportion of replicates. ",
    default = 0.66, type = "number", nargs = 1
)

# Parse the command line arguments
args <- parse_args(p)

# Check that the control given is present in the samplesheet
sample_sheet <- read.csv(args$samplesheet)
sample_sheet$Condition <- str_trim(sample_sheet$Condition)
if (!args$contrast %in% sample_sheet$Condition) {
    print("The contrast you have inserted is not present in the 'condition' column of the samplesheet")
    stop()
}

# Validate the summit given
contains_only_numbers <- function(x) !grepl("\\D", x)
summit <- str_to_lower(args$summit)
if ((!summit %in% c("false", "median")) && (!contains_only_numbers(summit))) {
    print(paste0("The --summit parameter given is not a valid option. Available 'false', 'median' or a whole number, given: ", summit))
    stop()
}

sample_sheet$bamReads <- str_trim(sample_sheet$bamReads)
sample_sheet$Peaks <- str_trim(sample_sheet$Peaks)

# Check if the BAM and peak files exist
for (bam in sample_sheet$bamReads) {
    if (!file.exists(bam)) {
        cat("The following BAM files does not exists:", bam, "\n")
        cat("Current working directory:", getwd(), "\n")
        quit()
    }
}

for (peak in sample_sheet$Peaks) {
    if (!file.exists(peak)) {
        cat("The following Peak files does not exists:", peak, "\n")
        cat("Current working directory:", getwd(), "\n")
        quit()
    }
}

# Remove one trailing slashes and create outdir
args$outdir <- gsub("/$", "", args$outdir)
args$outdir <- gsub("\\$", "", args$outdir)
dir.create(args$outdir, showWarnings = FALSE)

# Load heavy libraries later
print("Loading libraries, it may take a while")
suppressMessages(library(DiffBind, quietly = TRUE))
suppressMessages(library(profileplyr, quietly = TRUE))
suppressMessages(library(data.table, quietly = TRUE))


# Check BLACKLIST argument
if (str_to_lower(args$blacklist) == "true") {
    blacklist <- TRUE
} else if (str_to_lower(args$blacklist) == "false") {
    blacklist <- FALSE
} else if (str_to_lower(args$blacklist) == "dba_blacklist_hg19") {
    blacklist <- DBA_BLACKLIST_HG19
} else if (str_to_lower(args$blacklist) == "dba_blacklist_hg38") {
    blacklist <- DBA_BLACKLIST_HG38
} else if (str_to_lower(args$blacklist) == "dba_blacklist_grch37") {
    blacklist <- DBA_BLACKLIST_GRCH37
} else if (str_to_lower(args$blacklist) == "dba_blacklist_grch38") {
    blacklist <- DBA_BLACKLIST_GRCH38
} else if (str_to_lower(args$blacklist) == "dba_blacklist_mm9") {
    blacklist <- DBA_BLACKLIST_MM9
} else if (str_to_lower(args$blacklist) == "dba_blacklist_mm10") {
    blacklist <- DBA_BLACKLIST_MM10
} else if (str_to_lower(args$blacklist) == "dba_blacklist_ce10") {
    blacklist <- DBA_BLACKLIST_CE10
} else if (str_to_lower(args$blacklist) == "dba_blacklist_ce11") {
    blacklist <- DBA_BLACKLIST_CE11
} else if (str_to_lower(args$blacklist) == "dba_blacklist_dm3") {
    blacklist <- DBA_BLACKLIST_DM3
} else if (str_to_lower(args$blacklist) == "dba_blacklist_dm6") {
    blacklist <- DBA_BLACKLIST_DM6
} else {
    print(paste0("The --blacklist parameter given is not a valid option. Check the help page to know which options are available, given: ", args$blacklist))
    stop()
}


# ANALYSIS FLOW
# me1 <- dba(sampleSheet="tamoxifen.csv") %>%
#     + dba.blacklist()                         %>%
#     + dba.count()                             %>%
#     + dba.normalize()                         %>%
#     + dba.contrast()                          %>%
#     + dba.analyze()


print("Loading Samplesheet")
dba_samples <- dba(sampleSheet = sample_sheet)
dba_samples <- dba.blacklist(dba_samples, blacklist = blacklist)
dba_samples$config$greylist.pval <- 0.999


# Overlap rate of peaks between the replicates. It helps to decide how many replicate to use to define a consensus
print("Saving peaks overlap rates between the replicates as peak_overlap_rates_between_replicates.pdf")
olap.rate <- dba.overlap(dba_samples, mode = DBA_OLAP_RATE)
olap.rate
pdf(file.path(args$outdir, "peak_overlap_rates_between_replicates.pdf"))
plot(olap.rate, type = "b", ylab = "# peaks", xlab = "Overlap at least this many peaksets")
dev.off()

print("Saving correlation heatmap generated using the called peaks as correlation_hm_peaks_occupancy.pdf")
pdf(file.path(args$outdir, "correlation_hm_peaks_occupancy.pdf"))
plot(dba_samples)
dev.off()


if (summit == "false") {
    summit <- FALSE
} else if (summit == "median") {
    print("Computing the median length of all peaks")
    library(GenomicRanges)
    peak_lengths <- list()
    for (file_name in sample_sheet$Peaks) {
        file_name <- str_trim(file_name)
        bedfile <- read.table(file_name, header = FALSE, col.names = c("chr", "start", "end", "name", "score"), stringsAsFactors = FALSE)
        peaks <- makeGRangesFromDataFrame(bedfile, keep.extra.columns = TRUE)
        peaks$length <- width(peaks)
        peak_lengths <- append(peak_lengths, peaks$length)
    }
    summit <- median(peak_lengths)
    print(paste0("Peaks median length of all samples: ", summit))
} else {
    summit <- as.numeric(summit)
}

dba_consensus <- dba.peakset(dba_samples, consensus = -DBA_REPLICATE, minOverlap = args$minOverlap)
dba_consensus <- dba(dba_consensus, mask = dba_consensus$masks$Consensus, minOverlap = 1)
consensus_peaks <- dba.peakset(dba_consensus, bRetrieve = TRUE)
dba_samples <- dba.count(dba_samples, summit = summit, peaks = consensus_peaks)

# Save plot as pdf
print("Saving correlation heatmap generated using the affinity (n. of reads in consensous peaks, I think) as correlation_hm_peaks_affinity.pdf")
pdf(file.path(args$outdir, "correlation_hm_peaks_affinity.pdf"))
plot(dba_samples)
dev.off()


print("Showing the number of reads that are in consensous peaks: FRiP (Fraction of Reads in Peaks)")
info <- dba.show(dba_samples)
libsizes <- cbind(LibReads = info$Reads, FRiP = info$FRiP, PeakReads = round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
libsizes
fwrite(data.frame(libsizes, row.names = info$ID), file = file.path(args$outdir, "FRiP_per_sample.tsv"), sep = "\t", row.names = TRUE)


print("Normalizing based on sequencing depth")
dba_samples <- dba.normalize(dba_samples)


print("Modeling study design based on Condition and performing the differential analysis")
dba_samples <- dba.contrast(dba_samples, reorderMeta = list(Condition = args$contrast))
dba_samples <- dba.analyze(dba_samples)


# Save plot as pdf
print("Saving correlation heatmap generated using only significantly differentially bound sites as correlation_hm_peaks_significant.pdf")
pdf(file.path(args$outdir, "correlation_hm_peaks_significant.pdf"))
plot(dba_samples, contrast = 1)
dev.off()


print("Computing report of differentially bound sites...")
dba_samples.DB <- dba.report(dba_samples)
print(dba_samples.DB)


# Write differentially bound sites to a file in a format easy to read for HOMER annotation tool
db <- dba_samples.DB
db$PeakId <- names(db)
db <- data.frame(db)
names(db)[names(db) == "seqnames"] <- "Chr"

# Reorder columns to move PeakId at the start and have the first columns as: PeakId, Chr, start, end, strand
db <- db[, c(1, 2, 3, 12, 9, 5, 6, 7, 8, 10, 11, 4)]
# Rename the first column 'seqnames' to 'Chr'
fwrite(db, file = file.path(args$outdir, "differentially_bound_sites.tsv"), sep = "\t")

conditions <- unique(dba_samples$Condition[dba_samples$Condition != args$control])
num_db_sites <- sprintf(
    "Number of enriched sites in '%s' samples: %s",
    paste(conditions, collapse = ", "),
    sum(dba_samples.DB$Fold > 0)
)
num_db_sites <- sprintf(
    "Number of enriched sites in '%s' samples: %s",
    args$control,
    sum(dba_samples.DB$Fold < 0)
)
print(num_db_sites)
fwrite(list(num_db_sites), file = file.path(args$outdir, "num_enriched_sites.txt"), quote = FALSE)


print("Saving PCA plot of normalized read counts as pca_read_counts.pdf")
pdf(file.path(args$outdir, "pca_read_counts.pdf"))
dba.plotPCA(dba_samples, DBA_TISSUE, label = DBA_CONDITION)
dev.off()


print("Saving PCA plot of significantly differentially bound sites as pca_significant_sites.pdf")
pdf(file.path(args$outdir, "pca_significant_sites.pdf"))
dba.plotPCA(dba_samples, contrast = 1, label = DBA_TISSUE)
dev.off()


print("Saving PCA plot showing where the replicates for each of the unique tissues lies as pca_replicates_by_tissue.pdf")
pdf(file.path(args$outdir, "pca_replicates_by_tissue.pdf"))
dba.plotPCA(dba_samples, attributes = c(DBA_TISSUE, DBA_CONDITION), label = DBA_REPLICATE)
dev.off()

print("Saving MA plot of treated vs control samples as MA_treated_vs_control.pdf")
pdf(file.path(args$outdir, "MA_treated_vs_control.pdf"))
dba.plotMA(dba_samples)
dev.off()


print("Saving MA plot of treated vs control samples with concentration of each sample groups as MA_treated_vs_control_with_concentration.pdf")
pdf(file.path(args$outdir, "MA_treated_vs_control_with_concentration.pdf"))
dba.plotMA(dba_samples, bXY = TRUE)
dev.off()


print("Saving venn plot as venn_plot.pdf")
pdf(file.path(args$outdir, "venn_plot.pdf"))
dba.plotVenn(dba_samples, contrast = 1, bDB = TRUE, bGain = TRUE, bLoss = TRUE, bAll = FALSE)
dev.off()


print("Saving venn plot with all counts as venn_plot_all.pdf")
pdf(file.path(args$outdir, "venn_plot_all.pdf"))
dba.plotVenn(dba_samples, contrast = 1, bDB = TRUE, bGain = TRUE, bLoss = TRUE, bAll = TRUE)
dev.off()


print("Saving vulcano plot of log FoldChange as vulcano_plot.pdf")
pdf(file.path(args$outdir, "vulcano_plot.pdf"))
dba.plotVolcano(dba_samples)
dev.off()


print("Saving boxplot of concentration of treated vs control samples as boxplot_treated_vs_control.pdf")
pdf(file.path(args$outdir, "boxplot_treated_vs_control.pdf"))
pvals <- dba.plotBox(dba_samples)
dev.off()
fwrite(pvals, file = file.path(args$outdir, "pvals_boxplot.tsv"), sep = "\t", row.names = TRUE)


print("Saving correlation heatmap of each binding site as correlation_hm_binding_sites.pdf")
pdf(file.path(args$outdir, "correlation_hm_binding_sites.pdf"))
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
dba.plotHeatmap(dba_samples, contrast = 1, correlations = FALSE, scale = "row", colScheme = hmap)
dev.off()


print("Saving profile plot of each binding site as profile_plot_binding_sites.pdf")
pdf(file.path(args$outdir, "profile_plot_binding_sites.pdf"))
profiles <- dba.plotProfile(dba_samples)
dba.plotProfile(profiles)
dev.off()


print("Saving profile plot of each binding site of merged groups (TISSUE and REPLICATES) as profile_plot_group_binding_sites.pdf")
pdf(file.path(args$outdir, "profile_plot_group_binding_sites.pdf"))
profiles <- dba.plotProfile(dba_samples, merge = c(DBA_TISSUE, DBA_REPLICATE))
dba.plotProfile(profiles)
dev.off()


print("Saving profile plot of each binding site of each sample as profile_plot_sample_binding_sites.pdf")
pdf(file.path(args$outdir, "profile_plot_sample_binding_sites.pdf"))
profiles <- dba.plotProfile(dba_samples, merge = NULL)
dba.plotProfile(profiles)
dev.off()

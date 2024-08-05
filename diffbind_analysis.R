#!/usr/bin/env Rscript

library(argparser, quietly = TRUE)
library(stringr, quietly = TRUE)

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


# Parse the command line arguments
args <- parse_args(p)

sample_sheet <- read.csv(args$samplesheet)
sample_sheet$Condition <- str_trim(sample_sheet$Condition)
if (!args$contrast %in% sample_sheet$Condition) {
    print("The contrast you have inserted is not present in the 'condition' column of the samplesheet")
    stop()
}

contains_only_numbers <- function(x) !grepl("\\D", x)
summit <- args$summit
if (!str_to_lower(summit) %in% c("false", "median") && !contains_only_numbers(summit)) {
    print(paste0("The --summit parameter given is not a valid option. Available 'false', 'median' or a whole number, given: ", summit))
    stop()
}

library(DiffBind, quietly = TRUE)
library(profileplyr, quietly = TRUE)
library(data.table, quietly = TRUE)

me1_og <- dba(sampleSheet = args$samplesheet)

# ANALYSIS FLOW
# me1 <- dba(sampleSheet="tamoxifen.csv") %>%
#     + dba.blacklist()                         %>%
#     + dba.count()                             %>%
#     + dba.normalize()                         %>%
#     + dba.contrast()                          %>%
#     + dba.analyze()


# Remove one trailing slashes
args$outdir <- gsub("/$", "", args$outdir)
args$outdir <- gsub("\\$", "", args$outdir)
dir.create(args$outdir, showWarnings = FALSE)


print("Saving correlation heatmap generated using the called peaks as correlation_hm_peaks_occupancy.pdf")
pdf(file.path(args$outdir, "correlation_hm_peaks_occupancy.pdf"))
plot(me1_og)
dev.off()


if (str_to_lower(summit) == "false") {
    summit <- FALSE
} else if (str_to_lower(summit) == "median") {
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


me1 <- dba.count(me1_og, summit = summit)

# Save plot as pdf
print("Saving correlation heatmap generated using the affinity (n. of reads in consensous peaks, I think) as correlation_hm_peaks_affinity.pdf")
pdf(file.path(args$outdir, "correlation_hm_peaks_affinity.pdf"))
plot(me1)
dev.off()


print("Showing the number of reads that are in consensous peaks: FRiP (Fraction of Reads in Peaks)")
info <- dba.show(me1)
libsizes <- cbind(LibReads = info$Reads, FRiP = info$FRiP, PeakReads = round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
libsizes
fwrite(data.frame(libsizes, row.names = info$ID), file = file.path(args$outdir, "FRiP_per_sample.tsv"), sep = "\t", row.names = TRUE)


print("Normalizing based on sequencing depth")
me1 <- dba.normalize(me1)


print("Modeling study design based on Condition and performing the differential analysis")
me1 <- dba.contrast(me1, reorderMeta = list(Condition = args$contrast))
me1 <- dba.analyze(me1)


# Save plot as pdf
print("Saving correlation heatmap generated using only significantly differentially bound sites as correlation_hm_peaks_significant.pdf")
pdf(file.path(args$outdir, "correlation_hm_peaks_significant.pdf"))
plot(me1, contrast = 1)
dev.off()


print("Computing report of differentially bound sites...")
me1.DB <- dba.report(me1)
print(me1.DB)


# Write differentially bound sites to a file in a format easy to read for HOMER annotation tool
db <- me1.DB
db$PeakId <- names(db)
db <- data.frame(db)
names(db)[names(db) == "seqnames"] <- "Chr"

# Reorder columns to move PeakId at the start and have the first columns as: PeakId, Chr, start, end, strand
db <- db[, c(1, 2, 3, 12, 9, 5, 6, 7, 8, 10, 11, 4)]
# Rename the first column 'seqnames' to 'Chr'
fwrite(db, file = file.path(args$outdir, "differentially_bound_sites.tsv"), sep = "\t")


num_db_sites <- sprintf(
    "Number of enriched sites in control samples: %s \nNumber of enriched sites in treated samples: %s",
    sum(me1.DB$Fold > 0),
    sum(me1.DB$Fold < 0)
)
print(num_db_sites)
fwrite(list(num_db_sites), file = file.path(args$outdir, "num_enriched_sites.txt"), quote = FALSE)


print("Saving PCA plot of normalized read counts as pca_read_counts.pdf")
pdf(file.path(args$outdir, "pca_read_counts.pdf"))
dba.plotPCA(me1, DBA_TISSUE, label = "DBA_CONDITION")
dev.off()


print("Saving PCA plot of significantly differentially bound sites as pca_significant_sites.pdf")
pdf(file.path(args$outdir, "pca_significant_sites.pdf"))
dba.plotPCA(me1, contrast = 1, label = "DBA_CONDITION")
dev.off()


print("Saving MA plot of treated vs control samples as MA_treated_vs_control.pdf")
pdf(file.path(args$outdir, "MA_treated_vs_control.pdf"))
dba.plotMA(me1)
dev.off()


print("Saving vulcano plot of log FoldChange as vulcano_plot.pdf")
pdf(file.path(args$outdir, "vulcano_plot.pdf"))
dba.plotVolcano(me1)
dev.off()


print("Saving boxplot of concentration of treated vs control samples as boxplot_treated_vs_control.pdf")
pdf(file.path(args$outdir, "boxplot_treated_vs_control.pdf"))
pvals <- dba.plotBox(me1)
dev.off()
fwrite(pvals, file = file.path(args$outdir, "pvals_boxplot.tsv"), sep = "\t", row.names = TRUE)


print("Saving correlation heatmap of each binding site as correlation_hm_binding_sites.pdf")
pdf(file.path(args$outdir, "correlation_hm_binding_sites.pdf"))
hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
dba.plotHeatmap(me1, contrast = 1, correlations = FALSE, scale = "row", colScheme = hmap)
dev.off()


print("Saving profile plot of each binding site as profile_plot_binding_sites.pdf")
pdf(file.path(args$outdir, "profile_plot_binding_sites.pdf"))
profiles <- dba.plotProfile(me1)
dba.plotProfile(profiles)
dev.off()


print("Saving profile plot of each binding site of merged groups (TISSUE and REPLICATES) as profile_plot_group_binding_sites.pdf")
pdf(file.path(args$outdir, "profile_plot_group_binding_sites.pdf"))
profiles <- dba.plotProfile(me1, merge = c(DBA_TISSUE, DBA_REPLICATE))
dba.plotProfile(profiles)
dev.off()


print("Saving profile plot of each binding site of each sample as profile_plot_sample_binding_sites.pdf")
pdf(file.path(args$outdir, "profile_plot_sample_binding_sites.pdf"))
profiles <- dba.plotProfile(me1, merge = NULL)
dba.plotProfile(profiles)
dev.off()

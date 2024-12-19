#!/usr/bin/env python3

import argparse
import os
from pathlib import Path


### Parse cli arguments
argparser = argparse.ArgumentParser(
    prog="homer2igv",
    description="Transform the output of HOMER annotation to a bed file readable by IGV",
)
argparser.add_argument(
    "FILE",
    type=str,
    help="Output from HOMER annotation as .tsv file",
)

argparser.add_argument(
    "-d",
    "--diffbind",
    nargs="?",
    type=str,
    help="Table of differentially significant peaks given to HOMER. This will be used to get the differentially bound strength",
)
filters = argparser.add_argument_group("filters")

filters.add_argument(
    "-f",
    "--fold",
    help="Keep only peaks with abs(log10 fold change) >= FOLD",
    type=float,
)
filters.add_argument(
    "-r",
    "--fdr",
    help="Keep only peaks with -log10(fdr) >= FDR",
    type=float,
)
filters.add_argument(
    "-p",
    "--pvalue",
    help="Keep only peaks with -log10(pvalue) >= PVALUE",
    type=float,
)
filters.add_argument(
    "-t",
    "--tss",
    help="Keep only peaks with abs(distance from tss) <= TSS",
    type=int,
)

filters.add_argument(
    "--type",
    help="Keep only peaks that have type TYPE ['exon', 'intron', 'intergenic', ...]. \
        To add multiple types divide them with commas",
    type=str,
    nargs="*",
    action="append",
)

filters.add_argument(
    "--no-type",
    help="Filter out peaks that have type TYPE ['exon', 'intron', 'intergenic', ...]. \
        To add multiple types divide them with commas",
    type=str,
    nargs="*",
    action="append",
)

argparser.add_argument(
    "-o",
    "--output",
    default="./",
    type=str,
    help="Outputs name. Use a trailing slash to place the outputs in a directory. Default: 'homer2igv_results/'",
)


args = argparser.parse_args()


### Check arguments
input_file = Path(args.FILE)
diffbind_file = Path(args.diffbind) if args.diffbind else None

if not input_file.exists():
    print(f"The input FILE given does not exist: {input_file}")
    exit(1)

if diffbind_file and not diffbind_file.exists():
    print(f"The diffbind file given does not exist: {diffbind_file}")
    exit(1)

if not args.output:
    print("You have to either give an output or let the default one.")
    exit(1)

if args.output[-1] == os.sep:  # Directory
    output_name = Path(args.output, input_file)
else:
    output_name = Path(args.output)

output_name.parent.mkdir(exist_ok=True, parents=True)

tsv_output = output_name.name + ".tsv"
gtf_output = output_name.name + ".gtf"
if Path(tsv_output).exists():
    print(f'Output "{tsv_output}" already exists')
    exit(1)
if Path(gtf_output).exists():
    print(f'Output "{gtf_output}" already exists')
    exit(1)

# Defer the slow import of pandas and numpy to after the arguments are parsed
import pandas as pd
import numpy as np

### Read input files
df = pd.read_csv(input_file, delimiter="\t")
# The first column in the output from HOMER is `PeakID [COMMAND RUN]`, remove the cmd and keep only `peakid`
df = df.rename(columns=lambda x: "peakid" if "peakid" in x.lower() else x)
df = df.rename(columns=lambda x: x.strip().lower().replace(" ", "_"))


# Join HOMER output with DiffBind output on PeakID
if diffbind_file:
    # Read input file and cleanup column names
    db_df = pd.read_csv(diffbind_file, delimiter="\t")
    db_df = db_df.rename(columns=lambda x: x.strip().lower().replace(" ", "_"))
    db_df = db_df.set_index("peakid")
    db_df = db_df.drop(columns=["start", "end", "strand", "chr"])

    df = df.join(db_df, on="peakid", how="inner", rsuffix="_db")

# Drop columns with only NaN
df.dropna(axis=1, how="all", inplace=True)

if "score" not in df.columns:
    df["score"] = 0.0

if "source" not in df.columns:
    df["source"] = "diffbind_homer"

if "frame" not in df.columns:
    df["frame"] = "-"

if "distance_to_tss" not in df.columns:
    df["distance_to_tss"] = "-"

if "feature" not in df.columns:
    # The annotation is usually like "intron (ENST..., intron 1 of 4)"
    # I want to keep only the first word thus "intron", "exon", "CDS", ...
    df["feature"] = df["annotation"].apply(lambda x: str(x).strip().split(" ")[0])

if "name" not in df.columns:
    # Name to be displayed in IGV
    df["name"] = df["annotation"] + " - " + df["gene_name"]


# IGV automatically reads "gene_name" and displays it under the feature, but we want "name" to be displayed under,
# thus we need to rename "gene_name" adding an underscore to it. Then we remove the originl column
df["gene_name_"] = df["gene_name"]
df = df.drop(columns=["gene_name"])

### Filter the dataframe
if args.fold:
    df = df[df["fold"].abs() >= args.fold]

if args.fdr:
    df = df[-np.log10(df["fdr"]) >= args.fdr]

if args.pvalue:
    df = df[-np.log10(df["p.value"]) >= args.pvalue]

if args.tss:
    df = df[df["distance_to_tss"].abs() <= args.tss]

if args.type:
    df = df[df["feature"] in args.type.split(",")]

if args.no_type:
    df = df[df["feature"] not in args.no_type.split(",")]

print(f"TSV file saved: {tsv_output}")
df.to_csv(
    tsv_output,
    sep="\t",
    header=True,
    index=False,
    encoding="utf-8",
)


### Prepare database to be written
# get name of columns to be added in the `attribute` column: (all - list_below)
attribute_cols: pd.Index = df.columns.drop(
    [
        "chr",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
    ]
)
# Now write each column in the attribute column as "col_name1: value1;col_name2: value2 ..."
df["attribute"] = [
    "; ".join(f"{attribute_cols[i]} {value}" for i, value in enumerate(row))
    for row in df[attribute_cols].itertuples(index=False)
]

df = df[
    [
        "chr",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]
]
print(f"GTF file saved: {gtf_output}")
df.to_csv(gtf_output, sep="\t", header=False, index=False, encoding="utf-8")

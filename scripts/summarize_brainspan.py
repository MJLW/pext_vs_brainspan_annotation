import argparse
import pdb

import polars as pl
import gtfparse



def parse_args():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        "-b", "--brainspan",
        type=str,
        required=True,
        help="Path to the brainspan file."
    )
    parser.add_argument(
        "-g", "--gtf",
        type=str,
        required=True,
        help="Path to the Gencode v10 GTF file."
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        help="Path to the output file (TSV)."
    )

    return parser.parse_args()


def main():
    args = parse_args()

    df_gene_chr_map = gtfparse.read_gtf(args.gtf) \
        .filter(pl.col("feature") == "gene") \
        .select(["seqname", "gene_id"]) \
        .unique() \
        .with_columns(pl.col("gene_id").str.split(".").list.get(0).cast(pl.Categorical).alias("gene_id")) \
        .rename({"seqname": "chrom"})

    df_brainspan = pl.read_csv(args.brainspan, separator="\t", schema_overrides={"ensembl_gene_id": pl.Categorical}) \
        .drop(["gene_symbol", "age", "rpkm"])

    # Calculate max ratio over time span
    # Median over donor and structure
    # Write to TSV

    pdb.set_trace()


if __name__ == "__main__":
    main()

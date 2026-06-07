import argparse
import polars as pl



def parse_args():
    parser = argparse.ArgumentParser(description="Bin and sum PEXT/BrainSpan scores")

    parser.add_argument(
        "-i", "--input",
        type=str,
        required=True,
        help="Path to the input file."
    )
    parser.add_argument(
        "-t", "--tag",
        type=str,
        choices=["PEXT", "BrainSpan"],
        required=True,
        help="Tag to apply: 'PEXT' or 'BrainSpan'."
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        help="Path to the output file."
    )

    return parser.parse_args()


def main():
    args = parse_args()

    df = pl.read_csv(args.input, separator="\t", schema_overrides={"MR": pl.String, "DISTAL": pl.Categorical, args.tag: pl.String}) \
        .filter(pl.col("CONSEQUENCE").str.contains("NMD_transcript").not_()) \
        .with_columns(pl.col("CONSEQUENCE").str.split("&").list.get(0).alias("CONSEQUENCE")) \
        .drop(["TRANSCRIPT"]) \
        .unique() \
        .with_columns(
            pl.col("MR").cast(pl.Float16, strict=False).alias("MR"),
            pl.col(args.tag).cast(pl.Float16, strict=False).alias(args.tag),
        ) \
        .filter(pl.col("MR").is_not_null() & pl.col(args.tag).is_not_null()) \
        .with_columns(
            pl.col(args.tag).cut(breaks=[0.2, 0.4, 0.6, 0.8], labels=['0.0-0.2', '0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1.0']).alias("BIN"),
        ) \
        .group_by(["GENE", "BIN", "CONSEQUENCE", "DISTAL"]) \
        .agg(pl.sum("MR")) \
        .sort(["GENE", "BIN", "DISTAL"])

    df.write_csv(args.output, separator="\t")


if __name__ == "__main__":
    main()


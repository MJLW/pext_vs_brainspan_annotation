import pdb

import polars as pl



def merge_intervals(df: pl.DataFrame) -> pl.DataFrame:
    group_cols = ['chrom', 'transcript_id', 'gene_id']

    return (
        df
        .sort(group_cols + ['start'])
        .with_columns(
            pl.col('end')
            .shift(1)
            .cum_max()
            .over(group_cols)
            .alias('_prev_max_end')
        )
        .with_columns(
            (
                pl.col('_prev_max_end').is_null() |         # first row of each group
                (pl.col('start') > pl.col('_prev_max_end')) # gap detected
            )
            .cast(pl.Int32)
            .cum_sum()
            .over(group_cols)
            .alias('_group')
        )
        .group_by(group_cols + ['_group'])
        .agg(
            pl.col('start').min(),
            pl.col('end').max(),
            pl.col('distal').any(),
        )
        .drop('_group')
        .sort(group_cols + ['start'])
    )


def main():
    df = pl.read_csv("in/nmd_region_annotations.tsv", separator="\t") \
        .select(["chrom", "start", "end", "transcript_id", "gene_id", "distal"])

    pdb.set_trace()

    pass


if __name__ == "__main__":
    main()

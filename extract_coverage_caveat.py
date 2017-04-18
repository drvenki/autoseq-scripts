#!/usr/bin/env python
import click, logging, sys
import pandas as pd


def extract_qc_call(coverage_histogram_file,
                    high_thresh_fraction,
                    high_thresh_fold_cov,
                    low_thresh_fraction,
                    low_thresh_fold_cov):
    '''
    Generate a QC call on coverage for the specified coverage distribution,
    by considering the specified thresholds.

    If the sample has the designated fraction over the upper threshold, then "OK".
    Otherwise, if the sample has the designated fraction over the lower threshold,
    then "WARN". Otherwise, "FAIL".

    :param coverage_histogram_file: Output file from target_coverage_histogram.py
    :param high_thresh_fraction: Upper threshold fraction requirement
    :param high_thresh_fold_cov: Upper threshold fold change requirement
    :param low_thresh_fraction: Lower threshold fraction requirement
    :param low_thresh_fold_cov: Lower threshold fold change requirement
    :return:
    '''

    try:
        coverage_table = pd.read_csv(coverage_histogram_file, comment="#", sep="\t", header=None)
    except pd.io.common.EmptyDataError, e:
        # Will assume that no data in the coverage histogram indicates no coverage:
        return "FAIL"

    cum_fraction = coverage_table.cumsum(0)[3]
    coverage_table["cum_fraction"] = cum_fraction

    # Retrieve coverage at lowest observed coverage that has coverage
    # > high_thresh_fraction of the data:
    max_percent_data_under_high_thresh = 1 - high_thresh_fraction
    max_percent_data_under_low_thresh = 1 - low_thresh_fraction

    obs_cov_high_thresh = coverage_table.loc[coverage_table["cum_fraction"] >= max_percent_data_under_high_thresh, 0].iloc[0]
    if obs_cov_high_thresh > high_thresh_fold_cov:
        return "OK"

    obs_cov_high_thresh = coverage_table.loc[coverage_table["cum_fraction"] >= max_percent_data_under_low_thresh, 0].iloc[0]
    if obs_cov_high_thresh > low_thresh_fold_cov:
        return "WARN"

    else:
        return "FAIL"


def write_qc_json(output_file, qc_call):
    output = '''{
"CALL": "%s"
}''' % (qc_call)
    print >> output_file, output


@click.command()
@click.option('--high-thresh-fraction', default=0.95, help='Upper threshold, fraction requirement')
@click.option('--high-thresh-fold-cov', default=100, help='Upper threshold, fold coverage requirement')
@click.option('--low-thresh-fraction', default=0.95, help='Low threshold, fraction requirement')
@click.option('--low-thresh-fold-cov', default=50, help='Low threshold, fold coverage requirement')
@click.argument('coverage-histogram')
def main(coverage_histogram, high_thresh_fraction, high_thresh_fold_cov,
         low_thresh_fraction, low_thresh_fold_cov):
    """
Generate a call on contamination from contest output, and output it as a JSON file.

Example output contents: {"CALL": "OK"}
"""

    numeric_level = getattr(logging, "INFO", None)
    logging.basicConfig(level=numeric_level)

    coverage_histogram_file = open(coverage_histogram)
    qc_call = extract_qc_call(coverage_histogram_file,
                              high_thresh_fraction,
                              high_thresh_fold_cov,
                              low_thresh_fraction,
                              low_thresh_fold_cov)

    write_qc_json(sys.stdout, qc_call)


if __name__ == "__main__":
    sys.exit(main())

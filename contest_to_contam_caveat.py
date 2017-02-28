#!/usr/bin/env python
import click, logging, pdb, sys


def extract_qc_call(contest_output_file, max_contam):
    # Ignore header:
    contest_output_file.readline()

    # Extract the fourth field:
    contam_estimate = float(contest_output_file.readline().strip().split()[3])

    # Return QC call based on this estimate and the max acceptable contamination:
    if contam_estimate > max_contam:
        return "FAIL"
    else:
        return "OK"


def write_qc_json(output_file, qc_call):
    output = '''{
"CALL": "%s"
}''' % (qc_call)
    print >> output_file, output


@click.command()
@click.option('--loglevel', default='INFO', help='level of logging')
@click.option('--max-contam', default=1.0, help='level of logging')
@click.argument('contest-results')
def main(contest_results, loglevel, max_contam):
    """
Generate a call on contamination from contest output, and output it as a JSON file.

Example output contents: {"CALL": "OK"}
"""

    numeric_level = getattr(logging, "INFO", None)
    logging.basicConfig(level=numeric_level)

    logging.info("Parsing contest output...")
    contest_output_file = open(contest_results)
    qc_call = extract_qc_call(contest_output_file, max_contam)

    write_qc_json(sys.stdout, qc_call)


if __name__ == "__main__":
    sys.exit(main())

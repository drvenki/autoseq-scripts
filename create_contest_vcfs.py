#!/usr/bin/env python
import click, logging, sys


@click.command()
@click.option('--loglevel', default='INFO', help='level of logging')
@click.option('--tmpdir', default='/tmp', help='Temp directory for intermediate files')
@click.argument('genotype-target-bedfile')
@click.argument('eval-target-bedfile')
@click.argument('population-vcf')
def main(genotype_target_bedfile, eval_target_bedfile, population_vcf, loglevel, tmpdir):
    """
Create the contest VCF input file from the designated "genotype" and "eval" target bed region
files and population allele frequency VCF file.
"""

    numeric_level = getattr(logging, loglevel, None)
    logging.basicConfig(level=numeric_level)

    # Get the intersection of the two bed files:
    # XXX CONTINUE HERE; IMPLEMENT THESE STEPS (SEPARATE FUNCTIONS + THEIR UNIT TESTS)

    # Get the intersection of the resulting BED file with the population VCF file:
    # XXX

    # Reformat the resulting VCF file to be consistent with contest:
    # XXX


if __name__ == "__main__":
    sys.exit(main())

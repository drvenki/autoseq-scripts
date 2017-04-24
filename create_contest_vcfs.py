#!/usr/bin/env python
import click, logging, StringIO, sys
import pybedtools
import vcf


@click.command()
@click.option('--loglevel', default='INFO', help='level of logging')
@click.option('--tmpdir', default='/tmp', help='Temp directory for intermediate files')
@click.option('--max-maf', default=0.9, help='Maximum VCF record MAF')
@click.option('--min-maf', default=0.1, help='Minimum VCF record MAF')
@click.option('--output-filename', default='intersect_snps.vcf', help='Output filename')
@click.argument('genotype-target-bedfile')
@click.argument('eval-target-bedfile')
@click.argument('population-vcffile')
def main(genotype_target_bedfile, eval_target_bedfile, population_vcffile, output_filename, min_maf, max_maf, tmpdir, loglevel):
    """
Create the contest VCF input file from the designated "genotype" and "eval" target bed region
files and population allele frequency VCF file.
"""

    numeric_level = getattr(logging, loglevel, None)
    logging.basicConfig(level=numeric_level)

    # Get the intersection of the two bed files:
    logging.info("Getting genotype and eval bed intersections...")
    genotype_target_bed = pybedtools.BedTool(genotype_target_bedfile)
    eval_target_bed = pybedtools.BedTool(eval_target_bedfile)
    bed_intersection = genotype_target_bed.intersect(eval_target_bed)
    logging.info("Done.")

    # Get the intersection of the resulting BED file with the population VCF file:
    logging.info("Getting population vcf intersection...")
    pop_vcf = pybedtools.BedTool(population_vcffile)
    vcf_bed_intersect = pop_vcf.intersect(bed_intersection)
    logging.info("Done.")

    # Get a list of VCF strings from the resulting intersection information:
    # NOTE: Could consume a fair bit of memory if the above bed-vcf filter is not aggressive:
    logging.info("Preparing output stream...")
    vcf_bed_intersect_df = vcf_bed_intersect.to_dataframe()
    vcf_bed_intersect_list = list(\
        vcf_bed_intersect_df.apply(lambda elems: "\t".join(map(lambda elem: str(elem), elems)), axis=1))

    # Convert to a file-like object to facilitate filtering and output with pyvcf:
    vcf_bed_intersect_bigstr = "\n".join(vcf_bed_intersect_list)
    vcf_bed_intersect_filelike = StringIO.StringIO(vcf_bed_intersect_bigstr)
    logging.info("Done.")

    # Open outfile: FIXME: Improve this:
    logging.info("Filtering and writing output records...")
    with open(output_filename, 'w') as output_file, open(population_vcffile) as template_vcf:
        vcf_header_template = vcf.Reader(template_vcf)

        # Parse using pyvcf and apply filter, by explicit iteration over entries:
        vcf_bed_intersect_reader = vcf.Reader(vcf_bed_intersect_filelike)
        vcf_bed_intersect_filtered_writer = vcf.Writer(output_file, vcf_header_template)
        for record in vcf_bed_intersect_reader:
            if record.INFO['AF'][0] > min_maf and record.INFO['AF'][0] < max_maf:
                vcf_bed_intersect_filtered_writer.write_record(record)
    logging.info("Done.")


if __name__ == "__main__":
    sys.exit(main())

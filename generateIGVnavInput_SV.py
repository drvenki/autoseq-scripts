#!/usr/bin/env python

import subprocess
import argparse
import os
import glob


def parse_svaba(input_vcf, SDID, output, vcftype):
    """
    vawk '{print $1, $2, $2+1, "P-00356971_svaba", "BND", $5, S$*$AD, S$*$DP}'
    """
    header = "echo \"CHROM\tSTART\tEND\tSDID\tSVTYPE\tALT\tSUPPORT_normal\tSUPPORT_tumor \
             \tDPnormal\tDPtumor\tGENES\"" + " > " + output + "_" + vcftype + "_svaba.mut"
    svaba_cmd = "vawk '{print $1, $2, $2+1, \"" + SDID + "\",\"I$SVTYPE\", $5, S$*$AD,S$*$DP}'" + \
                " " + input_vcf + " >> " + output + "_" + vcftype + "_svaba.mut"
    
    # cmd = "awk 'NR>1 {OFS=\"\\t\";print $1, $2, $3, $5,\"svaba\", $4, \"" + vcftype + "\", $6, " + sup_reads + "}' " + output + \
    #       "_" + vcftype + "_svaba.mut >> " + output_dir + "/annotate_combined_sv.txt"
    subprocess.call(" && ".join([header, svaba_cmd]), shell=True)


def parse_svict(input_vcf, SDID, output, vcftype):
    """
    CHROM START   END SDID    SVTYPE  ALT SUPPORTING_READS    NOTES   GENES
    vawk '{if (I$SUPPORT>8 && I$SVTYPE ~ "INS" ) print $1, $2, I$END, "$SDID_svict",
                                                        I$SVTYPE, "<INS>", I$SUPPORT
    else if (I$SUPPORT>8 && I$SVTYPE ~ "INV" ) print $1, $2, I$END, "P-00356971_svict",
                                                         I$SVTYPE, "<INV>", I$SUPPORT
    else if (I$SUPPORT>8 && I$SVTYPE ~ "DEL" ) print $1, $2, I$END, "P-00356971_svict",
                                                         I$SVTYPE, "<DEL>", I$SUPPORT
    else if (I$SUPPORT>8 && I$SVTYPE ~ "BND" ) print $1, $2, $2+1, "$SDID_svict",
                                                         I$SVTYPE, $5, I$SUPPORT}' $input.vcf
    """
    header8 = "echo \"CHROM\tSTART\tEND\tSDID\tSVTYPE\tALT\tSUPPORT_READS\"" + \
              " > " + output + "_SR8.mut"

    header12 = "echo \"CHROM\tSTART\tEND\tSDID\tSVTYPE\tALT\tSUPPORT_READS\"" + \
               " > " + output + "_SR12.mut"

    sup_8 = "vawk '{if (I$SUPPORT>8 && I$SVTYPE ~ \"INS\" ) print $1, $2, I$END, \"" + \
            SDID + "\", I$SVTYPE, \"<INS>\", I$SUPPORT ;" + \
            " else if (I$SUPPORT>8 && I$SVTYPE ~ \"INV\" ) print $1, $2, I$END, \"" + \
            SDID + "\", I$SVTYPE, \"<INV>\", I$SUPPORT ; " + \
            " else if (I$SUPPORT>8 && I$SVTYPE ~ \"DEL\" ) print $1, $2, I$END, \"" + \
            SDID + "\", I$SVTYPE, \"<DEL>\", I$SUPPORT ; " + \
            " else if (I$SUPPORT>8 && I$SVTYPE ~ \"BND\" ) print $1, $2, $2+1, \"" + \
            SDID + "\", I$SVTYPE, $5, I$SUPPORT}' " + input_vcf + \
            " >> " + output + "_SR8.mut"

    sup_12 = "vawk '{if (I$SUPPORT>12 && I$SVTYPE ~ \"INS\" ) print $1, $2, I$END, \"" + \
             SDID + "\", I$SVTYPE, \"<INS>\", I$SUPPORT ; " + \
             " else if (I$SUPPORT>12 && I$SVTYPE ~ \"INV\" ) print $1, $2, I$END, \"" + \
             SDID + "\", I$SVTYPE, \"<INV>\", I$SUPPORT ; " + \
             " else if (I$SUPPORT>12 && I$SVTYPE ~ \"DEL\" ) print $1, $2, I$END, \"" + \
             SDID + "\", I$SVTYPE, \"<DEL>\", I$SUPPORT ;" + \
             " else if (I$SUPPORT>12 && I$SVTYPE ~ \"BND\" ) print $1, $2, $2+1, \"" + \
             SDID + "\", I$SVTYPE, $5, I$SUPPORT}' " + input_vcf + \
             " >> " + output + "_SR12.mut"

    # cmd = "awk ' NR>1 {OFS=\"\\t\"; print $1, $2, $3, $5,\"svict\", $4,\"" + vcftype + "\", $6, $7}' " + output + "_SR8.mut "\
    #       + " >> " + output_dir + "/annotate_combined_sv.txt"

    subprocess.call(" && ".join([header8, header12, sup_8, sup_12]), shell=True)


def parse_lumpy(input_vcf, SDID, output, vcftype):
    """
    ##CHROM START   END SDID    SVTYPE  ALT SVLENGTH    SUPPORT_READS    NOTES   GENES
    vawk '{if ((I$SVLEN>1000 || I$SVLEN<-1000) && $1 != "hs37d5" && $1 !~ "GL" && $5 !~ "hs37d5" &&
    I$SU>50 && I$SVTYPE ~ "BND") print $1, $2, $2+1, "P-00356971_lumpy", I$SVTYPE, $5, "NA", I$SU
    else if ((I$SVLEN>1000 || I$SVLEN<-1000) && $1 != "hs37d5" && $1 !~ "GL" && $5 !~ "hs37d5"
    && I$SU>50 && I$SVTYPE !~ "BND") print $1, $2, I$END, "P-00356971_lumpy", I$SVTYPE, $5, I$SVLEN, I$SU}
    """
    header1k_sup_50 = "echo \"CHROM\tSTART\tEND\tSDID\tSVTYPE\tALT\tSUPPORT_READS\"" + \
                      " > " + output + "_lumpy_len1k_SU50.mut"
    header500_sup_24 = "echo \"CHROM\tSTART\tEND\tSDID\tSVTYPE\tALT\tSUPPORT_READS\"" + \
                       " > " + output + "_lumpy_len500_SU24.mut"

    len1k_sup_50 = "vawk '{if ((I$SVLEN>1000 || I$SVLEN<-1000) && $1 != \"hs37d5\" && $1 !~ \"GL\" && $5 !~ \"hs37d5\" && I$SU>50 && I$SVTYPE ~ \"BND\") print $1, $2, $2+1, \""+ SDID +"\", I$SVTYPE, $5, \"NA\", I$SU ;"  + \
                " else if ((I$SVLEN>1000 || I$SVLEN<-1000) && $1 != \"hs37d5\" && $1 !~ \"GL\" && $5 !~ \"hs37d5\" && I$SU>50 && I$SVTYPE !~ \"BND\") print $1, $2, I$END, \""+ SDID +"\", I$SVTYPE, $5, I$SU}' " + input_vcf + \
                " >> " + output + "_lumpy_len1k_SU50.mut"

    len500_sup_24 = "vawk '{if ((I$SVLEN>500 || I$SVLEN<-500) && $1 != \"hs37d5\" && $1 !~ \"GL\" && $5 !~ \"hs37d5\" && I$SU>24 && I$SVTYPE ~ \"BND\") print $1, $2, $2+1, \"" + SDID + "\", I$SVTYPE, $5, I$SU ;" + \
                " else if ((I$SVLEN>500 || I$SVLEN<-500) && $1 != \"hs37d5\" && $1 !~ \"GL\" && $5 !~ \"hs37d5\" && I$SU>24 && I$SVTYPE !~ \"BND\") print $1, $2, I$END, \"" + SDID + "\", I$SVTYPE, $5, I$SU }' " + input_vcf + \
                " >> " + output + "_lumpy_len500_SU24.mut"
  
    # cmd = "awk 'NR>1 {OFS=\"\\t\";print $1, $2, $3, $5,\"lumpy\", $4, \"" + vcftype + "\", $6, $7}' " + output + "_lumpy_len500_SU24.mut" +  " >> " + output_dir + "/annotate_combined_sv.txt"

    subprocess.call(" && ".join([header1k_sup_50, header500_sup_24, len1k_sup_50, len500_sup_24]), shell=True)


def combine_mut(input_dir, output_dir):

    files = glob.glob(input_dir + "/*.mut")
    cmd = []

    if not os.path.exists(output_dir + "/annotate_combined_sv.txt"):
        header2 = "echo \"CHROM\tSTART\tEND\tSVTYPE\tTOOL\tSDID\tSAMPLE\tALT\tSUPPORT_READS\"" + " >> " + output_dir + "/annotate_combined_sv.txt"
        subprocess.call(header2, shell=True)

    for file in files:
        vcftype = ''
        sup_reads = ''
        if 'svict_SR8' in file:
            vcftype = 'cfdna' if 'CFDNA' in file else 'germline'
            cmd.append("awk ' NR>1 {OFS=\"\\t\"; print $1, $2, $3, $5,\"svict\", $4,\"" + vcftype + "\", $6, $7}' " \
                        + file + " >> " + output_dir + "/annotate_combined_sv.txt")
        elif 'lumpy_len500_SU24' in file:
            cmd.append("awk 'NR>1 {OFS=\"\\t\";print $1, $2, $3, $5,\"lumpy\", $4, \"somatic\", $6, $7}' " \
                        + file +  " >> " + output_dir + "/annotate_combined_sv.txt")
        elif 'svaba.mut' in file:
            vcftype = 'somatic' if 'somatic' in file else 'germline'
            sup_reads = '$8' if vcftype == 'SOMATIC' else '$7'
            cmd.append("awk 'NR>1 {OFS=\"\\t\";print $1, $2, $3, $5,\"svaba\", $4, \"" + vcftype + "\", $6, " + sup_reads + "}' " +\
                        file + " >> " + output_dir + "/annotate_combined_sv.txt")

    subprocess.call(" && ".join(cmd), shell=True)

    return output_dir + "/annotate_combined_sv.txt"


def load_bed(bed_file):
    genes = {}
    with open(bed_file, 'r') as genes_fh:
        genes_db = genes_fh.readlines()
        for each_entry in genes_db:
            data = each_entry.strip().split('\t')
            chrom = data[0]
            start = data[1]
            end = data[2]
            gene = data[3]

            if chrom in genes:
                genes[chrom].update({(start, end): gene})
            else:
                genes[chrom] = {(start, end): gene}
    return genes


def gene_annotation(chrom, start, end, genes, design):
    gene = ''
    in_gene = ''

    try:
        # annotate gene name
        for ranges, gene_name in genes[chrom].items():
            if int(ranges[0]) - 20 <= int(start) <= int(ranges[1]) + 20 or int(ranges[0]) - 20 <= int(end) <= int(ranges[1]) + 20:
                gene = gene_name
                break
        if not gene:
            gene = 'None'

        # annotate pancancer info
        for pan_ranges, gene_name in design[chrom].items():
            if int(pan_ranges[0]) <= int(start) <= int(pan_ranges[1]) or int(pan_ranges[0]) <= int(end) <= int(pan_ranges[1]):
                in_gene = 'YES'
                break
        if not in_gene:
            in_gene = 'NO'

        return (gene, in_gene)
    except KeyError:
        return ('NA', 'NA')
        print "Warning! chromosome {chrom} is not valid".format(chrom=chrom)


def annotate_combined_sv(combined_file, genes, pancancer, output):
    output_file = open(output, 'w')
    print (combined_file)
    with open(combined_file, 'r') as file:
        header = file.readline()
        output_file.write("\t".join(['CHROM_A', 'START_A', 'END_A', 'CHROM_B', 'START_B', 'END_B',
            'SVTYPE', 'SV_LENGTH', 'SUPPORT_READS', 'TOOL', 'SDID', 'SAMPLE', 'GENE_A', 'IN_DESIGN_A', 'GENE_B', 'IN_DESIGN_B', "GENE_A-GENE_B-sorted"]) + '\n')
        for line in file.readlines():
            data = line.strip().split('\t')
            chrom_a = data[0]
            start_a = data[1]
            end_a = data[2]
            svtype = data[3]
            tool = data[4]
            sdid = data[5].split('_')[0]
            sample = data[6]
            sup_reads = data[8] if len(data) == 9 else '.'
            #print (line)
            if ':' in data[7]:
                chrom_b = filter(str.isdigit, data[7].split(':')[0])
                start_b = int(filter(str.isdigit, data[7].split(':')[1]))
                end_b = start_b + 1
            else:
                chrom_b = 'NA'
                start_b = 'NA'
                end_b = 'NA'

            gene_a, in_gene_a = gene_annotation(chrom_a, start_a, end_a, genes, pancancer)

            if chrom_b != 'NA':
                gene_b, in_gene_b = gene_annotation(chrom_b, start_b, end_b, genes, pancancer)
            else:
                gene_b, in_gene_b = 'NA', 'NA'

            gene_a_b = [gene_a, gene_b]
            gene_a_b.sort()
            gene_a_b_sorted = ",".join(gene_a_b)

            # write output
            #print (data)
	    output_file.write("\t".join(map(str, [chrom_a, start_a, end_a, chrom_b, start_b,
                             end_b, svtype, int(end_a)-int(start_a), sup_reads, tool, sdid, sample, gene_a, in_gene_a, gene_b, in_gene_b, gene_a_b_sorted])) + '\n')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=
        'A MUT file (.mut) is a tab-delimited text file that lists mutations. \
        The first row contains column headings and each subsequent row identifies a mutation. \
        IGV ignores the column headings.It reads the first five columns as shown below and \
        ignores all subsequent columns:  \
        1. chromosome \
        2. start location (location of the first base pair in the mutated region) \
        3. end location (location of the last base pair in the mutated region) \
        4. sample or patient ID \
        5. mutation type (for example, Synonymous, Missense, Nonsense, Indel, etc.)')
    parser.add_argument('--input', required=True, help="Input VCF or tab-delimited file")
    parser.add_argument('--annotBed', help="UCSC hg19 genes bed file with chrom, start, \
                        end and genesymbol")
    parser.add_argument('--targetBed', help="Bed file to annotate targets eg: pancancer")
    parser.add_argument('--sdid', help="SDID from analysis")
    parser.add_argument('--vcftype', help="somatic (or) germline vcf (only for svaba)")
    parser.add_argument('--tool', help="Tool name - Variant callers")
    parser.add_argument('--output', required=True,
                        help="output tab delimited file for IGVNav, format=output.mut")
    args = parser.parse_args()


    vcftype = args.vcftype
    input_file = args.input
    annotBed = args.annotBed
    sv_caller = args.tool
    sdid = args.sdid
    output = args.output
    target_bed = args.targetBed


    output_dir = os.path.dirname(output)


    if sv_caller == 'svict':
        parse_svict(input_file, sdid, output, vcftype)
    elif sv_caller == 'lumpy':
        parse_lumpy(input_file, sdid, output, vcftype)
    elif sv_caller == 'svaba':
        parse_svaba(input_file, sdid, output, vcftype)

    if annotBed:
        combined_input = combine_mut(input_file, output_dir)
        genes = load_bed(annotBed)
        targets = load_bed(target_bed)
        annotate_combined_sv(combined_input, genes, targets, output)


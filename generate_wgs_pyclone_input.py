from argparse import ArgumentParser
import os
import pandas as pd
import sys
from collections import defaultdict
import numpy as np

"""


Usage:
    generate_wgs_pyclone_input.py 
        SNVs_OxoGFFPE/GCT001-TP-NT-SM-DPBZL-SM-DPBZK.ffpeBias.maf.annotated 
        Indels/GCT001-TP-NT-SM-DPBZL-SM-DPBZK.indel.wgs.maf.annotated 
        CNAs/GCT001-TP-NT-SM-DPBZL-SM-DPBZK.titan.optimalClust_ploidy2.segs.txt 
        output_folder

Can be used to generate the tab-delimited file used as input to the build_mutations_file script for PyClone.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
From TITAN instructions:

TSV Input File
The build_mutations_file takes a tab delimited file with a header as input and produces a YAML formatted file which can be used for running a PyClone analysis. Example files are contained in the examples/mixing/tsv folder which ships with the PyClone software.

The required fields in this file are:

mutation_id - A unique ID to identify the mutation. Good names are thing such a the genomic co-ordinates of the mutation i.e. chr22:12345. Gene names are not good IDs because one gene may have multiple mutations, in which case the ID is not unique and PyClone will fail to run or worse give unexpected results. If you want to include the gene name I suggest adding the genomic coordinates i.e. TP53_chr17:753342.

ref_counts - The number of reads covering the mutation which contain the reference (genome) allele.

var_counts - The number of reads covering the mutation which contain the variant allele.

normal_cn - The copy number of the cells in the normal population. For autosomal chromosomes this will be 2 and for sex chromosomes it could be either 1 or 2. For species besides human other values are possible.

minor_cn - The minor copy number of the cancer cells. Usually this value will be predicted from WGSS or array data.

major_cn - The major copy number of the cancer cells. Usually this value will be predicted from WGSS or array data.

If you do not major and minor copy number information you should set the minor copy number to 0, and the major copy number to the predicted total copy number. If you do this make sure to use the total_copy_number for the --var_prior flag of the build_mutations_file command. DO NOT use the parental copy number information method as it assumes you have knowledge of the minor and major copy number.

Any additional columns in the tsv file will be ignored so feel free to add additional annotation fields.

"""

DEFAULT_MAJOR_COPY_NUMBER = 1
DEFAULT_MINOR_COPY_NUMBER = 1
DEFAULT_CELL_COPY_NUMBER = 2


class TitanCols:
    CHR = "Chromosome"
    START = "Start_Position(bp)"
    END = "End_Position(bp)"
    MAJOR_CN = "MajorCN"
    MINOR_CN = "MinorCN"


class MafCols:
    START = "Start_position"
    END = "End_position"
    CHR = "Chromosome"
    CLASSIFICATION = "Variant_Classification"
    TYPE = "Variant_Type"
    GENE = "Hugo_Symbol"
    REF = "Reference_Allele"
    ALT = "Tumor_Seq_Allele2"
    REF_COUNT = "t_ref_count"
    ALT_COUNT = "t_alt_count"

    @classmethod
    def all(cls):
        return [cls.GENE, cls.CHR, cls.START, cls.END, cls.CLASSIFICATION, cls.TYPE, cls.GENE, cls.REF, cls.ALT,
                cls.REF_COUNT, cls.ALT_COUNT]


def mutation_id(row):
    """Generate a unique mutation ID"""
    return 'chr{}:{}({})'.format(row[MafCols.CHR], row[MafCols.START], row[MafCols.GENE].values[0])


def center_position(row):
    """Find the point halfway betwee the start and end of the mutation -- relevant for SNVs"""
    start_int = int(row[MafCols.START])
    end_int = int(row[MafCols.END])
    length = end_int - start_int
    return end_int - (length / 2)


def combine_snv_and_indel_mafs(snv_maf_df, indel_maf_df):
    """Combine data from snv and indel maf dataframes into one dataframe"""
    subset_snvs = snv_maf_df[MafCols.all()]
    subset_indels = indel_maf_df[MafCols.all()]
    combined_snvs_and_indels_df = pd.concat([subset_snvs, subset_indels])

    sys.stdout.write("Adding mutation ID...\n")
    combined_snvs_and_indels_df['mutation_id'] = combined_snvs_and_indels_df.apply(mutation_id, axis=1)
    sys.stdout.write("Added mutation ID...\n")

    sys.stdout.write("Adding section center...\n")
    combined_snvs_and_indels_df['Center_Position'] = combined_snvs_and_indels_df.apply(center_position, axis=1)
    sys.stdout.write("Added section center...\n")
    return combined_snvs_and_indels_df


def total_copy_number_based_on_sex_and_chrom(sex, chrom):
    """Based on whether looking at a sex chromosome or not, determine the normal cell copy number"""
    # Figure out normal copy number based on sex if this a sex chromosome
    if chrom in ['X', 'Y']:
        if sex == 'female':
            assert(chrom == 'X')
            return 2
        elif sex == 'male':
            return 1
    else:
        return DEFAULT_CELL_COPY_NUMBER


def build_titan_search_tree(titan_df, sex):
    """Given a dataframe of Titan data, build a dictionary where the chromosomes are the keys and the values
    are the segments"""
    chr_map = defaultdict(lambda: defaultdict(dict))
    for index, row in titan_df.iterrows():
        chrom = str(row[TitanCols.CHR])
        chr_map[chrom][row[TitanCols.START]] = {"end": row[TitanCols.END],
                                                             "major_cn": row[TitanCols.MAJOR_CN],
                                                             "minor_cn": row[TitanCols.MINOR_CN],
                                                             "normal_cn": total_copy_number_based_on_sex_and_chrom(sex, chrom)
                                                             }
    return chr_map


def add_cn_info_to_indel_snv_maf(indel_snv_maf, titan_search_tree, sex):
    """Add the titan allelic copy number information to the combined snv indel maf"""
    titan_chrs = titan_search_tree.keys()

    def major_and_minor_cns(combined_maf_row):
        # For all rows where the chromosome isn't even in the titan search tree, set the minor and major alleles to the
        # default, non-aneuploidy values
        chrom = combined_maf_row[MafCols.CHR]
        position = combined_maf_row["Center_Position"]

        if chrom not in titan_chrs:
            return DEFAULT_MAJOR_COPY_NUMBER, DEFAULT_MINOR_COPY_NUMBER, total_copy_number_based_on_sex_and_chrom(sex, chrom)
        else:
            starts = list(titan_search_tree[chrom].keys())
            start_index = np.searchsorted(starts, position)
            relevant_start = starts[start_index - 1]
            seg = titan_search_tree[chrom][relevant_start]
            if seg.get("end") >= position:
                return seg.get("major_cn"), seg.get("minor_cn"), seg.get("normal_cn")
            else:
                return DEFAULT_MAJOR_COPY_NUMBER, DEFAULT_MINOR_COPY_NUMBER, total_copy_number_based_on_sex_and_chrom(sex, chrom)
    cn_info = indel_snv_maf.apply(major_and_minor_cns, axis=1)

    # Get a series with a major_cn and minor_cn column of out this Series of tuples
    cn_info_cols = cn_info.apply(pd.Series, index=["major_cn", "minor_cn", "normal_cn"])

    cn_info_added_to_snv_and_maf_df = pd.concat([indel_snv_maf, cn_info_cols], axis=1)
    return cn_info_added_to_snv_and_maf_df


def get_sex(snv_maf):
    """Based on whether both a Y chromosome is present, make a sex determination"""
    chromosomes = snv_maf['Chromosome'].unique()
    y_in_sample = 'Y' in chromosomes
    if y_in_sample:
        return 'male'
    else:
        return 'female'


def main():
    parser = ArgumentParser(description='Generate .tsv input for PyClone from SNV maf,'
                                        ' GATK indel file, and TITAN output')
    parser.add_argument('snv_maf', metavar='snv_maf', type=str)
    parser.add_argument('indel_maf', metavar='indel_maf', type=str)
    parser.add_argument('titan_output', metavar='titan_output', type=str)
    parser.add_argument('output_dir', metavar='output_dir', type=str)
    parser.add_argument('--handle', metavar='handle', type=str)
    args = parser.parse_args()

    snv_maf = args.snv_maf
    indel_maf = args.indel_maf
    titan_output = args.titan_output
    handle = args.handle
    if not handle:
        # If not analysis handle is provided, simply name the file after the titan output file prefix
        handle = os.path.split(titan_output)[-1].split('.')[0]

    output_dir = args.output_dir

    snv_maf_df = pd.read_csv(snv_maf, delimiter='\t', comment='#', header='infer')
    indel_maf_df = pd.read_csv(indel_maf, delimiter='\t', comment='#', skiprows=3, encoding="ISO-8859-1", header='infer')

    sys.stdout.write("Combining SNV and indel maf file data...\n")
    indel_snv_maf = combine_snv_and_indel_mafs(snv_maf_df, indel_maf_df)

    sys.stdout.write("Loading Titan data...\n")
    titan_df = pd.read_csv(titan_output, delimiter='\t', comment='#', header='infer')

    # Build a data structure that allows for easier searching across the Titan data
    sex = get_sex(snv_maf_df)
    titan_search_tree = build_titan_search_tree(titan_df, sex)

    sys.stdout.write("Merging copy number, SNV, and indel information...\n")
    final_df = add_cn_info_to_indel_snv_maf(indel_snv_maf, titan_search_tree, sex)

    final_output_filepath = os.path.join(output_dir, '{}_pyclone_ready.tsv'.format(handle))
    sys.stdout.write("Writing output to {}...".format(final_output_filepath))

    final_df.to_csv(final_output_filepath, sep='\t',
                    columns=['mutation_id', MafCols.REF_COUNT, MafCols.ALT_COUNT, 'normal_cn', 'major_cn', 'minor_cn', MafCols.CHR, MafCols.START, MafCols.END, MafCols.CLASSIFICATION, MafCols.TYPE],
                    header=['mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'major_cn', 'minor_cn', MafCols.CHR, MafCols.START, MafCols.END, MafCols.CLASSIFICATION, MafCols.TYPE], index=False)

    sys.stdout.write("Done!")


if __name__ == '__main__':
    main()
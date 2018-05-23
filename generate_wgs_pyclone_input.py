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
        --facets_output/titan_output CNAs/GCT001-TP-NT-SM-DPBZL-SM-DPBZK.titan.optimalClust_ploidy2.segs.txt 
        output_folder
        
Generate the tab-delimited file used as input to the build_mutations_file script for PyClone.

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


class FacetsCols:
    CHR = "chrom"
    START = "start"
    END = "end"
    TOTAL_CN = "tcn.em"
    MINOR_CN = "lcn.em"


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
    """Find the point halfway between the start and end of the mutation -- relevant for SNVs and indels"""
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
            return int(2)
        elif sex == 'male':
            return int(1)
    else:
        return int(DEFAULT_CELL_COPY_NUMBER)


def build_chrom_map(titan_df, sex, input_type='titan'):
    """Given a dataframe of Titan data, build a dictionary where the chromosomes are the keys and the values
    are the segments"""
    sys.stdout.write("Building chromosome map for CNA data...\n")
    chr_map = defaultdict(lambda: defaultdict(dict))
    for index, row in titan_df.iterrows():
        if input_type == 'titan':
            chrom = str(row[TitanCols.CHR])
            major_cn = row[TitanCols.MAJOR_CN]
            minor_cn = row[TitanCols.MINOR_CN]
            start = row[TitanCols.START]
            end = row[TitanCols.END]
        elif input_type == 'facets':
            try:
                chrom = str(row[FacetsCols.CHR])
                minor_cn = int(row[FacetsCols.MINOR_CN])
                total_cn = int(row[FacetsCols.TOTAL_CN])
                major_cn = total_cn - minor_cn
                start = row[FacetsCols.START]
                end = row[FacetsCols.END]
            except:
                pass

        if not (chrom == 'Y' and sex == 'female'):
            chr_map[chrom][start] = {"end": end,
                                     "major_cn": major_cn,
                                     "minor_cn": minor_cn,
                                     "normal_cn": total_copy_number_based_on_sex_and_chrom(sex, chrom)
                                     }
    return chr_map


def add_cn_info_to_indel_snv_maf(indel_snv_maf, chrom_map, sex):
    """Add the titan allelic copy number information to the combined snv indel maf"""
    all_chrs = chrom_map.keys()

    def major_and_minor_cns(combined_maf_row):
        # For all rows where the chromosome isn't even in the search tree, set the minor and major alleles to the
        # default, non-aneuploidy values
        chrom = combined_maf_row[MafCols.CHR]
        position = combined_maf_row["Center_Position"]

        if chrom not in all_chrs:
            return int(DEFAULT_MAJOR_COPY_NUMBER), int(DEFAULT_MINOR_COPY_NUMBER), total_copy_number_based_on_sex_and_chrom(
                sex, chrom)
        else:
            starts = list(chrom_map[chrom].keys())
            start_index = np.searchsorted(starts, position)
            relevant_start = starts[start_index - 1]
            seg = chrom_map[chrom][relevant_start]
            if seg.get("end") >= position:
                return int(seg.get("major_cn")), int(seg.get("minor_cn")), int(seg.get("normal_cn"))
            else:
                return int(DEFAULT_MAJOR_COPY_NUMBER), int(DEFAULT_MINOR_COPY_NUMBER), total_copy_number_based_on_sex_and_chrom(
                    sex, chrom)
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
        sys.stdout.write('Y chromosome found in sample SNV file. Determining male sex. \n')
        return 'male'
    else:
        sys.stdout.write('Y chromosome not found in sample SNV. Determining female sex. Will ignore any Y chromosome '
                         'information from other files. \n')
        return 'female'


def main():
    parser = ArgumentParser(description='Generate .tsv input for PyClone from SNV maf,'
                                        ' GATK indel file, and TITAN output')
    parser.add_argument('snv_maf', metavar='snv_maf', type=str)
    parser.add_argument('indel_maf', metavar='indel_maf', type=str)
    parser.add_argument('--titan_output', metavar='titan_output', type=str)
    parser.add_argument('--facets_output', metavar='facets_output', type=str)
    parser.add_argument('output_dir', metavar='output_dir', type=str)
    parser.add_argument('--handle', metavar='handle', type=str)
    args = parser.parse_args()

    snv_maf = args.snv_maf
    indel_maf = args.indel_maf

    if args.titan_output:
        cn_output = args.titan_output
        input_type = 'titan'
    elif args.facets_output:
        cn_output = args.facets_output
        input_type = 'facets'
    else:
        sys.exit('Either a titan or facets output file must be provided with allelic copy number information')

    handle = args.handle
    if not handle:
        # If not analysis handle is provided, simply name the file after the titan output file prefix
        handle = os.path.split(cn_output)[-1].split('.')[0]

    output_dir = args.output_dir

    snv_maf_df = pd.read_csv(snv_maf, delimiter='\t', comment='#', header='infer')

    sex = get_sex(snv_maf_df)
    indel_maf_df = pd.read_csv(indel_maf, delimiter='\t', comment='#', skiprows=3, encoding="ISO-8859-1",
                               header='infer')
    if sex == 'female':
        # if there are no SNVs on the Y chromosome then let's ignore the data on Y chromosome for indels as well
        indel_maf_df = indel_maf_df[indel_maf_df.Chromosome != 'Y']

    sys.stdout.write("Combining SNV and indel maf file data...\n")
    indel_snv_maf = combine_snv_and_indel_mafs(snv_maf_df, indel_maf_df)

    sys.stdout.write("Loading CNA data...\n")
    cn_df = pd.read_csv(cn_output, delimiter='\t', comment='#', header='infer')

    # Build a data structure that allows for easier searching across the CNA data
    chrom_map = build_chrom_map(cn_df, sex, input_type)

    sys.stdout.write("Merging copy number, SNV, and indel information...\n")
    final_df = add_cn_info_to_indel_snv_maf(indel_snv_maf, chrom_map, sex)

    final_output_filepath = os.path.join(output_dir, '{}_pyclone_ready.tsv'.format(handle))
    sys.stdout.write("Writing output to {}...\n".format(final_output_filepath))

    sys.stdout.write("Filtering homozygous deletion sites out of output (major cn == 0)...\n")
    sys.stdout.write("Length before filtering: {}\n".format(len(final_df)))
    final_df = final_df[final_df['major_cn'] > 0]
    sys.stdout.write("Length after filtering: {}\n".format(len(final_df)))

    final_df[MafCols.REF_COUNT] = final_df[MafCols.REF_COUNT].fillna(0)
    final_df[MafCols.REF_COUNT] = final_df[MafCols.REF_COUNT].astype(int)
    final_df[MafCols.ALT_COUNT] = final_df[MafCols.ALT_COUNT].fillna(0)
    final_df[MafCols.ALT_COUNT] = final_df[MafCols.ALT_COUNT].astype(int)

    final_df.to_csv(final_output_filepath, sep='\t',
                    columns=['mutation_id', MafCols.REF_COUNT, MafCols.ALT_COUNT, 'normal_cn', 'major_cn', 'minor_cn',
                             MafCols.CHR, MafCols.START, MafCols.END, MafCols.CLASSIFICATION, MafCols.TYPE],
                    header=['mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'major_cn', 'minor_cn',
                            MafCols.CHR, MafCols.START, MafCols.END, MafCols.CLASSIFICATION, MafCols.TYPE], index=False)

    sys.stdout.write("Done!\n")


if __name__ == '__main__':
    main()
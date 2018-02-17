# PyCloneTSVGenerationWithTitan
Combine snv maf, indel maf, and titan allelic copy number calls to generate the tab-delimited file used as input to
PyClone's build_mutations_file script.

Usage:
    `generate_wgs_pyclone_input.py
        SNVs_OxoGFFPE/GCT001-TP-NT-SM-DPBZL-SM-DPBZK.ffpeBias.maf.annotated
        Indels/GCT001-TP-NT-SM-DPBZL-SM-DPBZK.indel.wgs.maf.annotated
        CNAs/GCT001-TP-NT-SM-DPBZL-SM-DPBZK.titan.optimalClust_ploidy2.segs.txt
        output_folder`

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
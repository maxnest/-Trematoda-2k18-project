"""
Author: Nesterenko Maksim
E-mail: maxnest@ro.ru
Organization: St.Petersburg State University, Russia
Department: Department of Invertebrate zoology
Data: 01.05.2018
Last modification: 07.05.2018
"""

try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

try:
    import numpy
except ImportError:
    print("Please check if module 'numpy' is installed")
    quit()

from operator import itemgetter
from math import *

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True,
                    help="table with correct groups of similar sequences "
                         "['append_annotation_and_select_correct_groups [Two Species]' output file]")
parser.add_argument('--species_1', type=str, required=True,
                    help="first species name tag [for instance: 'Psilo']")
parser.add_argument('--sp1_red_exp_1', type=argparse.FileType('r'), required=True,
                    help="first file (biological replicate) with TPM-values for first species rediae "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--sp1_red_exp_2', type=argparse.FileType('r'), required=True,
                    help="second file (biological replicate) with TPM-values for first species rediae "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--sp1_cer_exp_1', type=argparse.FileType('r'), required=True,
                    help="first file (biological replicate) with TPM-values for first species cercaria "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--sp1_cer_exp_2', type=argparse.FileType('r'), required=True,
                    help="second file (biological replicate) with TPM-values for first species cercaria "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--sp1_mar_exp_1', type=argparse.FileType('r'), required=True,
                    help="first file (biological replicate) with TPM-values for first species marita "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--sp1_mar_exp_2', type=argparse.FileType('r'), required=True,
                    help="second file (biological replicate) with TPM-values for first species marita "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--species_2', type=str, required=True,
                    help="second species name tag [for instance: 'Sphaer']")
parser.add_argument('--sp2_red_exp_1', type=argparse.FileType('r'), required=True,
                    help="first file (biological replicate) with TPM-values for second species rediae "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--sp2_red_exp_2', type=argparse.FileType('r'), required=True,
                    help="second file (biological replicate) with TPM-values for second species rediae "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--sp2_cer_exp_1', type=argparse.FileType('r'), required=True,
                    help="first file (biological replicate) with TPM-values for second species cercaria"
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--sp2_cer_exp_2', type=argparse.FileType('r'), required=True,
                    help="second file (biological replicate) with TPM-values for second species cercaria "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--sp2_mar_exp_1', type=argparse.FileType('r'), required=True,
                    help="first file (biological replicate) with TPM-values for second species marita "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--sp2_mar_exp_2', type=argparse.FileType('r'), required=True,
                    help="second file (biological replicate) with TPM-values for second species marita "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def read_table_with_groups(tab, group_dict, all_values):
    """
    the function reads table with groups of similar sequences and append them to the python dictionary.
    Last one have next structure:
    key = annotation, value for this key = python dictionary,
    where contigs is a keys and values for each of them is other dictionary with TPM (equal to zero at this step)
    @param tab: table with groups of homologs
    @param group_dict: python dictionary
    @param all_values: python list
    """
    head = tab.readline()
    for line in tab:
        group_description = line.strip().split("\t")
        group_name, values = group_description[0], group_description[1:]
        group_dict["{name}".format(name=group_name)] = \
            {value: {"rep_1": 0, "rep_2": 0} for value in values if len(value) > 1}
        all_values.extend([value for value in values if len(value) > 1])


def write_expression(group_dict, file_with_expression, all_values, rep_num):
    """
    The function append TPM values to each contig in dictionary
    @param group_dict: previously created python dictionary
    @param file_with_expression: table with TPM values
    @param all_values: python list with all contigs from 'group_dict'
    @param rep_num: number of biological replicate
    """
    TPM_dict = {}
    head = file_with_expression.readline()
    for line in file_with_expression:
        contig_description = line.strip().split("\t")
        contig_id, TPM = contig_description[0], contig_description[3]
        if contig_id in all_values:
            TPM_dict[contig_id] = float(TPM)

    keys = [key for key in TPM_dict.keys()]
    for group, contigs in group_dict.items():
        for contig in contigs:
            if contig in keys:
                group_dict[group][contig]["rep_{rep_num}".format(rep_num=rep_num)] += TPM_dict[contig]



def write_output(group_dict, tag, species_1, species_2):
    with open("{tag}.tab".format(tag=tag), 'a') as output_file:
        output_file.write("Annotation\t{species_1}_red\t{species_1}_cer\t{species_1}_mar\t"
                          "{species_2}_red\t{species_2}_cer\t{species_2}_mar\n".format(
                           species_1=species_1, species_2=species_2))
        for group, values in group_dict.items():
            species_1_red_value, species_1_cer_value, species_1_mar_value, \
            species_2_red_value, species_2_cer_value, species_2_mar_value = \
                0, 0, 0, 0, 0, 0
            for contig in values:
                if species_1 in contig:
                    if "{species}_red".format(species=species_1) in contig:
                        species_1_red_value += round(numpy.mean([float(group_dict[group][contig]["rep_1"]),
                                                           float(group_dict[group][contig]["rep_2"])]), 2)
                    elif "{species}_cer".format(species=species_1) in contig:
                        species_1_cer_value += round(numpy.mean([float(group_dict[group][contig]["rep_1"]),
                                                           float(group_dict[group][contig]["rep_2"])]), 2)
                    elif "{species}_mar".format(species=species_1) in contig:
                        species_1_mar_value += round(numpy.mean([float(group_dict[group][contig]["rep_1"]),
                                                           float(group_dict[group][contig]["rep_2"])]), 2)
                elif species_2 in contig:
                    if "{species}_red".format(species=species_2) in contig:
                        species_2_red_value += round(numpy.mean([float(group_dict[group][contig]["rep_1"]),
                                                           float(group_dict[group][contig]["rep_2"])]), 2)
                    elif "{species}_cer".format(species=species_2) in contig:
                        species_2_cer_value += round(numpy.mean([float(group_dict[group][contig]["rep_1"]),
                                                           float(group_dict[group][contig]["rep_2"])]), 2)
                    elif "{species}_mar".format(species=species_2) in contig:
                        species_2_mar_value += round(numpy.mean([float(group_dict[group][contig]["rep_1"]),
                                                           float(group_dict[group][contig]["rep_2"])]), 2)
            output_file.write("{ID}\t{sp1_red}\t{sp1_cer}\t{sp1_mar}\t"
                              "{sp2_red}\t{sp2_cer}\t{sp2_mar}\n".format(
                                ID=group, sp1_red=species_1_red_value, sp1_cer=species_1_cer_value,
                                sp1_mar=species_1_mar_value, sp2_red=species_2_red_value,
                                sp2_cer=species_2_cer_value, sp2_mar=species_2_mar_value))


if __name__ == "__main__":
    group_dict, all_values = {}, []
    print("***** Step 1: The script reads table with groups of similar sequences *****")
    read_table_with_groups(args.tab, group_dict, all_values)
    print("***** Step 2: The script reads first file with TPM-values of {species} rediae *****".format(
        species=args.species_1))
    write_expression(group_dict, args.sp1_red_exp_1, all_values, "1")
    print("***** Step 3: The script reads second file with TPM-values of {species} rediae *****".format(
        species=args.species_1))
    write_expression(group_dict, args.sp1_red_exp_2, all_values, "2")
    print("***** Step 4: The script reads first file with TPM-values of {species} cercariae *****".format(
        species=args.species_1))
    write_expression(group_dict, args.sp1_cer_exp_1, all_values, "1")
    print("***** Step 5: The script reads second file with TPM-values of {species} cercariae *****".format(
        species=args.species_1))
    write_expression(group_dict, args.sp1_cer_exp_2, all_values, "2")
    print("***** Step 6: The script reads first file with TPM-values of {species} marita *****".format(
        species=args.species_1))
    write_expression(group_dict, args.sp1_mar_exp_1, all_values, "1")
    print("***** Step 7: The script reads second file with TPM-values of {species} marita *****".format(
        species=args.species_1))
    write_expression(group_dict, args.sp1_mar_exp_2, all_values, "2")
    print("***** Step 8: The script reads first file with TPM-values of {species} rediae *****".format(
        species=args.species_2))
    write_expression(group_dict, args.sp2_red_exp_1, all_values, "1")
    print("***** Step 9: The script reads second file with TPM-values of {species} rediae *****".format(
        species=args.species_2))
    write_expression(group_dict, args.sp2_red_exp_2, all_values, "2")
    print("***** Step 10: The script reads first file with TPM-values of {species} cercariae *****".format(
        species=args.species_2))
    write_expression(group_dict, args.sp2_cer_exp_1, all_values, "1")
    print("***** Step 11: The script reads second file with TPM-values of {species} cercariae *****".format(
        species=args.species_2))
    write_expression(group_dict, args.sp2_cer_exp_2, all_values, "2")
    print("***** Step 12: The script reads first file with TPM-values of {species} marita *****".format(
        species=args.species_2))
    write_expression(group_dict, args.sp2_mar_exp_1, all_values, "1")
    print("***** Step 13: The script reads second file with TPM-values of {species} marita *****".format(
        species=args.species_2))
    write_expression(group_dict, args.sp2_mar_exp_2, all_values, "2")
    print("***** Step 14: {tag} output file creating *****".format(tag=args.tag))
    write_output(group_dict, args.tag, args.species_1, args.species_2)

"""
Author: Nesterenko Maksim
E-mail: maxnest@ro.ru
Organization: St.Petersburg State University, Russia
Department: Department of Invertebrate zoology
Data: 01.05.2018
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

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True,
                    help="table with correct groups of homologs or orthologs "
                         "['append_annotation_and_select_correct_groups.py' output file]")
parser.add_argument('--s1_first_exp', type=argparse.FileType('r'), required=True,
                    help="first file (biological replicate) with TPM-values for first sample "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--s1_second_exp', type=argparse.FileType('r'), required=True,
                    help="second file (biological replicate) with TPM-values for first sample "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--first_tag', type=str, required=True,
                    help="tag for first sample [for instance: 'Psilo_cer']")
parser.add_argument('--s2_first_exp', type=argparse.FileType('r'), required=True,
                    help="first file (biological replicate) with TPM-values for second sample "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--s2_second_exp', type=argparse.FileType('r'), required=True,
                    help="second file (biological replicate) with TPM-values for second sample "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--second_tag', type=str, required=True,
                    help="tag for second sample [for instance: 'Psilo_red']")
parser.add_argument('-t', '--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def read_table_with_groups(tab, group_dict, all_values):
    """
    The function reads table with groups of similar sequences and append them to the python dictionary.
    Last one have next structure:
    key = annotation, value for this key = python dictionary,
    where contigs is a keys and values for each of them is other dictionary with TPM (equal to zero at this step)
    @param tab: table with groups of similar sequences
    @param group_dict: python dictionary
    @param all_values: python list
    """
    head = tab.readline()
    for line in tab:
        group_description = line.strip().split("\t")
        group_name, values = group_description[0], group_description[1:]
        group_dict["{name}".format(name=group_name)] = \
            {value: {"rep_1": 0, "rep_2": 0} for value in values if len(value) > 1}
        all_values.extend(values)


def write_expression(group_dict, file_with_expression, all_values, rep_num):
    """
    The function append TPM values to each contig in dictionary
    @param group_dict: previously created python dictionary
    @param file_with_expression: table with TPM values
    @param all_values: python list with all contigs from 'group_dict'
    @param rep_num: number of biological replicate
    """
    head = file_with_expression.readline()
    for line in file_with_expression:
        contig_description = line.strip().split("\t")
        contig_id, TPM = contig_description[0], contig_description[3]
        if contig_id in all_values:
            for group, contigs in group_dict.items():
                if contig_id in contigs:
                    group_dict[group][contig_id]["rep_{rep_num}".format(rep_num=rep_num)] += float(TPM)
                    break


def write_output(group_dict, tag, first_tag, second_tag):
    with open("{tag}.tab".format(tag=tag), 'a') as output_file:
        output_file.write("Annotation\t{sample_1}\t{sample_2}\t\n".format(
            sample_1=first_tag, sample_2=second_tag
        ))
        for group, values in group_dict.items():
            sample_1, sample_2 = 0, 0
            for contig in values:
                if first_tag in contig:
                    sample_1 += round(numpy.mean([float(group_dict[group][contig]["rep_1"]),
                                                  float(group_dict[group][contig]["rep_2"])]), 2)
                elif second_tag in contig:
                    sample_2 += round(numpy.mean([float(group_dict[group][contig]["rep_1"]),
                                                  float(group_dict[group][contig]["rep_2"])]), 2)
            output_file.write("{ID}\t{sample_1}\t{sample_2}\n".format(ID=group, sample_1=sample_1, sample_2=sample_2))

if __name__ == "__main__":
    group_dict, all_values = {}, []
    print("***** Step 1: The script reads table with groups of similar sequences *****")
    read_table_with_groups(args.tab, group_dict, all_values)
    print("***** Step 2: The script reads first file with TPM-values for {sp1} *****".format(sp1=args.first_tag))
    write_expression(group_dict, args.s1_first_exp, all_values, "1")
    print("***** Step 3: The script reads second file with TPM-values for {sp1} *****".format(sp1=args.first_tag))
    write_expression(group_dict, args.s1_second_exp, all_values, "2")
    print("***** Step 4: The script reads first file with TPM-values for {sp2} *****".format(sp2=args.second_tag))
    write_expression(group_dict, args.s2_first_exp, all_values, "1")
    print("***** Step 5: The script reads second file with TPM-values for {sp2} *****".format(sp2=args.second_tag))
    write_expression(group_dict, args.s2_second_exp, all_values, "2")
    print("***** Step 6: {tag} output file creating *****".format(tag=args.tag))
    write_output(group_dict, args.tag, args.first_tag, args.second_tag)

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

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True,
                    help="table with groups of similar sequences "
                         "[Transcriptologs.py (Ambrosino and Chiusano, 2017) output file]")
parser.add_argument('--sample_1_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for first sample "
                         "[FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--sample_1_tag', type=str, required=True, help="tag for first sample [for instance 'Psilo_ref']")
parser.add_argument('--sample_2_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for second sample "
                         "[FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--sample_2_tag', type=str, required=True, help="tag for second sample [for instance 'Sphaer_ref']")
parser.add_argument('-t', '--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def read_table_with_groups(tab, group_dict, all_values):
    """
    The function reads table with groups of ortho|homologs contigs and write them to new python dictionary.
    Last one have next structure: key = group name;
    value for this key = python dictionary, where contigs is a keys and values for each of them is empty list.
    @param tab: table with groups of similar contigs
    @param group_dict: python dictionary
    @param all_values: python list
    """
    head = tab.readline()
    number = 1
    for line in tab:
        group_description = line.strip().split("\t")
        seq_1, seq_2 = group_description[0], group_description[1]
        group_dict["group_{num}".format(num=number)] = {seq_1: [], seq_2: []}
        number += 1
        all_values.extend([seq_1, seq_2])


def read_and_append_annotation(group_dict, file_with_annotation, all_values):
    """
    The function reads table with sample annotation
    and appends ID of potential protein in python dictionary with group of similar contigs
    @param group_dict: previously created python dictionary with groups
    @param file_with_annotation: table with annotation for transcriptome
    @param all_values: previously created python list with IDs of all contigs from dictionary with groups
    """
    head = file_with_annotation.readline()
    for line in file_with_annotation:
        contig_description = line.strip().split("\t")
        contig_id, best_hit = contig_description[0], contig_description[2]
        if contig_id in all_values:
            for group, contigs in group_dict.items():
                if contig_id in contigs:
                    group_dict[group][contig_id].append(best_hit)


def check_groups(group_dict, new_dict):
    """
    The function checks each group and select only groups where all contig have same annotation
    @param group_dict: python dictionary with groups and annotation
    @param new_dict: python dictionary with groups of correct ortho|homologs
    """
    unknown_hit = 1
    for group, contigs in group_dict.items():
        annotations = []
        for contig in contigs:
            annotations.append(group_dict[group][contig][0])
        if len(set(annotations)) == 1:
            if annotations[0] != "-":
                new_dict[annotations[0]] = {contig for contig in contigs}
            else:
                new_dict["unknown_hit_{num}".format(num=unknown_hit)] = {contig for contig in contigs}
                unknown_hit += 1


def write_output(new_dict, sample_1_tag, sample_2_tag, tag):
    with open("{tag}.tab".format(tag=tag), 'a') as output:
        output.write("Group_annotation\t{Sample_1_tag}\t{Sample_2_tag}\n".format(
                      Sample_1_tag=sample_1_tag, Sample_2_tag=sample_2_tag))
        for annotation, values in new_dict.items():
            sample_1_contig, sample_2_contig = "", ""
            for value in values:
                if sample_1_tag in value:
                    sample_1_contig += value
                elif sample_2_tag in value:
                    sample_2_contig += value
            output.write("{annotation}\t{sample_1_contig}\t{sample_2_contig}\n".format(
                          annotation=annotation, sample_1_contig=sample_1_contig, sample_2_contig=sample_2_contig))

if __name__ == "__main__":
    group_dict, new_dict, all_values = {}, {}, []
    print("***** Step 1: The script parses table with ortho|homologs groups *****")
    read_table_with_groups(args.tab, group_dict, all_values)
    print("***** Step 2: The script reads file with annotation for {sample_1} *****".format(sample_1=args.sample_1_tag))
    read_and_append_annotation(group_dict, args.sample_1_anno, all_values)
    print("***** Step 3: The script reads file with annotation for {sample_2} *****".format(sample_2=args.sample_2_tag))
    read_and_append_annotation(group_dict, args.sample_2_anno, all_values)
    print("***** Step 4: The script check annotations in groups of similar sequences*****")
    check_groups(group_dict, new_dict)
    print("***** Step 5: {tag} output file creating *****".format(tag=args.tag))
    write_output(new_dict, args.sample_1_tag, args.sample_2_tag, args.tag)

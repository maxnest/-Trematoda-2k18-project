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
                    help="table with groups of similar contigs "
                         "['contig_groups_recovery [Two species].py ' output file]")
parser.add_argument('--Pr_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for Psilotrema simillimum rediae transcriptome"
                         " [FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--Pc_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for P.simillium cercariae transcriptome"
                         " [FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--Pm_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for P.simillium marita transcriptome"
                         " [FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--Sr_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for Sphaeridiotrema pseudoglobulus rediae transcriptome"
                         " [FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--Sc_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for S.pseudoglobulus cercariae transcriptome"
                         " [FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--Sm_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for S.pseudoglobulus marita transcriptome"
                         " [FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('-t', '--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def read_table_with_groups(tab, group_dict, all_values):
    """
    The function reads table with groups of similar sequences and write them to new python dictionary.
    Last one have next structure: key = group name;
    value for this key = python dictionary, where contigs is a keys and values for each of them is empty list.
    @param tab: table with groups of similar contigs
    @param group_dict: python dictionary
    @param all_values: python list
    """
    head = tab.readline()
    for line in tab:
        group_description = line.strip().split("\t")
        group_name, values = group_description[0], group_description[1:]
        group_dict[group_name] = {value: [] for value in values if len(value) > 1}
        all_values.extend(values)


def read_and_append_annotation(group_dict, file_with_annotation, all_values):
    """
    The function reads table with transcriptome annotation
    and appends ID of potential protein in python dictionary with group of similar contigs
    @param group_dict: previously created python dictionary with groups
    @param file_with_annotation: table with transcriptome annotation
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
    @param group_dict: previously created python dictionary with groups
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


def write_output(new_dict, tag):
    with open("{tag}.tab".format(tag=tag), 'a') as output:
        output.write('Group_annotation\tPsilo_red\tPsilo_cer\tPsilo_mar\tSphaer_red\tSphaer_cer\tSphaer_mar\n')
        for group_key, contigs_values in new_dict.items():
            Psilo_red_contig, Psilo_cer_contig, Psilo_mar_contig, \
            Sphaer_red_contig, Sphaer_cer_contig, Sphaer_mar_contig = "", "", "", "", "", ""
            for contig in contigs_values:
                if "Psilo" in contig:
                    if "_red" in contig:
                        Psilo_red_contig += contig
                    elif "_cer" in contig:
                        Psilo_cer_contig += contig
                    elif "_mar" in contig:
                        Psilo_mar_contig += contig
                elif "Sphaer" in contig:
                    if "_red" in contig:
                        Sphaer_red_contig += contig
                    elif "_cer" in contig:
                        Sphaer_cer_contig += contig
                    elif "_mar" in contig:
                        Sphaer_mar_contig += contig
            output.write("{group_name}\t{Psilo_red}\t{Psilo_cer}\t{Psilo_mar}\t"
                         "{Sphaer_red}\t{Sphaer_cer}\t{Sphaer_mar}\n".format(
                          group_name=group_key, Psilo_red=Psilo_red_contig, Psilo_cer=Psilo_cer_contig,
                          Psilo_mar=Psilo_mar_contig, Sphaer_red=Sphaer_red_contig, Sphaer_cer=Sphaer_cer_contig,
                          Sphaer_mar=Sphaer_mar_contig))

if __name__ == "__main__":
    group_dict = {}
    new_dict = {}
    all_values = []
    print("***** Step 1: The script parses table with ortho|homologs groups *****")
    read_table_with_groups(args.tab, group_dict, all_values)
    print("***** Step 2: The script reads file with annotation for Psilo_red *****")
    read_and_append_annotation(group_dict, args.Pr_anno, all_values)
    print("***** Step 3: The script reads file with annotation for Psilo_cer *****")
    read_and_append_annotation(group_dict, args.Pc_anno, all_values)
    print("***** Step 4: The script reads file with annotation for Psilo_mar *****")
    read_and_append_annotation(group_dict, args.Pm_anno, all_values)
    print("***** Step 5: The script reads file with annotation for Sphaer_red *****")
    read_and_append_annotation(group_dict, args.Sr_anno, all_values)
    print("***** Step 6: The script reads file with annotation for Sphaer_cer *****")
    read_and_append_annotation(group_dict, args.Sc_anno, all_values)
    print("***** Step 7: The script reads file with annotation for Sphaer_mar *****")
    read_and_append_annotation(group_dict, args.Sm_anno, all_values)
    print("***** Step 8: The script checks annotations in groups *****")
    check_groups(group_dict, new_dict)
    print("***** Step 9: {tag} output file creating *****".format(tag=args.tag))
    write_output(new_dict, args.tag)

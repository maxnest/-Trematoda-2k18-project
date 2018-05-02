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
                    help="table with groups of orthologs or homologs "
                         "[Transcriptologs.py (Ambrosino and Chiusano, 2017) output file "
                         "with annotation (after 'append_annotation_and_select_correct_groups [For two samples]')]")
parser.add_argument('--sample_1_domains', type=argparse.FileType('r'), required=True,
                    help="table with possible domains architectures of predicted proteins from first sample"
                         "['hmmscan3_parser' output file]")
parser.add_argument('--sample_1_tag', type=str, required=True,
                    help="tag for first sample [for instance 'Psilo_cer']")
parser.add_argument('--sample_2_domains', type=argparse.FileType('r'), required=True,
                    help="table with possible domains architectures of predicted proteins from second sample"
                         "['hmmscan3_parser' output file]")
parser.add_argument('--sample_2_tag', type=str, required=True,
                    help="tag for second sample [for instance 'Psilo_red']")
parser.add_argument('-t', '--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def read_table_with_groups(tab, group_dict, all_values):
    """
    The function reads table with groups of ortho|homologs contigs and write them to new python dictionary.
    Last one have next structure: key = group name (annotation);
    value for this key = python dictionary, where contigs is a keys and values for each of them is empty list.
    @param tab: table with groups of similar contigs
    @param group_dict: python dictionary
    @param all_values: python list
    """
    head = tab.readline()
    for line in tab:
        group_description = line.strip().split("\t")
        annotation, seq_1, seq_2 = group_description[0], group_description[1], group_description[2]
        group_dict[annotation] = {seq_1: [], seq_2: []}
        all_values.extend([seq_1, seq_2])


def read_table_with_domains(table_with_domains, group_dict, all_values):
    """
    The function reads table with possible domain architectures of predicted protein
    @param table_with_domains: input file with domains
    @param group_dict: previously created python dictionary with groups of similar sequences
    @param all_values: python list
    """
    head = table_with_domains.readline()
    for line in table_with_domains:
        contig_description = line.strip().split("\t")
        contig_ID, domains = contig_description[0][:-3], [domain for domain in contig_description[2].split(" ")]
        if contig_ID in all_values:
            for annotation, contigs in group_dict.items():
                for contig in contigs:
                    if contig == contig_ID:
                        group_dict[annotation][contig].extend(domains)


def comparing_lists_of_domains(group_dict, without_splicing, with_splicing, without_common):
    """
    The function provide comparison of similar sequences domain architectures
    @param group_dict: previously created python dictionary with groups of similar sequences
    @param without_splicing: python list
    @param with_splicing: python list
    @param without_common: python list
    """
    for annotation, contigs in group_dict.items():
        contigs_keys = [contig for contig in contigs]

        first_contig, second_contig = contigs_keys[0], contigs_keys[1]
        domains_in_first, domains_in_second = group_dict[annotation][first_contig], \
                                              group_dict[annotation][second_contig]
        common_domains = set(domains_in_first) & set(domains_in_second)

        if len(domains_in_first) != 0 or len(domains_in_second) != 0:
            if len(common_domains) != 0:    # common domains exist!
                common_in_first, common_in_second = 0, 0
                for common_domain in common_domains:
                    if domains_in_first.count(common_domain) == domains_in_second.count(common_domain):
                        # copy number of common domains are equal
                        common_in_first += domains_in_first.count(common_domain)
                        common_in_second += domains_in_second.count(common_domain)
                    else:
                        with_splicing.append(annotation)
                        break
                if common_in_first == len(domains_in_first) and common_in_second == len(domains_in_second):
                    # another domains does not exist
                    without_splicing.append(annotation)
                else:
                    with_splicing.append(annotation)
            else:
                without_common.append(annotation)

    print("+++++ The script found: {num_with_splicing} proteins with variation in domain architecture "
          ",{num_without_splicing} proteins without it "
          "and {num_without_common} protein pairs without common domains +++++".format(
           num_with_splicing=len(set(with_splicing)), num_without_splicing=len(set(without_splicing)),
           num_without_common=len(set(without_common))))


def contigs_extraction(contig_1, first_tag, contig_2, second_tag, annotation, group_dict):
    """
    The function extracts contigs from python dictionary and append in correct python lists
    @param contig_1: python list for contig from first sample
    @param first_tag: tag for first sample
    @param contig_2: python list for contig from second sample
    @param second_tag: tag for second sample
    @param annotation: protein ID (python string)
    @param group_dict: previously created python dictionary with groups of similar sequences
    """
    for contig_key in group_dict[annotation].keys():
        if first_tag in contig_key:
            contig_1.clear()
            contig_1.append(contig_key)
        elif second_tag in contig_key:
            contig_2.clear()
            contig_2.append(contig_key)
    return contig_1, contig_2


def write_output(tag, with_splicing, without_splicing, without_common, sample_1_tag, sample_2_tag, group_dict):
    with open("{output_tag}_summary.tab".format(output_tag=tag), 'a') as output:
        output.write("Among similar proteins of {sp_1} and {sp_2} the script found: "
                     "{with_splicing} sequences with variation in domain architecture "
                     ",{without_splicing} without it and {without_common} protein pairs without common domains".format(
                      sp_1=sample_1_tag, sp_2=sample_2_tag, with_splicing=len(set(with_splicing)),
                      without_splicing=len(set(without_splicing)), without_common=len(set(without_common))))

    contig_from_first, contig_from_second = [], []

    if len(with_splicing) != 0:
        with open("{output_tag}_with_variation".format(output_tag=tag), 'a') as with_splicing_output:
            for annotation in set(with_splicing):
                contigs_extraction(contig_from_first, sample_1_tag, contig_from_second, sample_2_tag, annotation, group_dict)
                with_splicing_output.write("{annotation}\t{contig_1}\t{contig_2}\n".format(
                    annotation=annotation, contig_1=contig_from_first[0], contig_2=contig_from_second[0]))

    if len(without_splicing) != 0:
        with open("{output_tag}_without_variation".format(output_tag=tag), 'a') as without_splicing_output:
            for annotation in set(without_splicing):
                contigs_extraction(contig_from_first, sample_1_tag, contig_from_second, sample_2_tag, annotation, group_dict)
                without_splicing_output.write("{annotation}\t{contig_1}\t{contig_2}\n".format(
                    annotation=annotation, contig_1=contig_from_first[0], contig_2=contig_from_second[0]))

    if len(without_common) != 0:
        with open("{output_tag}_without_common".format(output_tag=tag), 'a') as without_common_output:
            for annotation in set(without_common):
                contigs_extraction(contig_from_first, sample_1_tag, contig_from_second, sample_2_tag, annotation, group_dict)
                without_common_output.write("{annotation}\t{contig_1}\t{contig_2}\n".format(
                    annotation=annotation, contig_1=contig_from_first[0], contig_2=contig_from_second[0]))


if __name__ == "__main__":
    group_dict, all_values, with_splicing, without_splicing, without_common = {}, [], [], [], []
    print("***** Step 1: The script reads table with groups *****")
    read_table_with_groups(args.tab, group_dict, all_values)
    print("***** Step 2: The script reads table with domains from {first_sample} *****".format(
        first_sample=args.sample_1_tag))
    read_table_with_domains(args.sample_1_domains, group_dict, all_values)
    print("***** Step 3: The script reads table with domains from {second_sample} *****".format(
        second_sample=args.sample_2_tag))
    read_table_with_domains(args.sample_2_domains, group_dict, all_values)
    print("***** Step 4: Lists of domains comparison *****")
    comparing_lists_of_domains(group_dict, without_splicing, with_splicing, without_common)
    print("***** Step 5: {tag} output file creating *****".format(tag=args.tag))
    write_output(args.tag, with_splicing, without_splicing, without_common,
                 args.sample_1_tag, args.sample_2_tag, group_dict)

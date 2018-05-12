"""
Author: Nesterenko Maksim
E-mail: maxnest@ro.ru
Organization: St.Petersburg State University, Russia
Department: Department of Invertebrate zoology
Data: 10.05.2018
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
                    help="table with groups of similar sequences and normalized TPM values"
                         "['create_table_with_TPM [Two sample]' + 'append_scaling_factors']")
parser.add_argument('--threshold', type=int, required=True,
                    help="the threshold of TPM (transcript per million) value, after which, in user opinion, "
                         "gene expression becomes biologically significant"
                         "[for instance: 5]")
parser.add_argument('--phase_1_tag', type=str, required=True,
                    help="tag for the first phase from the first compared pair [for instance: 'Psilo_cer']")
parser.add_argument('--phase_2_tag', type=str, required=True,
                    help="tag for the second phase from the first compared pair [for instance: 'Psilo_red']")
parser.add_argument('--phase_3_tag', type=str, required=True,
                    help="tag for the first phase from the second compared pair [for instance: 'Sphaer_cer']")
parser.add_argument('--phase_4_tag', type=str, required=True,
                    help="tag for the second phase from the second compared pair [for instance: 'Sphaer_red']")
parser.add_argument('--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def read_table_with_groups(tab, header, group_dict,
                           first_phase, second_phase, third_phase, fourth_phase, threshold, failing):
    """
    The function reads input file with normalized mean TPM-values for similar sequences from two compared samples
    @param tab: table with similar sequences
    @param header: python list
    @param group_dict: python dictionary
    @param first_phase: tag for first phase from 1 pair
    @param second_phase: tag for second phase from 1 pair
    @param third_phase: tag for first phase from 2 pair
    @param fourth_phase: tag for second phase from 2 pair
    @param threshold: user-defined threshold for minimal expression level
    @param threshold: python list
    """
    header.extend(tab.readline().strip().split("\t"))
    first_phase_index, second_phase_index, third_phase_index, fourth_phase_index = \
        header.index(first_phase), header.index(second_phase), header.index(third_phase), header.index(fourth_phase)
    count = 0
    for line in tab:
        group_description = line.strip().split("\t")
        annotation, phase_1_exp, phase_2_exp, phase_3_exp, phase_4_exp = \
            group_description[0], float(group_description[first_phase_index]), \
            float(group_description[second_phase_index]), float(group_description[third_phase_index]), \
            float(group_description[fourth_phase_index])

        if phase_1_exp >= threshold or phase_2_exp >= threshold or phase_3_exp >= threshold or phase_4_exp >= threshold:
            group_dict[annotation] = {"{phase_1}_contig_{num}".format(phase_1=first_phase, num=count): phase_1_exp,
                                      "{phase_2}_contig_{num}".format(phase_2=second_phase, num=count): phase_2_exp,
                                      "{phase_3}_contig_{num}".format(phase_3=third_phase, num=count): phase_3_exp,
                                      "{phase_4}_contig_{num}".format(phase_4=fourth_phase, num=count): phase_4_exp}
            count += 1
        else:
            failing.append(annotation)
    print("!!! {number} genes have expression level less, than custom threshold !!!".format(number=len(failing)))


def jaccard_similarity(one_list, other_list):
    """
    FROM: http://dataconomy.com/2015/04/implementing-the-five-most-popular-similarity-measures-in-python/
    The function measure the similarity between two sets of genes (Jaccard similarity index)
    @param one_list: list of genes
    @param other_list: list of genes
    """
    intersection_cardinality = len(set.intersection(*[set(one_list), set(other_list)]))
    union_cardinality = len(set.union(*[set(one_list), set(other_list)]))
    if union_cardinality != 0:
        return round(intersection_cardinality / float(union_cardinality), 2)
    else:
        return 0


def evaluation_of_TPM_difference(group_dict, first_phase, second_phase, phase_1_specific, phase_2_specific,
                                 first_pair_similar, first_phase_medium, second_phase_medium,
                                 first_phase_different, second_phase_different):
    """
    @param group_dict: previously created python dictionary
    @param first_phase: tag for first phase
    @param second_phase: tag for second phase
    @param phase_1_specific: python list for genes that expressed only in first phase
    @param phase_2_specific: python list for genes that expressed only in second phase
    @param first_pair_similar: python list for genes with similar expression level in first and second phases
    @param first_phase_medium: python list for genes:  2 < first_phase_TPM/second_phase_TPM < 4
    @param second_phase_medium: python list for genes: 2 < second_phase_TPM/first_phase_TPM < 4
    @param first_phase_different: python list for genes: first_phase_TPM/second_phase_TPM > 4
    @param second_phase_different: python list for genes: second_phase_TPM/first_phase_TPM > 4
    """
    for annotation, contigs in group_dict.items():
        first_contig, second_contig = [contig for contig in contigs.keys() if contig.startswith(first_phase)], \
                                      [contig for contig in contigs.keys() if contig.startswith(second_phase)]
        first_phase_exp, second_phase_exp = group_dict[annotation][first_contig[0]], \
                                            group_dict[annotation][second_contig[0]]
        if first_phase_exp != 0 and second_phase_exp == 0:
            phase_1_specific.append(annotation)
        elif first_phase_exp == 0 and second_phase_exp != 0:
            phase_2_specific.append(annotation)
        elif first_phase_exp != 0 and second_phase_exp != 0:
            difference = max([first_phase_exp, second_phase_exp]) / min([first_phase_exp, second_phase_exp])
            if first_phase_exp > second_phase_exp:
                if difference > 4:
                    first_phase_different.append(annotation)
                elif difference < 2:
                    first_pair_similar.append(annotation)
                else:
                    first_phase_medium.append(annotation)
            elif second_phase_exp > first_phase_exp:
                if difference > 4:
                    second_phase_different.append(annotation)
                elif difference < 2:
                    first_pair_similar.append(annotation)
                else:
                    second_phase_medium.append(annotation)


def comparison(first_pair_similar, second_pair_similar, first_phase_different, first_phase_medium,
               second_phase_different, second_phase_medium, third_phase_different, third_phase_medium,
               fourth_phase_different, fourth_phase_medium, dict_with_metrics):
    dict_with_metrics["Similar"] = jaccard_similarity(first_pair_similar, second_pair_similar)
    dict_with_metrics["Different: First_phase_vs_Third_phase"] = \
        jaccard_similarity(first_phase_different, third_phase_different)
    dict_with_metrics["Medium: First_phase_vs_Third_phase"] = \
        jaccard_similarity(first_phase_medium, third_phase_medium)
    dict_with_metrics["Different: Second_phase_vs_Fourth_phase"] = \
        jaccard_similarity(second_phase_different, fourth_phase_different)
    dict_with_metrics["Medium: Second_phase_vs_Fourth_phase"] = \
        jaccard_similarity(second_phase_medium, fourth_phase_medium)


def write_summary(tag, dict_with_metrics, first_phase, second_phase, third_phase, fourth_phase, threshold,
                  first_specific, second_specific, third_specific, fourth_specific):
    with open("{output_tag}_summary.txt".format(output_tag=tag), 'a') as summary:
        summary.write("[Input]   Comparison: {First_phase}-{Second_phase} versus {Third_phase}-{Fourth_phase}\n"
                      "[Options] Only genes pairs in which at least one gene have TPM >= {Threshold} were considered\n"
                      "*************************************************************************************************\n"
                      "[Metrics] 'Jaccard of Similar' [0 - 1] - \n"
                      "           evaluation of the similarity of gene lists,\n"
                      "           the expression of which in pairs of samples differs less than 2 times\n"
                      "           ({First_phase}/{Second_phase} < 2 and {Third_phase}/{Fourth_phase} < 2)\n"
                      "[Metrics] 'Jaccard of Medium: {First_phase} and {Third_phase}' [0 - 1] - \n"
                      "           evaluation of the similarity of gene lists,\n"
                      "           the expression of which in pairs of samples differs\n"
                      "           more than 2 times, but less than 4 times\n"
                      "           ((2 < {First_phase}/{Second_phase} < 4) and (2 < {Third_phase}/{Fourth_phase} < 4))\n"
                      "[Metrics] 'Jaccard of Different: {First_phase} and {Third_phase}' [0 - 1] - \n"
                      "           evaluation of the similarity of gene lists,\n"
                      "           the expression of which in pairs of samples differs more than 4 times\n"
                      "           ({First_phase}/{Second_phase} > 4 and {Third_phase}/{Fourth_phase} > 4)\n"
                      "[Metrics] 'Jaccard of Medium: {Second_phase} and {Fourth_phase}' [0 - 1] - \n"
                      "           evaluation of the similarity of gene lists,\n"
                      "           the expession of which in pairs of samples differes\n"
                      "           more than 2 times, but less than 4 times\n"
                      "           ((2 < {Second_phase}/{First_phase} < 4) and (2 < {Fourth_phase}/{Third_phase} < 4))\n"
                      "[Metrics] 'Jaccard of Different: {Second_phase} and {Fourth_phase}' [0 - 1] - \n"
                      "           evaluation of the the similarity of gene lists,\n"
                      "           the expression of which in pairs of samples differs more than 4 times\n"
                      "           ({Second_phase}/{First_phase} > 4 and {Fourth_phase}/{Third_phase} > 4)\n"
                      "*************************************************************************************************\n"
                      "[RESULTS] 'Expressed only in {First_phase}' : {first_specific}\n"
                      "[RESULTS] 'Expressed only in {Second_phase}' : {second_specific}\n"
                      "[RESULTS] 'Expressed only in {Third_phase}' : {third_specific}\n"
                      "[RESULTS] 'Expressed only in {Fourth_phase}' : {fourth_specific}\n"
                      "[RESULTS] 'Jaccard of Similar' : {Jaccard_of_similar}\n"
                      "[RESULTS] 'Jaccard of Medium: {First_phase} and {Third_phase}' : {Medium_first_third}\n"
                      "[RESULTS] 'Jaccard of Different: {First_phase} and {Third_phase}' : {Different_first_third}\n"
                      "[RESULTS] 'Jaccard of Medium: {Second_phase} and {Fourth_phase}' : {Medium_second_fourth}\n"
                      "[RESULTS] 'Jaccard of Different: {Second_phase} and {Fourth_phase}' : {Different_second_fourth}\n"
                      "*************************************************************************************************\n".format(
            First_phase=first_phase, Second_phase=second_phase, Third_phase=third_phase, Fourth_phase=fourth_phase,
            Threshold=threshold, first_specific=len(first_specific), second_specific=len(second_specific),
            third_specific=len(third_specific), fourth_specific=len(fourth_specific),
            Jaccard_of_similar=dict_with_metrics["Similar"],
            Medium_first_third=dict_with_metrics["Medium: First_phase_vs_Third_phase"],
            Different_first_third=dict_with_metrics["Different: First_phase_vs_Third_phase"],
            Medium_second_fourth=dict_with_metrics["Medium: Second_phase_vs_Fourth_phase"],
            Different_second_fourth=dict_with_metrics["Different: Second_phase_vs_Fourth_phase"])
        )


def write_in_file(output_name, one_list, other_list):
    with open(output_name, 'a') as output:
        common = set(one_list) & set(other_list)
        for annotation in common:
            output.write("{annotation}\n".format(annotation=annotation))

def write_lists(first_phase, second_phase, third_phase, fourth_phase,
                first_similar, first_medium, second_medium, first_different, second_different,
                second_similar, third_medium, fourth_medium, third_different, fourth_different):
    write_in_file("Similar_intersection.txt", first_similar, second_similar)
    write_in_file("{First}_and_{Third}_medium_intersection.txt".format(First=first_phase, Third=third_phase),
                  first_medium, third_medium)
    write_in_file("{First}_and_{Third}_different_intersection.txt".format(First=first_phase, Third=third_phase),
                  first_different, third_different)
    write_in_file("{Second}_and_{Fourth}_medium_intersection.txt".format(Second=second_phase, Fourth=fourth_phase),
                  second_medium, fourth_medium)
    write_in_file("{Second}_and_{Fourth}_different_intersection.txt".format(Second=second_phase, Fourth=fourth_phase),
                  second_different, fourth_different)


if __name__ == "__main__":
    group_dict, header, failing, first_similar, second_similar, first_medium, second_medium, third_medium, \
    fourth_medium, first_different, second_different, third_different, fourth_different = \
        {}, [], [], [], [], [], [], [], [], [], [], [], []
    first_specific, second_specific, third_specific, fourth_specific = [], [], [], []
    metrics = {}
    print("***** Step 1: The script reads table *****")
    read_table_with_groups(args.tab, header, group_dict, args.phase_1_tag, args.phase_2_tag, args.phase_3_tag,
                           args.phase_4_tag, args.threshold, failing)
    print("***** Step 2: TPM comparison *****")
    evaluation_of_TPM_difference(group_dict, args.phase_1_tag, args.phase_2_tag, first_specific, second_specific,
                                 first_similar, first_medium, second_medium, first_different, second_different)
    evaluation_of_TPM_difference(group_dict, args.phase_3_tag, args.phase_4_tag, third_specific, fourth_specific,
                                 second_similar, third_medium, fourth_medium, third_different, fourth_different)
    print("***** Step 3: Comparison *****")
    comparison(first_similar, second_similar, first_different, first_medium, second_different, second_medium,
               third_different, third_medium, fourth_different, fourth_medium, metrics)
    print("***** Step 4: The script writes summary in {tag}_summary.txt *****".format(tag=args.tag))
    write_summary(args.tag, metrics, args.phase_1_tag, args.phase_2_tag, args.phase_3_tag, args.phase_4_tag, args.threshold,
                  first_specific, second_specific, third_specific, fourth_specific)
    print("***** Step 5: The script writes sets of genes in different files *****")
    write_lists(args.phase_1_tag, args.phase_2_tag, args.phase_3_tag, args.phase_4_tag, first_similar,
                first_medium, second_medium, first_different, second_different, second_similar,
                third_medium, fourth_medium, third_different, fourth_different)

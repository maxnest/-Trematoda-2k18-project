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

from operator import itemgetter
from math import *



parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True,
                    help="table with groups of similar sequences and TPM values for two sample (with replications)"
                         "['create_table_with_TPM [Two sample]' + 'append_scaling_factors']")
parser.add_argument('--upper', type=int, required=True,
                    help="upper threshold (in percent) for separation based on genes expression ranks "
                         "[for instance: if '15' percent were selected, "
                         "15 percent of genes with highest ranks will be selected as 'highly expressed gene set'")
parser.add_argument('--lower', type=int, required=True,
                    help="lower threshold (in percent) for separation based on genes expression ranks "
                         "[for instance: if '15' percent were selected, "
                         "15 percent of genes with lowest ranks (with non-zero expression) "
                         "will be selected as 'lowly expressed gene set'")
parser.add_argument('--threshold', type=int, required=True,
                    help="the threshold of TPM (transcript per million) value, after which, in user opinion, "
                         "gene expression becomes biologically significant"
                         "[for instance: 5]")
parser.add_argument('--wobbling', type=int, required=True,
                    help="the number of position up- and downstream, "
                         "that the user considers it permissible to take into account "
                         "in comparison of ranks of genes in different transcriptomes "
                         "[for instance: user specifed '1', this means that  "
                         "if gene have rank 'x' in one sample and 'x-1', 'x' or 'x+1' in another "
                         "the script assume that ranks of considered gene are similar in both transcriptomes]")
parser.add_argument('--varying', type=argparse.FileType('r'),
                    help="file with IDs of predicted proteins with variation in domain architectures")
parser.add_argument('--sample_1_tag', type=str, required=True,
                    help="tag for first sample [for instance 'Psilo_cer']")
parser.add_argument('--sample_2_tag', type=str, required=True,
                    help="tag for second sample [for instance 'Psilo_red']")
parser.add_argument('--rep_num', type=int, required=True,
                    help="number of biological|technical replicates")
parser.add_argument('--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def read_table_with_groups(tab, group_dict, first_tag, second_tag, threshold):
    """
    The function reads input file with normalized TPM-values for similar sequences from two compared samples
    @param tab: table with similar sequences
    @param group_dict: python dictionary
    @param first_tag: tag for first sample
    @param second_tag: tag for second sample
    @param threshold: user-defined threshold for minimal expression level
    """
    head = tab.readline()
    count = 1
    less_than_threshold = 0
    for line in tab:
        group_description = line.strip().split("\t")
        annotation, average_sample_1_exp, average_sample_2_exp = \
            group_description[0], numpy.mean([float(el) for el in group_description[1:3]]), \
            numpy.mean([float(el) for el in group_description[3:5]])
        # how to take into account number of replicates ?
        if average_sample_1_exp >= threshold or average_sample_2_exp >= threshold:  # or 'and' ?
            group_dict[annotation] = {"{sample_1}_contig_{num}".format(sample_1=first_tag, num=count):
                                      {"exp": average_sample_1_exp, "rank": 0, "variation": 0},
                                      "{sample_2}_contig_{num}".format(sample_2=second_tag, num=count):
                                      {"exp": average_sample_2_exp, "rank": 0, "variation": 0}}
            count += 1
        else:
            less_than_threshold += 1
    print("!!! In {number} pairs of genes expression level less, than custom threshold !!!".format(
        number=less_than_threshold))


def sort_in_order(group_dict, sample_tag):
    """
    The functions sort genes in decreasing order (based on expression level)
    @param group_dict: previously created python dictionary with groups of similar sequences
    @param sample_tag: taf for sample
    """
    list_with_dicts = []

    for annotation, contigs in group_dict.items():
        for contig in contigs.keys():
            if sample_tag in contig:
                    list_with_dicts.append({"annotation": annotation, "contig": contig,
                                            "expression": group_dict[annotation][contig]["exp"]})

    order = sorted(list_with_dicts, key=lambda k: itemgetter("expression")(k), reverse=True)

    for dict in order:
        annotation, contig = dict["annotation"], dict["contig"]
        group_dict[annotation][contig]["rank"] = order.index(dict) + 1


def transcriptome_separation(group_dict, upper, lower, HEG_dict, LEG_dict, CORE_dict, sample_tag):
    """
    the function divides transcriptome for 3 sets: highly expressed genes (HEG),
                                                   lowly expressed genes (LEG) and others (CORE).
    Separation is based on the ranks.
    @param group_dict: python dictionary with groups of similar sequences
    @param upper: threshold for falling into the HEG
    @param lower: threshold for falling into the LEG
    @param HEG_dict: python dictionary
    @param LEG_dict: python dictionary
    @param CORE_dict: python dictionary
    @param sample_tag: tag for sample
    """
    # HEG = highly expressed genes
    # LEG = lowly expressed genes
    percent = len(group_dict.keys()) / 100
    highly, lowly = int(upper * percent), int((100 - lower) * percent)
    for annotation, contigs in group_dict.items():
        for contig in contigs:
            if sample_tag in contig:
                rank = group_dict[annotation][contig]["rank"]
                if rank <= highly:
                    HEG_dict[annotation] = rank
                elif rank >= lowly:
                    LEG_dict[annotation] = rank
                else:
                    CORE_dict[annotation] = rank


def jaccard_similarity(one_dict, other_dict):
    """
    FROM: http://dataconomy.com/2015/04/implementing-the-five-most-popular-similarity-measures-in-python/
    The function measure the similarity between two sets of genes (Jaccard similarity index)
    @param one_dict: set of genes (python dictionary)
    @param other_dict: set of genes (python dictionary)
    """
    intersection_cardinality = len(set.intersection(*[set(one_dict), set(other_dict)]))
    union_cardinality = len(set.union(*[set(one_dict), set(other_dict)]))
    if union_cardinality != 0:
        return round(intersection_cardinality / float(union_cardinality), 2)
    else:
        return 0


def ranks_similarity(one_dict, other_dict, wobbling):
    """
    The function measure similarity of ranks of similar sequences from different transcriptomes.
    Not only full coincidence is taken into account,
    but also small deviations (+ or - user-defined number of rank) (wobbling)
    @param one_dict: set of genes (python dictionary)
    @param other_dict: set of genes (python dictionary)
    @param wobbling: number
    """
    common_annotation = set(one_dict.keys()) & set(other_dict.keys())
    ranks_matches = 0
    for annotation in common_annotation:
        ### wobbling like soft borders ###
        if one_dict[annotation] in [other_el for other_el in
                                    range(other_dict[annotation] - wobbling, other_dict[annotation] + (wobbling + 1))] \
                and other_dict[annotation] in [one_el for one_el in range(one_dict[annotation] - wobbling,
                                                                          one_dict[annotation] + (wobbling + 1))]:

                    ranks_matches += 1
    if len(common_annotation) != 0:
        return round((ranks_matches/len(common_annotation)) * 100, 2)
    else:
        return 0


def proportion_of_varying(dict, varying):
    """
    The function estimates the proportion of proteins with variation in domain architecture
    from total number of such proteins in the set
    @param dict: set of genes (python dictionary)
    @param varying: set of protein with variation (python dictionary)
    """
    varying_in_dict = set(dict.keys()) & set(varying.keys())
    return round(len(varying_in_dict) / len(varying.keys()) * 100, 2)


def lists_comparing(one_dict, other_dict, varying, dict_with_metrics, wobbling):
    """
    The function collects results of comparison of sets
    @param one_dict: set of genes (python dictionary)
    @param other_dict: set of genes (python dictionary)
    @param varying: set of protein with varition (python dictionary)
    @param dict_with_metrics: python dictionary
    @param wobbling: number
    """
    # match between dictionaries (Jaccard)
    # match between ranks of genes (Ranks)
    # percent of proteins with variation in dictionaries
    common_annotation = set(one_dict.keys()) & set(other_dict.keys())
    dict_with_metrics["N_of_common_annotation"] = len(common_annotation)  # number
    dict_with_metrics["Varying_in_first"] = proportion_of_varying(one_dict, varying)  # percent
    dict_with_metrics["Jaccard_common_first"] = jaccard_similarity(common_annotation, one_dict.keys())  # [0 - 1]
    dict_with_metrics["Varying_in_second"] = proportion_of_varying(other_dict, varying)  # percent
    dict_with_metrics["Jaccard_common_second"] = jaccard_similarity(common_annotation, other_dict.keys())  # [0 - 1]
    dict_with_metrics["Jaccard_between"] = jaccard_similarity(one_dict.keys(), other_dict.keys())  # [0 - 1]
    dict_with_metrics["Ranks"] = ranks_similarity(one_dict, other_dict, wobbling)  # percent


def write_summary(tag, HEG_metrics, LEG_metrics, CORE_metrics, sample_1_tag, sample_2_tag, upper, lower, threshold):
    with open("{output_tag}_summary.txt".format(output_tag=tag), 'a') as summary:
        summary.write("[Input]   Table with {First_sample} and {Second_sample} transcriptomes\n"
                      "[Options] Only genes pairs in which at least one gene have TPM >= {Threshold} were considered\n"
                      "[Options] Rank variation ('wobbling') equal to {wobbling} "
                      "          is acceptable in rank similarity evaluation\n"
                      "[Options] 'Highly expressed genes' [HEG] = {Upper} percent of genes with highest ranks\n"
                      "[Options] 'Lowly expressed genes' [LEG] = {Lower} percent of genes with lowest ranks\n"
                      "[Options] 'Core genes' [CORE] = genes not included in HEG or LEG\n"
                      "********************************************************************************************\n"  
                      "[Metrics] 'N of common' [number] - number of proteins IDs detected both in "
                      "           HEG, LEG or CORE of {First_sample} and {Second_sample}\n"
                      "[Metrics] 'Varying in HEG|LEG|CORE of {First_sample}|{Second_sample}' [0-100%] - "
                      "           the percentage of proteins with a variation in domain architecture in set"
                      "           from the total number of such proteins\n"
                      "[Metrics] 'Jaccard [Common - HEG|LEG|CORE of {First_sample}|{Second_sample}]' [0-1] - "
                      "           similarity between sets of common proteins and HEG|LEG|CORE of one sample, respectively\n"
                      "[Metrics] 'Jaccard [between samples]' [0-1] - "
                      "           similarity between {First_sample} and {Second_sample} HEG|LEG|CORE sets\n"
                      "[Metrics] 'Coincidence of ranks' [0-100%] - "
                      "           the percentage of common proteins with same rank in both sample\n"
                      "********************************************************************************************\n"  
                      "[HEG-results] 'N_of_common' : {Common_in_HEG}\n"
                      "[HEG-results] 'Jaccard [between samples]' : {Jaccard_between_HEG}\n"
                      "[HEG-results] 'Coincidence of ranks of common' : {Coincidence_in_HEG}%\n"
                      "[HEG-results] 'Jaccard [Common <-> HEG of {First_sample}]' : {Jaccard_common_HEG_First}\n"
                      "[HEG-results] 'Jaccard [Common <-> HEG of {Second_sample}]' : {Jaccard_common_HEG_Second}\n"
                      "[HEG-results] 'Varying in HEG of {First_sample} : {Varying_in_HEG_of_First}%\n"
                      "[HEG-results] 'Varying in HEG of {Second_sample} : {Varying_in_HEG_of_Second}%\n"
                      "********************************************************************************************\n"
                      "[CORE-results] 'N_of_common' : {Common_in_CORE}\n"
                      "[CORE-results] 'Jaccard [between samples]' : {Jaccard_between_CORE}\n"
                      "[CORE-results] 'Coincidence of ranks of common' : {Coincidence_in_CORE}%\n"
                      "[CORE-results] 'Jaccard [Common <-> CORE of {First_sample}]' : {Jaccard_common_CORE_First}\n"
                      "[CORE-results] 'Jaccard [Common <-> CORE of {Second_sample}]' : {Jaccard_common_CORE_Second}\n"
                      "[CORE-results] 'Varying in CORE of {First_sample} : {Varying_in_CORE_of_First}%\n"
                      "[CORE-results] 'Varying in CORE of {Second_sample} : {Varying_in_CORE_of_Second}%\n"
                      "********************************************************************************************\n"
                      "[LEG-results] 'N_of_common' : {Common_in_LEG}\n"
                      "[LEG-results] 'Jaccard [between samples]' : {Jaccard_between_LEG}\n"
                      "[LEG-results] 'Coincidence of ranks of common' : {Coincidence_in_LEG}%\n"
                      "[LEG-results] 'Jaccard [Common <-> LEG of {First_sample}]' : {Jaccard_common_LEG_First}\n"
                      "[LEG-results] 'Jaccard [Common <-> LEG of {Second_sample}]' : {Jaccard_common_LEG_Second}\n"
                      "[LEG-results] 'Varying in LEG of {First_sample} : {Varying_in_LEG_of_First}%\n"
                      "[LEG-results] 'Varying in LEG of {Second_sample} : {Varying_in_LEG_of_Second}%\n"
                      "********************************************************************************************\n".format(
            First_sample=sample_1_tag, Second_sample=sample_2_tag, Threshold=threshold, wobbling=args.wobbling,
            Upper=upper, Lower=lower, Common_in_HEG=HEG_metrics["N_of_common_annotation"],
            Jaccard_between_HEG=HEG_metrics["Jaccard_between"],
            Coincidence_in_HEG=HEG_metrics["Ranks"], Jaccard_common_HEG_First=HEG_metrics["Jaccard_common_first"],
            Jaccard_common_HEG_Second=HEG_metrics["Jaccard_common_second"],
            Varying_in_HEG_of_First=HEG_metrics["Varying_in_first"],
            Varying_in_HEG_of_Second=HEG_metrics["Varying_in_second"],
            Common_in_CORE=CORE_metrics["N_of_common_annotation"], Jaccard_between_CORE=CORE_metrics["Jaccard_between"],
            Coincidence_in_CORE=CORE_metrics["Ranks"], Jaccard_common_CORE_First=CORE_metrics["Jaccard_common_first"],
            Jaccard_common_CORE_Second=CORE_metrics["Jaccard_common_second"],
            Varying_in_CORE_of_First=CORE_metrics["Varying_in_first"],
            Varying_in_CORE_of_Second=CORE_metrics["Varying_in_second"],
            Common_in_LEG=LEG_metrics["N_of_common_annotation"], Jaccard_between_LEG=LEG_metrics["Jaccard_between"],
            Coincidence_in_LEG=LEG_metrics["Ranks"], Jaccard_common_LEG_First=LEG_metrics["Jaccard_common_first"],
            Jaccard_common_LEG_Second=LEG_metrics["Jaccard_common_second"],
            Varying_in_LEG_of_First=LEG_metrics["Varying_in_first"],
            Varying_in_LEG_of_Second=LEG_metrics["Varying_in_second"])
            )


def write_lists(sample_tag, HEG_dict, LEG_dict, CORE_dict):
    with open("{sample_tag}_HEG.txt".format(sample_tag=sample_tag), 'a') as HEG_output:
        for annotation_key in HEG_dict.keys():
            HEG_output.write("{annotation}\n".format(annotation=annotation_key))

    with open("{sample_tag}_CORE.txt".format(sample_tag=sample_tag), 'a') as CORE_output:
        for annotation_key in CORE_dict.keys():
            CORE_output.write("{annotation}\n".format(annotation=annotation_key))

    with open("{sample_tag}_LEG.txt".format(sample_tag=sample_tag), 'a') as LEG_output:
        for annotation_key in LEG_dict.keys():
            LEG_output.write("{annotation}\n".format(annotation=annotation_key))


if __name__ == "__main__":
    group_dict, HEG_first, HEG_second, CORE_first, CORE_second, LEG_first, LEG_second = {}, {}, {}, {}, {}, {}, {}
    varying_dict = {}
    HEG_metrics, CORE_metrics, LEG_metrics = {}, {}, {}
    print("***** Step 1: The script reads file with ID of proteins with variation *****")
    for line in args.varying:
        description = line.strip().split("\t")
        annotation, contig_1, contig_2 = description[0], description[1], description[2]
        varying_dict[annotation] = [contig_1, contig_2]
    print("***** Step 2: The script reads table with TPM-values *****")
    read_table_with_groups(args.tab, group_dict, args.sample_1_tag, args.sample_2_tag, args.threshold)
    print("***** Step 3: The script assign rank for sequences from {sample_1} *****".format(sample_1=args.sample_1_tag))
    sort_in_order(group_dict, args.sample_1_tag)
    print("***** Step 4: The script assign rank for sequences from {sample_2} *****".format(sample_2=args.sample_2_tag))
    sort_in_order(group_dict, args.sample_2_tag)
    print("***** Step 5: The script divides transcriptomes *****")
    transcriptome_separation(group_dict, args.upper, args.lower, HEG_first, LEG_first, CORE_first, args.sample_1_tag)
    transcriptome_separation(group_dict, args.upper, args.lower, HEG_second, LEG_second, CORE_second, args.sample_2_tag)
    print("***** Step 6: Comparison of lists with highly expressed genes *****")
    lists_comparing(HEG_first, HEG_second, varying_dict, HEG_metrics, args.wobbling)
    print("***** Step 7: Comparison of lists with 'core' genes *****")
    lists_comparing(CORE_first, CORE_second, varying_dict, CORE_metrics, args.wobbling)
    print("***** Step 8: Comparison of lists with lowly expressed genes *****")
    lists_comparing(LEG_first, LEG_second, varying_dict, LEG_metrics, args.wobbling)
    print("***** Step 9: The script writes summary in {tag}_summary.txt *****".format(tag=args.tag))
    write_summary(args.tag, HEG_metrics, LEG_metrics, CORE_metrics, args.sample_1_tag, args.sample_2_tag, args.upper,
                  args.lower, args.threshold)
    print("***** Step 10: The script writes results of transcriptome division in separate files *****")
    write_lists(args.sample_1_tag, HEG_first, LEG_first, CORE_first)
    write_lists(args.sample_2_tag, HEG_second, LEG_second, CORE_second)
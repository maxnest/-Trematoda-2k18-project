"""
Author: Nesterenko Maksim
E-mail: maxnest@ro.ru
Organization: St.Petersburg State University, Russia
Department: Department of Invertebrate zoology
Data: 09.05.2018
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
parser.add_argument('--varying', type=argparse.FileType('r'),
                    help="file with IDs of predicted proteins with variation in domain architectures")
parser.add_argument('--sample_1_tag', type=str, required=True,
                    help="tag for first sample [for instance 'Psilo_cer']")
parser.add_argument('--sample_2_tag', type=str, required=True,
                    help="tag for second sample [for instance 'Psilo_red']")
parser.add_argument('--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def read_table_with_groups(tab, group_dict, first_tag, second_tag, threshold, failing):
    """
    The function reads input file with normalized mean TPM-values for similar sequences from two compared samples
    @param tab: table with similar sequences
    @param group_dict: python dictionary
    @param first_tag: tag for first sample
    @param second_tag: tag for second sample
    @param threshold: user-defined threshold for minimal expression level
    """
    head = tab.readline()
    count = 1
    for line in tab:
        group_description = line.strip().split("\t")
        annotation, average_sample_1_exp, average_sample_2_exp = \
            group_description[0], float(group_description[1]), float(group_description[2])
        if average_sample_1_exp >= threshold or average_sample_2_exp >= threshold:
            group_dict[annotation] = {"{sample_1}_contig_{num}".format(sample_1=first_tag, num=count):
                                      {"exp": average_sample_1_exp, "rank": 0},
                                      "{sample_2}_contig_{num}".format(sample_2=second_tag, num=count):
                                      {"exp": average_sample_2_exp, "rank": 0}}
            count += 1
        else:
            failing.append(annotation)
    print("!!! In {number} pairs of genes expression level less, than custom threshold !!!".format(
           number=len(failing)))


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


def proportion_of_varying(dict, varying):
    """
    The function estimates the proportion of proteins with variation in domain architecture
    from total number of such proteins in the set
    @param dict: set of genes (python dictionary)
    @param varying: set of protein with variation (python dictionary)
    """
    varying_in_dict = set(dict) & set(varying.keys())
    return round(len(varying_in_dict) / len(dict) * 100, 2)


def evaluation_of_TPM_difference(group_dict, first_tag, second_tag, first_specific, second_specific,
                                 similar, medium, different):
    for annotation, value in group_dict.items():
        first_TPM, second_TPM = 0, 0
        for key in value.keys():
            if first_tag in key:
                first_TPM += group_dict[annotation][key]["exp"]
            elif second_tag in key:
                second_TPM += group_dict[annotation][key]["exp"]

        if first_TPM == 0:
            first_specific.append(annotation)
        elif second_TPM == 0:
            second_specific.append(annotation)
        else:
            difference = max([first_TPM, second_TPM]) / min([first_TPM, second_TPM])
            if difference > 4:
                different.append(annotation)
            elif difference < 2:   # threshold
                similar.append(annotation)
            else:
                medium.append(annotation)


def lists_comparing(group_dict, HEG_first, HEG_second, CORE_first, CORE_second, LEG_first, LEG_second,
                    varying, similar, medium, different, dict_with_metrics):
    # match between dictionaries (Jaccard)
    # percent of proteins with variation in dictionaries
    common_expressed = []
    for annotation, value in group_dict.items():
        contigs = [contig for contig in value]
        if group_dict[annotation][contigs[0]]["exp"] > 0 and group_dict[annotation][contigs[1]]["exp"] > 0:
            common_expressed.append(annotation)
    dict_with_metrics["Jaccard:HEG"] = jaccard_similarity(HEG_first.keys(), HEG_second.keys())
    dict_with_metrics["Jaccard:CORE"] = jaccard_similarity(CORE_first.keys(), CORE_second.keys())
    dict_with_metrics["Jaccard:LEG"] = jaccard_similarity(LEG_first.keys(), LEG_second.keys())
    dict_with_metrics["N_with_similar"] = len(similar)
    dict_with_metrics["N_with_medium"] = len(medium)
    dict_with_metrics["N_with_different"] = len(different)
    dict_with_metrics["Similar-Common"] = round((len(similar)/len(common_expressed)) * 100, 2)  # percent
    dict_with_metrics["Medium-Common"] = round((len(medium)/len(common_expressed)) * 100, 2)  # percent
    dict_with_metrics["Different-Common"] = round((len(different)/len(common_expressed)) * 100, 2)  # percent
    dict_with_metrics["Varying in Similar"] = proportion_of_varying(similar, varying)
    dict_with_metrics["Varying in Medium"] = proportion_of_varying(medium, varying)
    dict_with_metrics["Varying in Different"] = proportion_of_varying(different, varying)


def write_summary(tag, dict_with_metrics, sample_1_tag, sample_2_tag, upper, lower, threshold, failing,
                  first_specific, second_specific):
    with open("{output_tag}_summary.txt".format(output_tag=tag), 'a') as summary:
        summary.write("[Input]   Table with {First_sample} and {Second_sample} transcriptomes\n"
                      "[Options] Only genes pairs in which at least one gene have TPM >= {Threshold} were considered\n"
                      "[Options] 'Highly expressed genes' [HEG] = {Upper} percent of genes with highest ranks\n"
                      "[Options] 'Lowly expressed genes' [LEG] = {Lower} percent of genes with lowest ranks\n"
                      "[Options] 'Core genes' [CORE] = genes not included in HEG or LEG\n"
                      "*************************************************************************************************\n"
                      "[Metrics] 'Jaccard: HEG' [0 - 1] - \n"
                      "           similarity between HEG sets of {First_sample} and {Second_sample}\n"
                      "[Metrics] 'Jaccard: CORE' [0 - 1] - \n"
                      "           similarity between CORE sets of {First_sample} and {Second_sample}\n"
                      "[Metrics] 'Jaccard: LEG' [0 - 1] - \n"
                      "           similarity between LEG sets of {First_sample} and {Second_sample}\n"
                      "[Metrics] 'N with similar' [number] - number of genes, \n"
                      "           for which difference in expression less than two times is observed\n"
                      "[Metrics] 'N with medium' [number] - "
                      "           number of genes, for which difference in expression \n"
                      "           more than two times and less than four times is observed\n"
                      "[Metrics] 'N with different' [number] - number of genes, \n"
                      "           for which difference in expression more than four time\n"
                      "[Metrics] 'Similar/Common' - [0 - 100%] - "
                      "           the percentage of genes simultaneously working in both sample,\n"
                      "           with difference in expression levels less than two times\n"
                      "[Metrics] 'Medium/Common' - [0 - 100%] - \n"
                      "           the percentage of genes simultaneously working in both samples,\n"
                      "           with difference in expression levels more than two times and less than four times\n"
                      "[Metrics] 'Different/Common' - [0 - 100%] - \n"
                      "           the percentage of genes simultaneously working in both samples,\n"
                      "           with difference in expression levels more than four times\n"
                      "[Metrics] 'Varying in Similar' - [0 - 100%] - \n"
                      "           the percentage of proteins with a variation in domain architecture\n"
                      "           in set of genes with similar level of expression in both samples\n"
                      "[Metrics] 'Varying in Medium' - [0 - 100%] - \n"
                      "           the percentage of proteins with a variation in domain architecture\n"
                      "           in set of genes with medium level of difference in expression between both samples\n"
                      "[Metrics] 'Varying in Different' - [0 - 100%] - \n"
                      "           the percentage of proteins with a variation in domain architecture\n"
                      "           in set of genes with very different level of expression in both samples\n"
                      "[Metrics] 'Failing' [number] - \n"
                      "           number of pairs of genes which expression level less than custom threshold\n"
                      "[Metrics] 'Expressed only in {First_sample}|{Second_sample}' [number] - \n"
                      "           number of genes with expression only in one sample \n"
                      "*************************************************************************************************\n"
                      "[RESULTS] 'Expressed only in {First_sample}' : {first_specific}\n"
                      "[RESULTS] 'Expressed only in {Second_sample}' : {second_specific}\n"
                      "[RESULTS] 'Failing' : {failing}\n"
                      "[RESULTS] 'Jaccard: HEG' : {Jaccard_HEG}\n"
                      "[RESULTS] 'Jaccard: CORE' : {Jaccard_CORE}\n"
                      "[RESULTS] 'Jaccard: LEG' : {Jaccard_LEG}\n"
                      "[RESULTS] 'N with similar' : {N_with_similar}\n"
                      "[RESULTS] 'N with medium' : {N_with_medium}\n"
                      "[RESULTS] 'N with different' : {N_with_different}\n"
                      "[RESULTS] 'Similar/Common' : {Similar_Common}%\n"
                      "[RESULTS] 'Medium/Common' : {Medium_Common}%\n"
                      "[RESULTS] 'Different/Common' : {Different_Common}%\n"
                      "[RESULTS] 'Varying in Similar' : {Varying_in_similar}%\n"
                      "[RESULTS] 'Varying in Medium' : {Varying_in_medium}%\n"
                      "[RESULTS] 'Varying in Different' : {Varying_in_different}%\n"
                      "*************************************************************************************************\n".format(
            First_sample=sample_1_tag, Second_sample=sample_2_tag, Threshold=threshold, Upper=upper, Lower=lower,
            first_specific=len(first_specific), second_specific=len(second_specific), failing=len(failing),
            Jaccard_HEG=dict_with_metrics["Jaccard:HEG"], Jaccard_CORE=dict_with_metrics["Jaccard:CORE"],
            Jaccard_LEG=dict_with_metrics["Jaccard:LEG"], N_with_similar=dict_with_metrics["N_with_similar"],
            N_with_medium=dict_with_metrics["N_with_medium"], N_with_different=dict_with_metrics["N_with_different"],
            Similar_Common=dict_with_metrics["Similar-Common"], Medium_Common=dict_with_metrics["Medium-Common"],
            Different_Common=dict_with_metrics["Different-Common"], Varying_in_similar=dict_with_metrics["Varying in Similar"],
            Varying_in_medium=dict_with_metrics["Varying in Medium"],
            Varying_in_different=dict_with_metrics["Varying in Different"])
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
    varying_dict, dict_with_metrics, failing, first_specific, second_specific = {}, {}, [], [], []
    similar, medium, different = [], [], []
    print("***** Step 1: The script reads file with ID of proteins with variation *****")
    for line in args.varying:
        description = line.strip().split("\t")
        annotation, contig_1, contig_2 = description[0], description[1], description[2]
        varying_dict[annotation] = [contig_1, contig_2]
    print("***** Step 2: The script reads table with TPM-values *****")
    read_table_with_groups(args.tab, group_dict, args.sample_1_tag, args.sample_2_tag, args.threshold, failing)
    print("***** Step 3: The script assign rank for sequences from {sample_1} *****".format(sample_1=args.sample_1_tag))
    sort_in_order(group_dict, args.sample_1_tag)
    print("***** Step 4: The script assign rank for sequences from {sample_2} *****".format(sample_2=args.sample_2_tag))
    sort_in_order(group_dict, args.sample_2_tag)
    print("***** Step 5: The script divides transcriptomes *****")
    transcriptome_separation(group_dict, args.upper, args.lower, HEG_first, LEG_first, CORE_first, args.sample_1_tag)
    transcriptome_separation(group_dict, args.upper, args.lower, HEG_second, LEG_second, CORE_second, args.sample_2_tag)
    print("***** Step 6: Evaluation of TPM difference *****")
    evaluation_of_TPM_difference(group_dict, args.sample_1_tag, args.sample_2_tag, first_specific, second_specific,
                                 similar, medium, different)
    print("***** Step 7: Comparison ******")
    lists_comparing(group_dict, HEG_first, HEG_second, CORE_first, CORE_second, LEG_first, LEG_second, varying_dict,
                    similar, medium, different, dict_with_metrics)
    print("***** Step 8: The script writes summary in {tag}_summary.txt *****".format(tag=args.tag))
    write_summary(args.tag, dict_with_metrics, args.sample_1_tag, args.sample_2_tag, args.upper, args.lower,
                  args.threshold, failing, first_specific, second_specific)
    print("***** Step 9: The script writes results of transcriptome division in separate files *****")
    write_lists(args.sample_1_tag, HEG_first, LEG_first, CORE_first)
    write_lists(args.sample_2_tag, HEG_second, LEG_second, CORE_second)

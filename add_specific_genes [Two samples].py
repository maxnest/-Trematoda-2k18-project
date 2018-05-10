"""
Author: Nesterenko Maksim
E-mail: maxnest@ro.ru
Organization: St.Petersburg State University, Russia
Department: Department of Invertebrate zoology
Data: 7.05.2018
"""

try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True,
                    help="table with correct groups of homologs or orthologs"
                         "['append_annotation_and_select_correct_groups' output file]")
parser.add_argument('--s1_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for first sample"
                         "[FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--s1_scores', type=argparse.FileType('r'), required=True,
                    help="file with scores for contigs from first sample"
                         "[TransRate (Smith-Unna et al., 2016) output file: 'contig.csv']")
parser.add_argument('--s1_tag', type=str, required=True,
                    help="first sample name tag ['Psilo_red', for instance]")
parser.add_argument('--s2_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for second sample"
                         "[FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--s2_scores', type=argparse.FileType('r'), required=True,
                    help="file with scores for contigs from second sample"
                         "[TransRate (Smith-Unna et al., 2016) output file: 'contig.csv']")
parser.add_argument('--s2_tag', type=str, required=True,
                    help="second sample name tag ['Sphaer_red', for instance]")
parser.add_argument('--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def read_table_with_groups(tab, group_dict, all_values):
    """
    The function reads table with groups of homologs and orthologs and append them to the python dictionary.
    @param tab: table with groups of similar sequences
    @param group_dict: python dictionary
    @param all_values: python list
    """
    head = tab.readline()
    for line in tab:
        group_description = line.strip().split("\t")
        group_name, values = group_description[0], group_description[1:]
        group_dict[group_name] = [contig for contig in values if len(contig) > 1]
        all_values.extend(group_dict[group_name])


def adding_additional_contigs(group_dict, file_with_annotation, sample_tag, all_values):
    """
    The function adds additional contigs in previously created python dictionary with groups of similar sequences
    """
    head = file_with_annotation.readline()
    contig_without_hits = 1
    for line in file_with_annotation:
        contig_description = line.strip().split("\t")
        contig_id, best_hit = contig_description[0], contig_description[2]
        if contig_id not in all_values:
            if best_hit == "-":
                name = "{sample_tag}_unknown_hit_{num}".format(sample_tag=sample_tag, num=contig_without_hits)
                contig_without_hits += 1
                group_dict[name] = [contig_id]
            else:
                name = "{hit}/{sample_tag}".format(hit=best_hit, sample_tag=sample_tag)
                if name not in group_dict.keys():
                    group_dict[name] = [contig_id]
                else:
                    group_dict[name].append(contig_id)


def contig_selection(group_dict, file_with_scores, sample_tag):
    """
    If several contigs have the same annotation, the function selects only one of them by means of TransRate score
    @param group_dict: previously created python dictionary
    @param file_with_scores: table with contigs TransRate scores
    @param phase: life cycle phase name tag
    @param all_values: python list with all contigs from 'group_dict'
    """
    dict_with_repeats = {}
    contigs = []
    name = "/{sample_tag}".format(sample_tag=sample_tag)
    for key, value in group_dict.items():
        if key.endswith(name) and len(value) > 1:
            dict_with_repeats[key] = {contig: 0 for contig in value}
            contigs.extend(value)

    # parsing csv file with TransRate scores for contigs
    head_csv = file_with_scores.readline()
    for line in file_with_scores:
        contig_metrics = line.strip().split(',')
        contig_id, score = contig_metrics[0], contig_metrics[8]
        if contig_id in contigs:
            for key, value in dict_with_repeats.items():
                if contig_id in value:
                    dict_with_repeats[key][contig_id] += float(score)

    # selection
    for key in dict_with_repeats.keys():
        for contig_key, score in dict_with_repeats[key].items():
            if score == max([score for score in dict_with_repeats[key].values()]):
                group_dict[key] = [contig_key]


def write_output(group_dict, tag, species_1_tag, species_2_tag):
    with open("{output_tag}.tab".format(output_tag=tag), 'a') as output:
        output.write('Annotation\t{s1_tag}\t{s2_tag}\n'.format(
            s1_tag=species_1_tag, s2_tag=species_2_tag))
        for group_key, contigs_values in group_dict.items():
            s1_contig, s2_contig = "", ""
            for contig in contigs_values:
                if species_1_tag in contig:
                    s1_contig += contig
                elif species_2_tag in contig:
                    s2_contig += contig
            output.write("{group_name}\t{s1_contig}\t{s2_contig}\n".format(
                          group_name=group_key, s1_contig=s1_contig, s2_contig=s2_contig))


if __name__ == "__main__":
    group_dict, all_values = {}, []
    print("***** Step 1: The script reads table with ortho|homologs groups *****")
    read_table_with_groups(args.tab, group_dict, all_values)
    print("***** Step 2: The script checks file with annotation for {sample} *****".format(sample=args.s1_tag))
    adding_additional_contigs(group_dict, args.s1_anno, args.s1_tag, all_values)
    print("***** Step 3: The script checks file with annotation for {sample} *****".format(sample=args.s2_tag))
    adding_additional_contigs(group_dict, args.s2_anno, args.s2_tag, all_values)
    print("***** Step 4: Selection of additional {sample} contigs based on TransRate-score *****".format(
        sample=args.s1_tag))
    contig_selection(group_dict, args.s1_scores, args.s1_tag)
    print("***** Step 5: Selection of additional {sample} contigs based on TransRate-score *****".format(
        sample=args.s2_tag))
    contig_selection(group_dict, args.s2_scores, args.s2_tag)
    print("***** Step 6: {tag} output file creating *****".format(tag=args.tag))
    write_output(group_dict, args.tag, args.s1_tag, args.s2_tag)

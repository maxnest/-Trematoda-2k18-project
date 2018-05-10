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
                    help="table with correct groups of homologs and orthologs"
                         "['append_annotation_and_select_correct_groups' output file]")
parser.add_argument('--sp1_red_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for first species rediae transcriptome "
                         "[FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--sp1_red_scores', type=argparse.FileType('r'), required=True,
                    help="file with scores for contigs from first species rediae transcriptome "
                         "[TransRate (Smith-Unna et al., 2016) output file: 'contig.csv']")
parser.add_argument('--sp1_cer_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for first species cercaria transcriptome "
                         "[FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--sp1_cer_scores', type=argparse.FileType('r'), required=True,
                    help="file with scores for contigs from first species cercaria transcriptome "
                         "[TransRate (Smith-Unna et al., 2016) output file: 'contig.csv']")
parser.add_argument('--sp1_mar_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for first species marita transcriptome "
                         "[FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--sp1_mar_scores', type=argparse.FileType('r'), required=True,
                    help="file with scores for contigs from first species marita transcriptome "
                         "[TransRate (Smith-Unna et al., 2016) output file: 'contig.csv']")
parser.add_argument('--sp1_tag', type=str, required=True,
                    help="first species name tag ['Psilo', for instance]")
parser.add_argument('--sp2_red_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for second species rediae transcriptome "
                         "[FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--sp2_red_scores', type=argparse.FileType('r'), required=True,
                    help="file with scores for contigs from second species rediae transcriptome "
                         "[TransRate (Smith-Unna et al., 2016) output file: 'contig.csv']")
parser.add_argument('--sp2_cer_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for second species cercaria transcriptome "
                         "[FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--sp2_cer_scores', type=argparse.FileType('r'), required=True,
                    help="file with scores for contigs from second species cercaria transcriptome "
                         "[TransRate (Smith-Unna et al., 2016) output file: 'contig.csv']")
parser.add_argument('--sp2_mar_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for second species marita transcriptome "
                         "[FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--sp2_mar_scores', type=argparse.FileType('r'), required=True,
                    help="file with scores for contigs from second species marita transcriptome "
                         "[TransRate (Smith-Unna et al., 2016) output file: 'contig.csv']")
parser.add_argument('--sp2_tag', type=str, required=True,
                    help="second species name tag ['Sphaer', for instance]")
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


def adding_additional_contigs(group_dict, file_with_annotation, species, phase, all_values):
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
                name = "{species}_{phase}_unknown_hit_{num}".format(species=species,
                                                                    phase=phase, num=contig_without_hits)
                contig_without_hits += 1
                group_dict[name] = [contig_id]
            else:
                name = "{hit}/{species}_{phase}".format(hit=best_hit, species=species, phase=phase)
                if name not in group_dict.keys():
                    group_dict[name] = [contig_id]
                else:
                    group_dict[name].append(contig_id)


def contig_selection(group_dict, file_with_scores, species, phase):
    """
    If several contigs have the same annotation, the function selects only one of them by means of TransRate score
    @param group_dict: previously created python dictionary
    @param file_with_scores: table with contigs TransRate scores
    @param phase: life cycle phase name tag
    @param all_values: python list with all contigs from 'group_dict'
    """
    dict_with_repeats = {}
    contigs = []
    name = "/{species}_{phase}".format(species=species, phase=phase)
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
        output.write('Annotation\t{sp1}_red\t{sp1}_cer\t{sp1}_mar\t{sp2}_red\t{sp2}_cer\t{sp2}_mar\n'.format(
            sp1=species_1_tag, sp2=species_2_tag))
        for group_key, contigs_values in group_dict.items():
            sp1_red_contig, sp1_cer_contig, sp1_mar_contig, \
            sp2_red_contig, sp2_cer_contig, sp2_mar_contig = "", "", "", "", "", ""
            for contig in contigs_values:
                if species_1_tag in contig:
                    if "_red" in contig:
                        sp1_red_contig += contig
                    elif "_cer" in contig:
                        sp1_cer_contig += contig
                    elif "_mar" in contig:
                        sp1_mar_contig += contig
                elif species_2_tag in contig:
                    if "_red" in contig:
                        sp2_red_contig += contig
                    elif "_cer" in contig:
                        sp2_cer_contig += contig
                    elif "_mar" in contig:
                        sp2_mar_contig += contig
            output.write("{group_name}\t{sp1_red}\t{sp1_cer}\t{sp1_mar}\t"
                         "{sp2_red}\t{sp2_cer}\t{sp2_mar}\n".format(
                          group_name=group_key, sp1_red=sp1_red_contig, sp1_cer=sp1_cer_contig,
                          sp1_mar=sp1_mar_contig, sp2_red=sp2_red_contig, sp2_cer=sp2_cer_contig,
                          sp2_mar=sp2_mar_contig))


if __name__ == "__main__":
    group_dict, all_values = {}, []
    print("***** Step 1: The script reads table with ortho|homologs groups *****")
    read_table_with_groups(args.tab, group_dict, all_values)
    print("***** Step 2: The script checks file with annotation for {species} rediae *****".format(species=args.sp1_tag))
    adding_additional_contigs(group_dict, args.sp1_red_anno, args.sp1_tag, "red", all_values)
    print("***** Step 3: The script checks file with annotation for {species} cercaria *****".format(species=args.sp1_tag))
    adding_additional_contigs(group_dict, args.sp1_cer_anno, args.sp1_tag, "cer", all_values)
    print("***** Step 4: The script checks file with annotation for {species} marita *****".format(species=args.sp1_tag))
    adding_additional_contigs(group_dict, args.sp1_mar_anno, args.sp1_tag, "mar", all_values)
    print("***** Step 5: The script checks file with annotation for {species} rediae *****".format(species=args.sp2_tag))
    adding_additional_contigs(group_dict, args.sp2_red_anno, args.sp2_tag, "red", all_values)
    print("***** Step 6: The script checks file with annotation for {species} cercaria *****".format(species=args.sp2_tag))
    adding_additional_contigs(group_dict, args.sp2_cer_anno, args.sp2_tag, "cer", all_values)
    print("***** Step 7: The script checks file with annotation for {species} marita *****".format(species=args.sp2_tag))
    adding_additional_contigs(group_dict, args.sp2_mar_anno, args.sp2_tag, "mar", all_values)
    print("***** Step 8: Selection of additional {species} rediae contigs based on TransRate-scores *****".format(
        species=args.sp1_tag))
    contig_selection(group_dict, args.sp1_red_scores, args.sp1_tag, "red")
    print("***** Step 9: Selection of additional {species} cercaria contigs based on TransRate-scores *****".format(
        species=args.sp1_tag))
    contig_selection(group_dict, args.sp1_cer_scores, args.sp1_tag, "cer")
    print("***** Step 10: Selection of additional {species} marita contigs based on TransRate-scores *****".format(
        species=args.sp1_tag))
    contig_selection(group_dict, args.sp1_mar_scores, args.sp1_tag, "mar")
    print("***** Step 11: Selection of additional {species} rediae contigs based on TransRate-scores *****".format(
        species=args.sp2_tag))
    contig_selection(group_dict, args.sp2_red_scores, args.sp2_tag, "red")
    print("***** Step 12: Selection of additional {species} cercaria contigs based on TransRate-scores *****".format(
        species=args.sp2_tag))
    contig_selection(group_dict, args.sp2_cer_scores, args.sp2_tag, "cer")
    print("***** Step 13: Selection of additional {species} marita contigs based on TransRate-scores *****".format(
        species=args.sp2_tag))
    contig_selection(group_dict, args.sp2_mar_scores, args.sp2_tag, "mar")
    print("***** Step 14: {tag} output file creating *****".format(tag=args.tag))
    write_output(group_dict, args.tag, args.sp1_tag, args.sp2_tag)

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
                    help="table with correct groups of homologs and orthologs "
                         "['append_annotation_and_select_correct_groups' output file]")
parser.add_argument('--red_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for rediae transcriptome "
                         "[FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--red_exp_1', type=argparse.FileType('r'), required=True,
                    help="first file (biological replicate) with TPM-values for rediae "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--red_exp_2', type=argparse.FileType('r'), required=True,
                    help="second file (biological replicate) with TPM-values for rediae "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--red_scores', type=argparse.FileType('r'), required=True,
                    help="file with scores for contigs from rediae transcriptome "
                         "[TransRate (Smith-Unna et al., 2016) output file: 'contig.csv']")
parser.add_argument('--cer_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for cercaria transcriptome "
                         "[FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--cer_exp_1', type=argparse.FileType('r'), required=True,
                    help="first file (biological replicate) with TPM-values for cercaria "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--cer_exp_2', type=argparse.FileType('r'), required=True,
                    help="second file (biological replicate) with TPM-values for cercaria "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--cer_scores', type=argparse.FileType('r'), required=True,
                    help="file with scores for contigs from cercaria transcriptome "
                         "[TransRate (Smith-Unna et al., 2016) output file: 'contig.csv']")
parser.add_argument('--mar_anno', type=argparse.FileType('r'), required=True,
                    help="file with annotation for marita transcriptome "
                         "[FunctionAnnotator (Chen et al., 2017) output file: 'AnnotationTable']")
parser.add_argument('--mar_exp_1', type=argparse.FileType('r'), required=True,
                    help="first file (biological replicate) with TPM-values for marita "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--mar_exp_2', type=argparse.FileType('r'), required=True,
                    help="second file (biological replicate) with TPM-values for marita "
                         "[salmon output file: 'quant.sf']")
parser.add_argument('--mar_scores', type=argparse.FileType('r'), required=True,
                    help="file with scores for contigs from marita transcriptome "
                         "[TransRate (Smith-Unna et al., 2016) output file: 'contig.csv']")
parser.add_argument('--species', type=str, required=True,
                    help="species name tag ['Psilo' or 'Sphaer', for instance]")
parser.add_argument('--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def read_table_with_groups(tab, group_dict, species, all_values):
    """
    The function reads table with groups of homologs and orthologs and append them to the python dictionary.
    Last one have next structure:
    key = annotation/species_name,
    value for this key is python dictionary,
    where contigs is a keys and values for each of them is other dictionary with TPM values (equal to zero at this step)
    @param tab: table with groups of similar sequences
    @param group_dict: python dictionary
    @param species: species name tag
    @param all_values: python list
    """
    head = tab.readline()
    for line in tab:
        group_description = line.strip().split("\t")
        group_name, values = group_description[0], group_description[1:]
        contigs_from_species = []
        for contig in values:
            if species in contig:
                contigs_from_species.append(contig)
        if len(contigs_from_species) != 0:
            if len(contigs_from_species) == 1:
                phase = contigs_from_species[0].split("_")[1]
                group_dict["{name}/{species}_{phase}".format(name=group_name, species=species, phase=phase)] \
                    = {contigs_from_species[0]: {"rep_1": 0, "rep_2": 0}}
            else:
                group_dict["{name}/{species}".format(name=group_name, species=species)] \
                    = {contig_ID: {"rep_1": 0, "rep_2": 0} for contig_ID in contigs_from_species}
            all_values.extend(contigs_from_species)


def adding_additional_contigs(group_dict, file_with_annotation, species, phase, all_values):
    """
    The function append additional contigs in previously created python dictionary with groups of homologs.
    As additional(specific) sequences we consider contigs without similar sequences in other transcriptomes
    @param group_dict: previously created python dictionary
    @param file_with_annotation: table with transcriptome annotation
    @param species: species name tag
    @param phase: life cycle phase name tag
    @param all_values: python list with all contigs from 'group_dict'
    """
    head = file_with_annotation.readline()
    contig_without_hits = 0
    for line in file_with_annotation:
        contig_description = line.strip().split("\t")
        contig_id, best_hit = contig_description[0], contig_description[2]
        if contig_id not in all_values:
            if best_hit == "-":
                name = "{species}_{phase}_unknown_hit_{num}".format(species=species,
                                                                    phase=phase, num=contig_without_hits)
                contig_without_hits += 1
                group_dict[name] = {contig_id: {"rep_1": 0, "rep_2": 0}}
                all_values.append(contig_id)
            else:
                name = "{best_hit}/{species}".format(best_hit=best_hit, species=species)
                if name not in group_dict.keys() \
                        and "{name}_{phase}".format(name=name, phase=phase) not in group_dict.keys():
                    new_name = "{name}_{phase}_specific".format(name=name, phase=phase)
                    if new_name not in group_dict.keys():
                        group_dict[new_name] = [contig_id]
                    else:
                        group_dict[new_name].append(contig_id)


def contig_selection(group_dict, file_with_scores, phase, all_values):
    """
    If several contigs have the same annotation, the function selects only one of them by means of TransRate score
    @param group_dict: previously created python dictionary
    @param file_with_scores: table with contigs TransRate scores
    @param phase: life cycle phase name tag
    @param all_values: python list with all contigs from 'group_dict'
    """
    dict_with_repeats = {}
    contigs = []
    name = "{phase}_specific".format(phase=phase)
    for key, value in group_dict.items():
        if key.endswith(name) and len(value) > 1:
            dict_with_repeats[key] = {contig: 0 for contig in value}
            contigs.extend(value)
        elif key.endswith(name) and len(value) == 1:
            for contig in value:
                group_dict[key] = {contig: {"rep_1": 0, "rep_2": 0}}
                all_values.append(contig)

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
                group_dict[key] = {contig_key: {"rep_1": 0, "rep_2": 0}}
                all_values.append(contig_key)


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


def write_output(group_dict, tag, species):
    with open("{tag}.tab".format(tag=tag), 'a') as output_file:
        output_file.write("Annotation\t{species}_red_1\t{species}_red_2\t{species}_cer_1\t"
                          "{species}_cer_2\t{species}_mar_1\t{species}_mar_2\n".format(species=species))
        for group, values in group_dict.items():
            red_value_1, red_value_2, cer_value_1, cer_value_2, mar_value_1, mar_value_2 = 0, 0, 0, 0, 0, 0
            for contig in values:
                if "{species}_red".format(species=species) in contig:
                        red_value_1 += group_dict[group][contig]["rep_1"]
                        red_value_2 += group_dict[group][contig]["rep_2"]
                elif "{species}_cer".format(species=species) in contig:
                        cer_value_1 += group_dict[group][contig]["rep_1"]
                        cer_value_2 += group_dict[group][contig]["rep_2"]
                elif "{species}_mar".format(species=species) in contig:
                        mar_value_1 += group_dict[group][contig]["rep_1"]
                        mar_value_2 += group_dict[group][contig]["rep_2"]
            output_file.write("{ID}\t{red_1}\t{red_2}\t{cer_1}\t{cer_2}\t{mar_1}\t{mar_2}\n".format(ID=group,
                              red_1=red_value_1, red_2=red_value_2, cer_1=cer_value_1, cer_2=cer_value_2,
                              mar_1=mar_value_1, mar_2=mar_value_2))


if __name__ == "__main__":
    group_dict, all_values = {}, []
    print("***** Step 1: The script parses table with ortho|homologs groups *****")
    read_table_with_groups(args.tab, group_dict, args.species, all_values)
    print("***** Step 2: The script check file with annotation for {species} rediae *****".format(species=args.species))
    adding_additional_contigs(group_dict, args.red_anno, args.species, "red", all_values)
    print("***** Step 3: The script check file with annotation for {species} cercaria *****".format(species=args.species))
    adding_additional_contigs(group_dict, args.cer_anno, args.species, "cer", all_values)
    print("***** Step 4: The script check file with annotation for {species} marita *****".format(species=args.species))
    adding_additional_contigs(group_dict, args.mar_anno, args.species, "mar", all_values)
    print("***** Step 5: Selection of rediae contigs based on TransRate-score *****")
    contig_selection(group_dict, args.red_scores, "red", all_values)
    print("***** Step 6: Selection of cercaria contigs based on TransRate-scores *****")
    contig_selection(group_dict, args.cer_scores, "cer", all_values)
    print("***** Step 7: Selection of marita contigs based on TransRate-scores *****")
    contig_selection(group_dict, args.mar_scores, "mar", all_values)
    print("***** Step 8: The script reads first file with TPM-values of {species} rediae *****".format(
        species=args.species))
    write_expression(group_dict, args.red_exp_1, all_values, "1")
    print("***** Step 9: The script reads second file with TPM-values of {species} rediae *****".format(
        species=args.species))
    write_expression(group_dict, args.red_exp_2, all_values, "2")
    print("***** Step 10: The script reads first file with TPM-values of {species} cercariae *****".format(
        species=args.species))
    write_expression(group_dict, args.cer_exp_1, all_values, "1")
    print("***** Step 11: The script reads second file with TPM-values of {species} cercariae *****".format(
        species=args.species))
    write_expression(group_dict, args.cer_exp_2, all_values, "2")
    print("***** Step 12: The script reads first file with TPM-values of {species} marita *****".format(
        species=args.species))
    write_expression(group_dict, args.mar_exp_1, all_values, "1")
    print("***** Step 13: The script reads second file with TPM-values of {species} marita *****".format(
        species=args.species))
    write_expression(group_dict, args.mar_exp_2, all_values, "2")
    print("***** Step 14: {tag} output file creating *****".format(tag=args.tag))
    write_output(group_dict, args.tag, args.species)

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
parser.add_argument('--contig_tab', type=argparse.FileType('r'), required=True,
                    help="table with correct groups of similar sequences "
                         "['append_annotation_and_select_correct_groups [Two Species]' output file]")
parser.add_argument('--TPM_tab', type=argparse.FileType('r'), required=True,
                    help="table with normalized TPM-values for contigs in groups")
parser.add_argument('--module_tab', type=argparse.FileType('r'), required=True,
                    help="table with genes ID and module colors [WGCNA R package]")
parser.add_argument('--sp1_red_KAAS', type=argparse.FileType('r'), required=True,
                    help="table with results of annotation by means of KAAS "
                         "[URL: http://www.genome.jp/kaas-bin/kaas_main] "
                         "for rediae phase of first species [for instance: Psilo_red]")
parser.add_argument('--sp1_cer_KAAS', type=argparse.FileType('r'), required=True,
                    help="table with results of annotation by means of KAAS "
                         "[URL: http://www.genome.jp/kaas-bin/kaas_main] "
                         "for cercariae phase of first species [for instance: Psilo_cer]")
parser.add_argument('--sp1_mar_KAAS', type=argparse.FileType('r'), required=True,
                    help="table with results of annotation by means of KAAS "
                         "[URL: http://www.genome.jp/kaas-bin/kaas_main] "
                         "for marita phase of first species [for instance: Psilo_mar]")
parser.add_argument('--sp2_red_KAAS', type=argparse.FileType('r'), required=True,
                    help="table with results of annotation by means of KAAS "
                         "[URL: http://www.genome.jp/kaas-bin/kaas_main] "
                         "for rediae phase of second species [for instance: Sphaer_red]")
parser.add_argument('--sp2_cer_KAAS', type=argparse.FileType('r'), required=True,
                    help="table with results of annotation by means of KAAS "
                         "[URL: http://www.genome.jp/kaas-bin/kaas_main] "
                         "for cercariae phase of second species [for instance: Sphaer_cer]")
parser.add_argument('--sp2_mar_KAAS', type=argparse.FileType('r'), required=True,
                    help="table with results of annotation by means of KAAS "
                         "[URL: http://www.genome.jp/kaas-bin/kaas_main] "
                         "for marita phase of second species [for instance: Sphaer_mar]")
parser.add_argument('--sp1_red_anno', type=argparse.FileType('r'), required=True,
                    help="table with first species rediae transcriptome annotation "
                         "[FunctionAnnotator (Chen et al., 2017) output file]")
parser.add_argument('--sp1_cer_anno', type=argparse.FileType('r'), required=True,
                    help="table with first species cercariae transcriptome annotation "
                         "[FunctionAnnotator (Chen et al., 2017) output file]")
parser.add_argument('--sp1_mar_anno', type=argparse.FileType('r'), required=True,
                    help="table with first species marita transcriptome annotation "
                         "[FunctionAnnotator (Chen et al., 2017) output file]")
parser.add_argument('--sp2_red_anno', type=argparse.FileType('r'), required=True,
                    help="table with second species rediae transcriptome annotation "
                         "[FunctionAnnotator (Chen et al., 2017) output file]")
parser.add_argument('--sp2_cer_anno', type=argparse.FileType('r'), required=True,
                    help="table with second species cercariae transcriptome annotation "
                         "[FunctionAnnotator (Chen et al., 2017) output file]")
parser.add_argument('--sp2_mar_anno', type=argparse.FileType('r'), required=True,
                    help="table with second species marita transcriptome annotation "
                         "[FunctionAnnotator (Chen et al., 2017) output file]")
parser.add_argument('--species_1_tag', type=str, required=True,
                    help="first species name tag [for instance: Psilo]")
parser.add_argument('--species_2_tag', type=str, required=True,
                    help="second species name tag [for instance: Sphaer]")
parser.add_argument('-t', '--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def read_table_with_colors(module_tab, group_dict, species_1_tag, species_2_tag):
    """
    The function reads table with genes and colors of modules in which they were assigned [WGCNA packages for R] and
    creates for each gene key in new dictionary. As value for this key next information will be stores:
    module color, contigs from each phase of both species, expression level of this gene in TPM (normalized) values,
    result of gene ontology (GO) (biological and molecular level) annotation and results of annotation by means of KAAS
    @param module_tab: table with genes and colors
    @param group_dict: python dictionary
    @param species_1_tag: tag for first species
    @param species_2_tag: tag for second species
    """
    head = module_tab.readline()
    for line in module_tab:
        description = line.strip().split("\t")
        if len(description[0]) > 1:
            annotation, color = description[0][3: -3], description[1]
            group_dict[annotation] = {"color": color,
                                      "{sp_1}_red".format(sp_1=species_1_tag): {"contig": [], "TPM1": 0, "TPM2": 0},
                                      "{sp_1}_cer".format(sp_1=species_1_tag): {"contig": [], "TPM1": 0, "TPM2": 0},
                                      "{sp_1}_mar".format(sp_1=species_1_tag): {"contig": [], "TPM1": 0, "TPM2": 0},
                                      "{sp_2}_red".format(sp_2=species_2_tag): {"contig": [], "TPM1": 0, "TPM2": 0},
                                      "{sp_2}_cer".format(sp_2=species_2_tag): {"contig": [], "TPM1": 0, "TPM2": 0},
                                      "{sp_2}_mar".format(sp_2=species_2_tag): {"contig": [], "TPM1": 0, "TPM2": 0},
                                      "GO_bio": [], "GO_mol": [], "KAAS": []}


def append_contig(group_dict, annotation, species_tag, contig):
    """
    The function append contig in correct list
    @param group_dict: python dictionary
    @param annotation: key for python dictionary (string)
    @param species_tag: tag for species
    @param contig: string
    """
    if "{sp}_red_".format(sp=species_tag) in contig:
        group_dict[annotation]["{sp}_red".format(sp=species_tag)]["contig"].append(contig)
    elif "{sp}_cer_".format(sp=species_tag) in contig:
        group_dict[annotation]["{sp}_cer".format(sp=species_tag)]["contig"].append(contig)
    elif "{sp}_mar_".format(sp=species_tag) in contig:
        group_dict[annotation]["{sp}_mar".format(sp=species_tag)]["contig"].append(contig)


def read_table_with_contigs(contigs_tab, group_dict, species_1_tag, species_2_tag, all_contigs):
    """
    The function reads table with contigs and append each of them in previously created python dictionary
    @param contigs_tab: table with contigs
    @param group_dict: previously created python dictionary
    @param species_1_tag: tag for first species
    @param species_2_tag: tag for second species
    @param all_contigs: python list
    """
    head = contigs_tab.readline()
    for line in contigs_tab:
        description = line.strip().split("\t")
        annotation, contigs = description[0], description[1:]
        if annotation in group_dict.keys():
            for contig in contigs:
                if species_1_tag in contig:
                    append_contig(group_dict, annotation, species_1_tag, contig)
                elif species_2_tag in contig:
                    append_contig(group_dict, annotation, species_2_tag, contig)
            all_contigs.extend([contig for contig in contigs if len(contig) > 1])

            for key in group_dict[annotation].keys():
                if species_1_tag in key or species_2_tag in key:
                    if len(group_dict[annotation][key]["contig"]) == 0:
                        group_dict[annotation][key]["contig"].append("No contig")


def read_table_with_expression(TPM_table, group_dict):
    """
    The function reads table with normalized TPM-values and append them in previously created python dictionary
    @param TPM_table: table with TPMs
    @param group_dict: previously created python dictionary
    """
    header, dict_with_TPMs = [], {}
    header.extend(TPM_table.readline().strip().split("\t"))
    for line in TPM_table:
        description = line.strip().split("\t")
        annotation, TPM_values = description[0], description[1:]
        if annotation in group_dict.keys():
            dict_with_TPMs[annotation] = {sample: float(TPM_values[header.index(sample) - 1]) for sample in header[1:]}
            for sample_key in dict_with_TPMs[annotation].keys():
                if sample_key[:-2] in group_dict[annotation].keys():
                    if sample_key[-1:] == "1":
                        group_dict[annotation][sample_key[:-2]]["TPM1"] = dict_with_TPMs[annotation][sample_key]
                    elif sample_key[-1:] == "2":
                        group_dict[annotation][sample_key[:-2]]["TPM2"] = dict_with_TPMs[annotation][sample_key]


def read_table_with_annotation(anno_table, group_dict, sample_tag, all_contigs):
    """
    The function reads table with results of annotation by means of FunctionAnnotator and append GO results in
    previously created python dictionary
    @param anno_table: table with annotation results
    @param group_dict: previously created python dictionary
    @param sample_tag: tag for sample
    @param all_contigs: python list with all contigs from group_dict dictionary
    """
    head = anno_table.readline()
    bio_counter, mol_counter = 0, 0
    for line in anno_table:
        contig_description = line.strip().split("\t")
        contig_ID, bioGO, molGO = contig_description[0], \
                                  [GO for GO in contig_description[6].split(" | ") if GO != "-"], \
                                  [GO for GO in contig_description[8].split(" | ") if GO != "-"]

        if contig_ID in all_contigs:
            for annotation in group_dict.keys():
                contig_in_dict = group_dict[annotation][sample_tag]["contig"][0]
                if contig_in_dict != "No contig" and contig_ID == contig_in_dict:
                    if len(bioGO) != 0:
                        group_dict[annotation]["GO_bio"].extend(bioGO)
                        bio_counter += 1
                    elif len(bioGO) == 0:
                        group_dict[annotation]["GO_bio"].append("None")

                    if len(molGO) != 0:
                        group_dict[annotation]["GO_mol"].extend(molGO)
                        mol_counter += 1
                    elif len(molGO) == 0:
                        group_dict[annotation]["GO_mol"].append("None")
    print("{sample} : bioGO added: {bio} ; molGO added: {mol} ".format(sample=sample_tag,
                                                                       bio=bio_counter, mol=mol_counter))


def read_table_with_KAAS(table_with_KAAS, group_dict, all_contigs, sample_tag):
    """
    The function reads table with results of KAAS annotation
    @param table_with_KAAS: table with annotation
    @param group_dict: python dictionary
    @param all_contigs: python list with all contigs from group_dict dictionary
    @param sample_tag: tag for sample
    """
    for line in table_with_KAAS:
        description = line.strip().split("\t")
        if len(description) > 1:
            contig_ID, KAAS_ID = description[0][:-3], description[1]
        else:
            contig_ID, KAAS_ID = description[0][:-3], "None"

        if contig_ID in all_contigs:
            for annotation in group_dict.keys():
                contig_in_dict = group_dict[annotation][sample_tag]["contig"][0]
                if contig_in_dict != "No contig" and contig_ID == contig_in_dict:
                    group_dict[annotation]["KAAS"].append(KAAS_ID)


def check_long_list(long_list, frequent_hit):
    """
    The function checks list of values and select only one of them, that has maximal frequency.
    If maximal frequency equal to 1, the function 'return' list of values.
    @param long_list: python list
    @param frequent_hit: python dictionary
    """
    max_frequency = 0
    for hit in set(long_list):
        if long_list.count(hit) > max_frequency and hit != "None":
            max_frequency = long_list.count(hit)
            frequent_hit.clear()
            frequent_hit.append(hit)

    if max_frequency == 1:
        frequent_hit.clear()
        frequent_hit.extend([hit for hit in long_list if hit != "None"])


def check_GO(group_dict, annotation, frequent_GO, frequent_GO_key, GO_cat):
    """
    The function checks results of GO annotation
    @param group_dict: previously created python dictionary
    @param annotation: key for python dictionaries (group_dict and frequent_GO) (string)
    @param frequent_GO: python dictionary
    @param frequent_GO_key: second key for python dictionary 'frequent_GO' (string)
    @param GO_cat: gene ontology category (string)
    """
    link = group_dict[annotation][GO_cat]

    if len(set(link)) == 1:
        frequent_GO[annotation][frequent_GO_key].append(link[0])
    elif len(set(link)) > 1:
        frequent_hit = []
        check_long_list(link, frequent_hit)
        frequent_GO[annotation][frequent_GO_key].append(frequent_hit[0])


def find_frequent_GO(group_dict, frequent_GO):
    """
    The function searches hit of GO annotation with maximal frequency for each gene
    @param group_dict: previously created python dictionary
    @param frequent_GO: python dictionary
    """
    for annotation in group_dict.keys():
        frequent_GO[annotation] = {"frequent_GO_bio": [], "frequent_GO_mol": []}
        check_GO(group_dict, annotation, frequent_GO, "frequent_GO_bio", "GO_bio")
        check_GO(group_dict, annotation, frequent_GO, "frequent_GO_mol", "GO_mol")


def find_frequent_KAAS(group_dict, frequent_KAAS):
    """
    The function searches KAAS hit with maximal frequency for each gene
    @param group_dict: previously created python dictionary
    @param frequent_KAAS: python dictionary
    """
    for annotation in group_dict.keys():
        link = group_dict[annotation]["KAAS"]
        frequent_KAAS[annotation] = []

        if len(set(link)) == 1:
            frequent_KAAS[annotation].append(link[0])
        elif len(set(link)) > 1:
            frequent_hit = []
            check_long_list(link, frequent_hit)
            frequent_KAAS[annotation].append(frequent_hit[0])

        if len(frequent_KAAS[annotation]) == 0:
            # this gene was not annotated by means of KAAS
            frequent_KAAS[annotation].append("No info")


def write_output(group_dict, frequent_GO, frequent_KAAS, tag, species_1_tag, species_2_tag):
    with open("{tag}.tab".format(tag=tag), 'a') as output_file:
        output_file.write("Gene_ID\tModule_color\tFrequent_bio_GO\tFrequent_mol_GO\tFrequent_KAAS\t"
                          "{species_1}_red_contig\t{species_1}_cer_contig\t{species_1}_mar_contig\t"
                          "{species_2}_red_contig\t{species_2}_cer_contig\t{species_2}_mar_contig\t"
                          "{species_1}_red_1_normTPM\t{species_1}_red_2_normTPM\t{species_1}_cer_1_normTPM\t"
                          "{species_1}_cer_2_normTPM\t{species_1}_mar_1_normTPM\t{species_1}_mar_2_normTPM\t"
                          "{species_2}_red_1_normTPM\t{species_2}_red_2_normTPM\t{species_2}_cer_1_normTPM\t"
                          "{species_2}_cer_2_normTPM\t{species_2}_mar_1_normTPM\t{species_2}_mar_2_normTPM\n".format(
                            species_1=species_1_tag, species_2=species_2_tag))
        for annotation in group_dict.keys():
            output_file.write("{anno}\t{color}\t{bioGO}\t{molGO}\t{KAAS}\t{sp1_red}\t{sp1_cer}\t{sp1_mar}\t"
                              "{sp2_red}\t{sp2_cer}\t{sp2_mar}\t{sp1_red_1}\t{sp1_red_2}\t{sp1_cer_1}\t{sp1_cer_2}\t"
                              "{sp1_mar_1}\t{sp1_mar_2}\t{sp2_red_1}\t{sp2_red_2}\t{sp2_cer_1}\t{sp2_cer_2}\t"
                              "{sp2_mar_1}\t{sp2_mar_2}\n".format(
                                anno=annotation, color=group_dict[annotation]["color"],
                                bioGO=frequent_GO[annotation]["frequent_GO_bio"][0],
                                molGO=frequent_GO[annotation]["frequent_GO_mol"][0],
                                KAAS=frequent_KAAS[annotation][0],

                                sp1_red=group_dict[annotation]["{sp1}_red".format(sp1=species_1_tag)]["contig"][0],
                                sp1_cer=group_dict[annotation]["{sp1}_cer".format(sp1=species_1_tag)]["contig"][0],
                                sp1_mar=group_dict[annotation]["{sp1}_mar".format(sp1=species_1_tag)]["contig"][0],
                                sp2_red=group_dict[annotation]["{sp2}_red".format(sp2=species_2_tag)]["contig"][0],
                                sp2_cer=group_dict[annotation]["{sp2}_cer".format(sp2=species_2_tag)]["contig"][0],
                                sp2_mar=group_dict[annotation]["{sp2}_mar".format(sp2=species_2_tag)]["contig"][0],

                                sp1_red_1=group_dict[annotation]["{sp1}_red".format(sp1=species_1_tag)]["TPM1"],
                                sp1_red_2=group_dict[annotation]["{sp1}_red".format(sp1=species_1_tag)]["TPM2"],
                                sp1_cer_1=group_dict[annotation]["{sp1}_cer".format(sp1=species_1_tag)]["TPM1"],
                                sp1_cer_2=group_dict[annotation]["{sp1}_cer".format(sp1=species_1_tag)]["TPM2"],
                                sp1_mar_1=group_dict[annotation]["{sp1}_mar".format(sp1=species_1_tag)]["TPM1"],
                                sp1_mar_2=group_dict[annotation]["{sp1}_mar".format(sp1=species_1_tag)]["TPM2"],
                                sp2_red_1=group_dict[annotation]["{sp2}_red".format(sp2=species_2_tag)]["TPM1"],
                                sp2_red_2=group_dict[annotation]["{sp2}_red".format(sp2=species_2_tag)]["TPM2"],
                                sp2_cer_1=group_dict[annotation]["{sp2}_cer".format(sp2=species_2_tag)]["TPM1"],
                                sp2_cer_2=group_dict[annotation]["{sp2}_cer".format(sp2=species_2_tag)]["TPM2"],
                                sp2_mar_1=group_dict[annotation]["{sp2}_mar".format(sp2=species_2_tag)]["TPM1"],
                                sp2_mar_2=group_dict[annotation]["{sp2}_mar".format(sp2=species_2_tag)]["TPM2"]))

if __name__ == "__main__":
    group_dict, all_contigs, frequent_GO, frequent_KAAS = {}, [], {}, {}
    print("***** Step 1: The script reads table with module colors *****")
    read_table_with_colors(args.module_tab, group_dict, args.species_1_tag, args.species_2_tag)
    print("***** Step 2: The script reads table with contigs *****")
    read_table_with_contigs(args.contig_tab, group_dict, args.species_1_tag, args.species_2_tag, all_contigs)
    print("***** Step 3: The script reads table with expression data *****")
    read_table_with_expression(args.TPM_tab, group_dict)
    print("***** Step 4: The script reads tables with transcriptomes annotation *****")
    read_table_with_annotation(args.sp1_red_anno, group_dict, "{sp1}_red".format(sp1=args.species_1_tag), all_contigs)
    read_table_with_annotation(args.sp1_cer_anno, group_dict, "{sp1}_cer".format(sp1=args.species_1_tag), all_contigs)
    read_table_with_annotation(args.sp1_mar_anno, group_dict, "{sp1}_mar".format(sp1=args.species_1_tag), all_contigs)
    read_table_with_annotation(args.sp2_red_anno, group_dict, "{sp2}_red".format(sp2=args.species_2_tag), all_contigs)
    read_table_with_annotation(args.sp2_cer_anno, group_dict, "{sp2}_cer".format(sp2=args.species_2_tag), all_contigs)
    read_table_with_annotation(args.sp2_mar_anno, group_dict, "{sp2}_mar".format(sp2=args.species_2_tag), all_contigs)
    print("***** Step 5: The script reads tables with KAAS annotation *****")
    read_table_with_KAAS(args.sp1_red_KAAS, group_dict, all_contigs, "{sp1}_red".format(sp1=args.species_1_tag))
    read_table_with_KAAS(args.sp1_cer_KAAS, group_dict, all_contigs, "{sp1}_cer".format(sp1=args.species_1_tag))
    read_table_with_KAAS(args.sp1_mar_KAAS, group_dict, all_contigs, "{sp1}_mar".format(sp1=args.species_1_tag))
    read_table_with_KAAS(args.sp2_red_KAAS, group_dict, all_contigs, "{sp2}_red".format(sp2=args.species_2_tag))
    read_table_with_KAAS(args.sp2_cer_KAAS, group_dict, all_contigs, "{sp2}_cer".format(sp2=args.species_2_tag))
    read_table_with_KAAS(args.sp2_mar_KAAS, group_dict, all_contigs, "{sp2}_mar".format(sp2=args.species_2_tag))
    print("***** Step 6: The script check GO annotation *****")
    find_frequent_GO(group_dict, frequent_GO)
    print("***** Step 7: The script check KAAS annotation *****")
    find_frequent_KAAS(group_dict, frequent_KAAS)
    print("***** Step 8: {tag} output file creating *****".format(tag=args.tag))
    write_output(group_dict, frequent_GO, frequent_KAAS, args.tag, args.species_1_tag, args.species_2_tag)
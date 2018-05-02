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
#parser.add_argument('--tab', type=argparse.FileType('r'), required=True,
#                    help="file with all ortho|homologs groups ['create_groups_of_similar_seq' output file]")
parser.add_argument('--Pc', type=argparse.FileType('r'), required=True,
                    help="file with selected contigs from Psilotrema simillimum cercariae transcriptome "
                         "['selection_of_one_contig_from_phase.py' output file]")
parser.add_argument('--Pr', type=argparse.FileType('r'), required=True,
                    help="file with selected contigs from P.simillum rediae transcriptome "
                         "['selection_of_one_contig_from_phase.py' output file]")
parser.add_argument('--Pm', type=argparse.FileType('r'), required=True,
                    help="file with selected contigs from P.simillum marita transcriptome "
                         "['selection_of_one_contig_from_phase.py' output file]")
parser.add_argument('--Sc', type=argparse.FileType('r'), required=True,
                    help="file with selected contigs from Sphaeridiotrema pseudoglobulus cercariae transcriptome "
                         "['selection_of_one_contig_from_phase.py' output file]")
parser.add_argument('--Sr', type=argparse.FileType('r'), required=True,
                    help="file with selected contigs from S.pseudoglobulus rediae transcriptomes "
                         "['selection_of_one_contig_from_phase.py' output file]")
parser.add_argument('--Sm', type=argparse.FileType('r'), required=True,
                    help="file with selected contigs from S.pseudoglobulus marita transcriptomes "
                         "['selection_of_one_contig_from_phase.py' output file]")
parser.add_argument('--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def dict_without_repeats(file, dict):
    """
    Function append previously selected contigs in new python dictionary
    @param file: file with selected contigs from one life cycle phase
    @param dict: python dictionary
    """
    for line in file:
        group_description = line.strip().split("\t")
        group_name, selected_contig = group_description[0], group_description[1]
        if group_name not in dict.keys():
            dict[group_name] = [selected_contig]
        else:
            dict[group_name].append(selected_contig)


def create_output(new_dict, tag):
    with open("{output_tag}.tab".format(output_tag=tag), 'a') as output:
        output.write('Group_ID\tPsilo_red\tPsilo_cer\tPsilo_mar\tSphaer_red\tSphaer_cer\tSphaer_mar\n')
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
    without_repeats = {}
    print("***** Step 1: The script reads file with contigs from P.simillimum cercariae *****")
    dict_without_repeats(args.Pc, without_repeats)
    print("***** Step 2: The script reads file with contigs from P.simillimum rediae *****")
    dict_without_repeats(args.Pr, without_repeats)
    print("***** Step 3: The script reads file with contigs from P.simillimum marita *****")
    dict_without_repeats(args.Pm, without_repeats)
    print("***** Step 4: The script reads file with contigs from S.pseudoglobulus cercariae *****")
    dict_without_repeats(args.Sc, without_repeats)
    print("***** Step 5: The script reads file with contigs from S.pseudoglobulus rediae *****")
    dict_without_repeats(args.Sr, without_repeats)
    print("***** Step 6: The script reads file with contigs from S.pseudoglobulus marita *****")
    dict_without_repeats(args.Sm, without_repeats)
    print("***** Step 7: {tag} output file creating *****".format(tag=args.tag))
    create_output(without_repeats, args.tag)


# =====================================================================================================================
# TransRate [URL:http://hibberdlab.com/transrate/index.html] 'contig.csv' reminder:
#
# p_seq_true (TransRate v1.0.1) == sCnuc(TransRate v1.0.3) ==
# == "A measure of whether each base has been called correctly"
#
# p_bases_covered (TransRate v1.0.1) = sCcov (TransRate v1.0.3) ==
# == "A measure of whether each base is truly part of the transcript"
#
# p_good (TransRate v1.0.1) == sCord (TransRate v1.0.3) ==
# == "The probability that the contig is derived from a single transcript"
#
# p_not_segmented (TransRate v1.0.1) == sCseg (TransRate v1.0.3) ==
# == "The probability that the contig is structurally complete and correct score"
#
# -- The contig score is the product of the components --
# =====================================================================================================================

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
parser.add_argument('--input', type=argparse.FileType('r'), required=True,
                    help="table with groups of similar sequences ['create_groups_of_similar_seq.py' output file] ")
parser.add_argument('--contigs', type=argparse.FileType('r'), required=True,
                    help="file with metrics for each contigs of one sample "
                         "[TransRate (v1.0.1) (Smith-Unna et al., 2016) output file: 'contig.csv'")
parser.add_argument('--sample_tag', type=str, required=True, help="sample name tag [for instance: 'Psilo_cer']")
parser.add_argument('--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def read_table_with_sequence(tab, dict, sample_tag):
    """
    The function reads table with groups of similar sequences and writes contigs IDs in python dictionary.
    Dictionary have next structure: key = group name;
    value for this key = dictionary with sequences and their scores, which at this step equal to zero
    @param tab: table with groups of similar sequences
    @param dict: python dictionary
    @param sample_tag: tag for sample
    """
    for line in tab:
        group = line.strip().split("\t")
        name, orthologs = group[0], group[1:]
        dict[name] = {contig: 0 for contig in orthologs if sample_tag in contig}


def parsing_contigs_csv(transrate, dict):
    """
    The function parses TransRate output file and
    append contig scores in previously created dictionary with groups of similar sequences
    @param transrate: TransRate output file ('contig.csv')
    @param dict: python dictionary with groups
    """
    desired_contigs = []
    for group_key in dict.keys():
        desired_contigs.extend(dict[group_key].keys())

    transrate_dict = {}
    header_csv = transrate.readline().split(',')
    for line in transrate:
        contig_metrics = line.strip().split(',')
        name, score = contig_metrics[0], float(contig_metrics[8])
        if name in desired_contigs:
            transrate_dict[name] = score

    for contig_key in transrate_dict.keys():
        for group_key, contigs in dict.items():
            for contig in contigs.keys():
                if contig_key == contig:
                    dict[group_key][contig] += transrate_dict[contig_key]


def selection(dict, tag):
    """
    Based on scores best contig per group are selected and this information writes in output file.
    This step help solve the problem when in group we have more than one sequence from one life cycle phase.
    @param dict: python dictionary with groups
    @param tag: tag for output file
    """
    for group_key in dict.keys():
        for contig_key, score in dict[group_key].items():
            if score == max([score for score in dict[group_key].values()]):
                with open("{output_tag}.tab".format(output_tag=tag), "a") as output:
                    output.write("{group_name}\t{best_contig}\n".format(group_name=group_key, best_contig=contig_key))
                break


if __name__ == "__main__":
    contigs_dict = {}
    print("***** Step 1: The script reads table with ortho|homologs groups *****")
    read_table_with_sequence(args.input, contigs_dict, args.sample_tag)
    print("***** Step 2: The script parses TransRate contigs.csv file *****")
    parsing_contigs_csv(args.contigs, contigs_dict)
    print("***** Step 3: The script selects of contigs based on TransRate-score *****")
    selection(contigs_dict, args.tag)

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
    from Bio import SeqIO
except ImportError:
    print("Please check if 'biopython' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True,
                    help="table with groups of similar sequences "
                         "['append_annotation_and_select_correct_groups.py' output file']")
parser.add_argument('--sample_tag', type=str, required=True,
                    help="tag for sample [phase] which aminoacid sequences you want to extract "
                         "[for instance: 'Psilo_cer']")
parser.add_argument('--fasta', type=argparse.FileType('r'), required=True,
                    help="file with nucleotide sequences [fasta file with assembled contigs]")
parser.add_argument('--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def read_table_with_groups(tab, sample_tag, dict_with_similar):
    """
    The function reads table with groups of similar sequences
    and write only contigs from desired phase to new python dictionary.
    @param tab: table with groups of similar sequences
    @param sample_tag: tag for desired sample [phase]
    @param dict_with_similar: python dictionary
    """
    head = tab.readline()
    for line in tab:
        group_description = line.strip().split("\t")
        group_name, values = group_description[0], group_description[1:]
        for contig in values:
            if sample_tag in contig:
                dict_with_similar[group_name] = contig
                break


def fasta_file_parsing(fasta_file, dict_with_similar, tag):
    """
    Function use SeqIO from BioPython for parsing input file, that include nucleotide sequences of assembled contigs.
    In addition function write nucleotide sequences of selected contigs in new fasta-format file
    @param fasta_file: input file with contigs
    @param dict_with_similar: previously created python dictionary with contigs from sample
    @param tag: tag for the created output fasta file
    """
    contigs = SeqIO.parse(fasta_file, 'fasta')
    wanted_contigs = []
    for contig in dict_with_similar.values():
        wanted_contigs.append(contig)
    SeqIO.write((seq for seq in contigs if seq.id in wanted_contigs), "{tag}.fasta".format(tag=tag), "fasta")

if __name__ == "__main__":
    dict_with_similar = {}
    print("***** Step 1: The script reads table with groups of similar sequences *****")
    read_table_with_groups(args.tab, args.sample_tag, dict_with_similar)
    print("***** Step 2: Parsing and creating of {tag}.fasta file *****".format(tag=args.tag))
    fasta_file_parsing(args.fasta, dict_with_similar, args.tag)

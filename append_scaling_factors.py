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
                    help="table with correct groups with TPM-values of similar sequences (before normalization) "
                         "['create_table_with_TPM' output file without normalization]")
parser.add_argument('--scaling_factors', type=argparse.FileType('r'), required=True,
                    help="file with scaling factors "
                         "[normalize.pl (Glusman et al., 2013) (--method net) stdout output file]")
parser.add_argument('-t', '--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def read_table_with_groups(tab, header, group_dict):
    """
    The function reads table with TPM-values of ortho|homologs contigs and write them to new python dictionary
    @param tab: table with TPM-values of similar sequences
    @param header: python list
    @param group_dict: python dictionary
    """
    header.extend(tab.readline().strip().split("\t"))
    for line in tab:
        group_description = line.strip().split("\t")
        annotation, TPM_values = group_description[0], group_description[1:]
        group_dict[annotation] = {sample: float(TPM_values[header.index(sample) - 1]) for sample in header[1:]}


def read_file_with_scaling_factors(file, scaling_factors):
    """
    The function reads file with scaling factors for TPM-values, which will be used for normalization
    @param file: file with scaling factors
    @param scaling_factors: python dictionary
    """
    for line in file:
        description = line.strip().split("\t")
        sample, scaling_factor = description[0], description[1]
        scaling_factors[sample] = float(scaling_factor)


def normalization(group_dict, scaling_factors):
    """
    The function used scaling factors for TPM-values normalization
    @param group_dict: previously created python dictionary with annotation and TPM-values
    @param scaling_factors: python dictionary with factors
    """
    for annotation_key, samples in group_dict.items():
        for sample in samples.keys():
            TPM_value = group_dict[annotation_key][sample]
            if TPM_value != 0:
                group_dict[annotation_key][sample] = round(TPM_value * scaling_factors[sample], 2)


def write_output(group_dict, tag, header):
    with open("{tag}.tab".format(tag=tag), 'a') as output_file:
        output_file.write("{samples}\n".format(samples="\t".join(header)))
        for annotation in group_dict.keys():
            TPM_values = [str(group_dict[annotation][sample_name]) for sample_name in header[1:]]
            output_file.write("{annotation}\t{TPMs}\n".format(annotation=annotation, TPMs="\t".join(TPM_values)))


if __name__ == "__main__":
    group_dict, header, scaling_factors = {}, [], {}
    print("***** Step 1: The script parses table with ortho|homologs groups *****")
    read_table_with_groups(args.tab, header, group_dict)
    print("***** Step 2: The script parses file with scaling factors *****")
    read_file_with_scaling_factors(args.scaling_factors, scaling_factors)
    print("***** Step 3: Normalization *****")
    normalization(group_dict, scaling_factors)
    print("***** Step 4: {tag} output file creating *****".format(tag=args.tag))
    write_output(group_dict, args.tag, header)

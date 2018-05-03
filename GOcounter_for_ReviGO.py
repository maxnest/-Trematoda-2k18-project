"""
Author: Nesterenko Maksim
E-mail: maxnest@ro.ru
Organization: St.Petersburg State University, Russia
Department: Department of Invertebrate zoology
Data: 03.05.2018
"""

try:
    import argparse
except ImportError:
    print("Please check if module 'argparse' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True,
                    help="table with genes, module colors, GOterms (bio and mol): "
                         "genes ID = 1 column, module color = 2 column, "
                         "GO bio and GO mol = 3 and 4 column, respectively"
                         "[for instance: 'genes_in_module' output file]")
parser.add_argument('-t', '--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


def read_table(table, group_dict):
    """
    The function reads table and adds GO terms in python dictionary.
    As a key in it module color is used
    @param table: table with module colors anr gene ontology terms
    @param group_dict: python dictionary
    """
    head = table.readline()
    for line in table:
        description = line.strip().split("\t")
        geneID, color, bioGO, molGO = description[0], description[1], \
                                      description[2].split(" "), description[3].split(" ")
        if color not in group_dict.keys():
            group_dict[color] = {"bio": [bioGO[0]], "mol": [molGO[0]], "bio_count": {}, "mol_count": {}}
        else:
            group_dict[color]["bio"].append(bioGO[0])
            group_dict[color]["mol"].append(molGO[0])


def GO_counter(group_dict, color):
    """
    The function counts the number of each GO term in previously created python dictionary
    @param group_dict: previously created python dictionary
    @param color: module color (key in "group_dict")
    """
    bio_in_module = [GOterm for GOterm in group_dict[color]["bio"] if GOterm != "None"]
    mol_in_module = [GOterm for GOterm in group_dict[color]["mol"] if GOterm != "None"]

    group_dict[color]["bio_count"] = {term: bio_in_module.count(term) for term in set(bio_in_module)}
    group_dict[color]["mol_count"] = {term: mol_in_module.count(term) for term in set(mol_in_module)}


def write_output(group_dict, tag):
    colors = [color for color in group_dict.keys()]
    for color in colors:
        with open("{tag}_{color}_bioGO.tab".format(tag=tag, color=color), 'a') as bio_output_file:
            for bio_term, bio_count in group_dict[color]["bio_count"].items():
                bio_output_file.write("{term}\t{count}\n".format(term=bio_term, count=bio_count))

        with open("{tag}_{color}_molGO.tab".format(tag=tag, color=color), 'a') as mol_output_file:
            for mol_term, mol_count in group_dict[color]["mol_count"].items():
                mol_output_file.write("{term}\t{count}\n".format(term=mol_term, count=mol_count))


if __name__ == "__main__":
    group_dict = {}
    print("***** Step 1: The script reads table with genes, module colors and results of GO annotation *****")
    read_table(args.tab, group_dict)
    print("***** Step 2: The script counts the number of each GO term in each module *****")
    for color_key in group_dict.keys():
        GO_counter(group_dict, color_key)
    print("***** Step 3: {tag} output files creating *****".format(tag=args.tag))
    write_output(group_dict, args.tag)




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

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('--nucl', type=argparse.FileType('r'), required=True)
parser.add_argument('--amino', type=argparse.FileType('r'), required=True)
parser.add_argument('--transrate', type=argparse.FileType('r'), required=True)
# parser.add_argument('--PLEK_model', type=argparse.FileType('r'), required=True)
parser.add_argument('--PLEK_Smansoni', type=argparse.FileType('r'), required=True)
parser.add_argument('--CPC2', type=argparse.FileType('r'), required=True)
parser.add_argument('--GO', type=argparse.FileType('r'), required=True)
parser.add_argument('--Ortho', type=argparse.FileType('r'), required=True)
parser.add_argument('--BLAST', type=argparse.FileType('r'), required=True)
parser.add_argument('--domains', type=argparse.FileType('r'), required=True)
parser.add_argument('--KOBAS', type=argparse.FileType('r'), required=True)
parser.add_argument('--red_exp_1', type=argparse.FileType('r'), required=True)
parser.add_argument('--red_exp_2', type=argparse.FileType('r'), required=True)
parser.add_argument('--cer_exp_1', type=argparse.FileType('r'), required=True)
parser.add_argument('--cer_exp_2', type=argparse.FileType('r'), required=True)
parser.add_argument('--mar_exp_1', type=argparse.FileType('r'), required=True)
parser.add_argument('--mar_exp_2', type=argparse.FileType('r'), required=True)
parser.add_argument('--scaling_factor', type=argparse.FileType('r'), required=True)
parser.add_argument('--threshold', type=int, default=5)     # + expression level in all samples (in TPM)
parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()


contig_dict = {}
orthogroups = {}
scaling_factors = {}


def nucl_parsing(contig_dict, nucl_fasta):
    contigs_nucl = SeqIO.parse(nucl_fasta, 'fasta')
    for seq in contigs_nucl:
        contig_dict[seq.id] = {"ID": [seq.id], "nucl": seq.seq, "amino": "", "score": 0,
                               "domains": [], "PLEK_Smansoni": [], "CPC2": [], "GO_bio": [],
                               "GO_mol": [], "GO_cell": [], "Ortho": [], "KEGG_Smansoni": [],
                               "red_exp_1": 0, "red_exp_2": 0,
                               "cer_exp_1": 0, "cer_exp_2": 0,
                               "mar_exp_1": 0, "mar_exp_2": 0,
                               "red_exp_1_norm": 0, "red_exp_2_norm": 0,
                               "cer_exp_1_norm": 0, "cer_exp_2_norm": 0,
                               "mar_exp_1_norm": 0, "mar_exp_2_norm": 0,
                               "red_mean": 0, "cer_mean": 0, "mar_mean": 0,
                               "pattern": [],
                               "Anno_blast": {"hit": [], "identity": []}}


def amino_parsing(contig_dict, amino_fasta):
    contigs_amino = SeqIO.parse(amino_fasta, 'fasta')
    for seq in contigs_amino:
        #amino_ID = seq.id[:-3]
        description = seq.id.strip().split(" ")
        amino_ID = description[0].split(".")[0]
        contig_dict[amino_ID]["amino"] += seq.seq

    for contig in contig_dict.keys():
        if len(contig_dict[contig]["amino"]) == 0:
            contig_dict[contig]["amino"] += "without ORF"

def transrate_parsing(contig_dict, transrate):

    # ===================================================================================================================
    # p_seq_true == sCnuc(Transrate v1.0.3) == "A measure of whether each base has been called correctly"
    # p_bases_covered = sCcov (Transrate v1.0.3) == "A measure of whether each base is truly part of the transcript"
    # p_good == sCord (Transrate v1.0.3) == "The probability that the contig is derived from a single transcript"
    # p_not_segmented == sCseg (Transrate v1.0.3) == "The probability that the contig is structurally complete and correct"
    # score -- The contig score is the product of the components
    # ===================================================================================================================

    header_csv = transrate.readline()
    for line in transrate:
        metrics = line.strip().split(',')
        name, score = metrics[0], metrics[8]
        if name in contig_dict.keys():
            contig_dict[name]["score"] += float(score)


def plek_parsing(contig_dict, tab, tag):
    for line in tab:
        description = line.strip().split("\t")
        contig_ID, molecule_type = description[-1][1:], description[0]
        contig_dict[contig_ID][tag].append(molecule_type)

    for contig in contig_dict.keys():
        if len(contig_dict[contig][tag]) == 0:
            contig_dict[contig][tag].append("Too Short")


def cpc2_parsing(contig_dict, tab):
    head = tab.readline()
    for line in tab:
        description = line.strip().split("\t")
        contig_ID, molecule_type = description[0], description[-1]
        if molecule_type == "noncoding":
            contig_dict[contig_ID]["CPC2"].append("Non-coding")
        else:
            contig_dict[contig_ID]["CPC2"].append("Coding")

    for contig in contig_dict.keys():
        if len(contig_dict[contig]["CPC2"]) == 0:
            contig_dict[contig]["CPC2"].append("Too Short")


def GO_parsing(contig_dict, GO):
    head = GO.readline()
    for line in GO:
        description = line.strip().split("\t")
        contig_ID, category, ID, rank, anno = description[0].strip().split(".")[0], description[1], \
                                              "GO:{ID}".format(ID=description[2]), description[-1], description[3]
        if rank == "1":
            if category == "MF":
                contig_dict[contig_ID]["GO_mol"].append("{ID}|{Anno}".format(ID=ID, Anno=anno))
            elif category == "CC":
                contig_dict[contig_ID]["GO_cell"].append("{ID}|{Anno}".format(ID=ID, Anno=anno))
            elif category == "BP":
                contig_dict[contig_ID]["GO_bio"].append("{ID}|{Anno}".format(ID=ID, Anno=anno))

    for contig, values in contig_dict.items():
        if len(values["GO_mol"]) == 0:
            contig_dict[contig]["GO_mol"].append("-")

        if len(values["GO_cell"]) == 0:
            contig_dict[contig]["GO_cell"].append("-")

        if len(values["GO_bio"]) == 0:
            contig_dict[contig]["GO_bio"].append("-")


def orthogroups_parsing(contig_dict, orthogroups, Ortho):
    head = Ortho.readline()
    for line in Ortho:
        #description = line.strip().split(",")
        description = line.strip().split("\t")
        group_ID, proteins = description[0], description[1:]
        #orthogroups[group_ID] = [].extend([protein for protein in proteins if protein != ""])
        orthogroups[group_ID] = proteins

    for orthogroup, proteins in orthogroups.items():
        for protein in proteins:
            #if protein[:-3] in contig_dict.keys():
            if protein[:-3] in contig_dict.keys():
                contig_dict[protein[:-3]]["Ortho"].append(orthogroup)

    for contig in contig_dict.keys():
        if len(contig_dict[contig]["Ortho"]) == 0:
            contig_dict[contig]["Ortho"].append("-")


def blast_parsing(contig_dict, BLAST):
    head = BLAST.readline()
    for line in BLAST:
        description = line.strip().split("\t")
        #ID, hit_name, hit_desc, identity = description[0][:-3], description[3], description[4], description[8]
        ID, hit_name = description[0].strip().split(".")[0], description[3]
        if hit_name == "no hits":
            contig_dict[ID]["Anno_blast"]["hit"].append("No hit")
            contig_dict[ID]["Anno_blast"]["identity"].append("No hit")
        else:
            contig_dict[ID]["Anno_blast"]["hit"].append(description[4])
            contig_dict[ID]["Anno_blast"]["identity"].append(description[8])

    for contig in contig_dict.keys():
        if contig_dict[contig]["amino"] == "without ORF":
            contig_dict[contig]["Anno_blast"]["hit"].append("-")
            contig_dict[contig]["Anno_blast"]["identity"].append("-")


def domains_parsing(contig_dict, domains):
    head = domains.readline()
    for line in domains:
        description = line.strip().split("\t")
        ID, domain_name = description[0].strip().split(".")[0], description[1]
        contig_dict[ID]["domains"].append(domain_name)

    for contig in contig_dict.keys():
        if len(contig_dict[contig]["domains"]) == 0:
            contig_dict[contig]["domains"].append("-")


def KEGG_parsing(contig_dict, KEGG):
    for number in range(5):
        head = KEGG.readline()

    for line in KEGG:
        description = line.strip().split("\t")
        if len(description) > 1:
            ID, hit = description[0].strip().split(".")[0], description[1]
            #print(description)
            if hit == "None":
                contig_dict[ID]["KEGG_Smansoni"].append("No hit")
            else:
                KEGG_hit = hit.split("||")
                contig_dict[ID]["KEGG_Smansoni"].append(KEGG_hit[0])
        else:
            break

    for contig in contig_dict.keys():
        if len(contig_dict[contig]["KEGG_Smansoni"]) == 0:
            contig_dict[contig]["KEGG_Smansoni"].append("-")

def expression_parsing(contig_dict, file_with_expression, sample_tag):
    head = file_with_expression.readline()
    for line in file_with_expression:
        contig_description = line.strip().split("\t")
        contig_id, TPM = contig_description[0], contig_description[3]
        contig_dict[contig_id][sample_tag] += float(TPM)


def scaling_factor_parsing(file_with_factors, scaling_factors):
    """
    The function reads file with scaling factors for TPM-values, which will be used for normalization
    @param file: file with scaling factors
    @param scaling_factors: python dictionary
    """
    for line in file_with_factors:
        description = line.strip().split("\t")
        sample, scaling_factor = description[0], description[1]
        phase, rep = sample.split("_")[1], sample.split("_")[2]
        scaling_factors["{phase}_{rep}".format(phase=phase, rep=rep)] = float(scaling_factor)


def normalization(contig_dict, scaling_factors):
    for contig, values in contig_dict.items():
        red_1_norm, red_2_norm, cer_1_norm, cer_2_norm, mar_1_norm, mar_2_norm = \
            numpy.round(values["red_exp_1"] * scaling_factors["red_1"]), \
            numpy.round(values["red_exp_2"] * scaling_factors["red_2"]), \
            numpy.round(values["cer_exp_1"] * scaling_factors["cer_1"]), \
            numpy.round(values["cer_exp_2"] * scaling_factors["cer_2"]), \
            numpy.round(values["mar_exp_1"] * scaling_factors["mar_1"]), \
            numpy.round(values["mar_exp_2"] * scaling_factors["mar_2"])

        contig_dict[contig]["red_mean"] += float(numpy.round(numpy.mean([red_1_norm, red_2_norm])))
        contig_dict[contig]["cer_mean"] += float(numpy.round(numpy.mean([cer_1_norm, cer_2_norm])))
        contig_dict[contig]["mar_mean"] += float(numpy.round(numpy.mean([mar_1_norm, mar_2_norm])))


def pattern(contig_dict):
    """ Glusman et al, 2013:
        We defined genes specific to a sample as those with a positive value for Jongeneel’s specificity measure,
        i.e., the gene was observed in one sample at a level higher than the sum of all other samples combined"""
    for contig, values in contig_dict.items():
        red_exp, cer_exp, mar_exp = values["red_mean"], values["cer_mean"], values["mar_mean"]
        if float(red_exp) > float(cer_exp + mar_exp):
            contig_dict[contig]["pattern"].append("rediae-specific")
        elif float(cer_exp) > float(red_exp + mar_exp):
            contig_dict[contig]["pattern"].append("cercaria-specific")
        elif float(mar_exp) > float(red_exp + cer_exp):
            contig_dict[contig]["pattern"].append("marita-specific")
        else:
            # shared <- common for all & generation-specific & parasites-specific
            contig_dict[contig]["pattern"].append("shared")


def write_table(contig_dict, tag, threshold):
    with open("{output}.tab".format(output=tag), 'a') as output:
        output.write("Contig_ID\tTransRate\tPLEK_Smansoni\tCPC2\t"
                     "Orthogroup_ID\tBLAST-hit\tBLAST-identity\tDomains\tKEGG_Smansoni\tGO_bio\tGO_mol\tGO_cell\t"
                     "Rediae_average_norm_TPM\tCercaria_average_norm_TPM\tMarita_average_norm_TPM\t"
                     "Pattern\tNucl_seq\tAmino_seq\n")
        for contig, values in contig_dict.items():
            # FILTER
            if values["red_mean"] >= threshold or values["cer_mean"] >= threshold or values["mar_mean"] >= threshold:
                output.write("{Id}\t{Score}\t{PLEK_Smansoni}\t{CPC2}\t"
                            "{Ortho}\t{Hit}\t{Identity}\t{Domain}\t{KEGG}\t{GO_bio}\t{GO_mol}\t{GO_cell}\t"
                            "{Red_mean}\t{Cer_mean}\t{Mar_mean}\t{Pattern}\t{Nucl_seq}\t{Amino_seq}\n".format(Id=contig,
                            Score=values["score"], PLEK_Smansoni=values["PLEK_Smansoni"][0],
                            CPC2=values["CPC2"][0], Ortho=values["Ortho"][0], Hit=values["Anno_blast"]["hit"][0],
                            Identity=values["Anno_blast"]["identity"][0], Domain="|".join(values["domains"]),
                            KEGG=values["KEGG_Smansoni"][0], GO_bio=values["GO_bio"][0],
                            GO_mol=values["GO_mol"][0], GO_cell=values["GO_cell"][0], Red_mean=values["red_mean"],
                            Cer_mean=values["cer_mean"], Mar_mean=values["mar_mean"], Pattern=values["pattern"][0],
                            Nucl_seq=values["nucl"], Amino_seq=values["amino"]))

if __name__ == "__main__":
    print("*** Step 1: Files with sequences parsing ***")
    nucl_parsing(contig_dict, args.nucl)
    amino_parsing(contig_dict, args.amino)
    print("*** Step 2: TransRate score adding ***")
    transrate_parsing(contig_dict, args.transrate)
    print("*** Step 3: Classification adding ***")
    plek_parsing(contig_dict, args.PLEK_Smansoni, "PLEK_Smansoni")
    cpc2_parsing(contig_dict, args.CPC2)
    print("*** Step 4: Annotation adding ***")
    GO_parsing(contig_dict, args.GO)
    orthogroups_parsing(contig_dict, orthogroups, args.Ortho)
    blast_parsing(contig_dict, args.BLAST)
    domains_parsing(contig_dict, args.domains)
    KEGG_parsing(contig_dict, args.KOBAS)
    print("*** Step 5: Expression level adding ***")
    expression_parsing(contig_dict, args.red_exp_1, "red_exp_1")
    expression_parsing(contig_dict, args.red_exp_2, "red_exp_2")
    expression_parsing(contig_dict, args.cer_exp_1, "cer_exp_1")
    expression_parsing(contig_dict, args.cer_exp_2, "cer_exp_2")
    expression_parsing(contig_dict, args.mar_exp_1, "mar_exp_1")
    expression_parsing(contig_dict, args.mar_exp_2, "mar_exp_2")
    scaling_factor_parsing(args.scaling_factor, scaling_factors)
    print("*** Step 6: Normalization ***")
    normalization(contig_dict, scaling_factors)
    pattern(contig_dict)
    print("*** Step 7: Output table creating ***")
    write_table(contig_dict, args.output, args.threshold)
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
parser.add_argument('--psilo_cer_mar', type=argparse.FileType('r'), required=True,
                    help="table with homologs between transcriptomes of "
                         "Psilotrema simillimum cercaria and marita "
                         "[Transcriptologs.py (Ambrosino and Chiusano, 2017) output file]")
parser.add_argument('--psilo_cer_red', type=argparse.FileType('r'), required=True,
                    help="table with homologs between transcriptomes of "
                         "P.simillium cercaria and rediae "
                         "[Transcriptologs.py (Ambrosino and Chiusano, 2017) output file]")
parser.add_argument('--psilo_mar_red', type=argparse.FileType('r'), required=True,
                    help="table with homologs between transcriptomes of "
                         "P.simillium marita and rediae "
                         "[Transcriptologs.py (Ambrosino and Chiusano, 2017) output file]")
parser.add_argument('--sphaer_cer_mar', type=argparse.FileType('r'), required=True,
                    help="table with homologs between transcriptomes of "
                         "Sphaeridiotrema pseudoglobulus cercaria and marita "
                         "[Transcriptologs.py (Ambrosino and Chiusano, 2017) output file]")
parser.add_argument('--sphaer_cer_red', type=argparse.FileType('r'), required=True,
                    help="table with homologs between transcriptomes of "
                         "S.pseudoglobulus cercaria and rediae "
                         "[Transcriptologs.py (Ambrosino and Chiusano, 2017) output file]")
parser.add_argument('--sphaer_red_mar', type=argparse.FileType('r'), required=True,
                    help="table with homologs between transcriptomes of "
                         "S.pseudoglobulus rediae and marite "
                         "[Transcriptologs.py (Ambrosino and Chiusano, 2017) output file]")
parser.add_argument('--psilo_sphaer_cer', type=argparse.FileType('r'), required=True,
                    help="table with orthologs between transcriptomes of "
                         "P.simillium and S.pseudoglobulus cercaria "
                         "[Transcriptologs.py (Ambrosino and Chiusano, 2017) output file]")
parser.add_argument('--psilo_sphaer_red', type=argparse.FileType('r'), required=True,
                    help="table with orthologs between transcriptomes of "
                         "P.simillium and S.pseudoglobulus rediae "
                         "[Transcriptologs.py (Ambrosino and Chiusano, 2017) output file]")
parser.add_argument('--psilo_sphaer_mar', type=argparse.FileType('r'), required=True,
                    help="table with orthologs between transcriptomes of "
                         "P.simillium and S.pseudoglobulus marita "
                         "[Transcriptologs.py (Ambrosino and Chiusano, 2017) output file]")
parser.add_argument('--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


class undirected_graph:

    def __init__(self):
        self.number_of_vertex = 0
        self.nodes = {}

    def add_vertex(self, sequence):
        """
        The function adds sequence to the undirected graph
        @param sequence: contig ID
        """
        self.sequence = sequence
        self.nodes[sequence] = []
        self.number_of_vertex += 1

    def add_edge(self, seq, other_seq):
        """
        The function create link between two vertex in undirected graph
        @param seq: first sequence (contig ID)
        @param other_seq: second sequence (contig ID)
        """
        self.add_vertex(seq) if seq not in self.nodes.keys() else None
        self.add_vertex(other_seq) if other_seq not in self.nodes.keys() else None

        if seq not in self.nodes[other_seq]:
            self.nodes[other_seq].append(seq)

        if other_seq not in self.nodes[seq]:
            self.nodes[seq].append(other_seq)

    def num_of_vertex(self):
        return self.number_of_vertex

    def vertices(self):
        return self.nodes.keys()

    def neighbourhood(self, seq):
        """
        The function return all vertices which connected with requested sequence (contig ID)
        @param seq: requested sequence (contig ID)
        """
        if seq in self.nodes.keys():
            return self.nodes[seq]


def dfs(vertex, graph, visited_vertex, web):
    """
    dfs = depth-first search;
    The function writes all connected vertices in one set ("web")
    @param vertex: one of the graph vertex
    @param graph: undirected graph
    @param visited_vertex: list of previously visited vertex
    @param web: set of the connected vertices
    """
    web.add(vertex)
    visited_vertex.add(vertex)
    for neighbourhood in graph.neighbourhood(vertex):
        if neighbourhood not in visited_vertex:
            dfs(neighbourhood, graph, visited_vertex, web)


def dfs_plus(graph, output):
    """
    The function used dfs to find all connected component and writes them as group in output file
    @param graph: undirected graph
    @output: output file
    """
    count = 0
    visited_vertex = set()

    while len(visited_vertex) != graph.num_of_vertex():
        web = set()

        for vertex in graph.vertices():
            if vertex not in visited_vertex:
                spider = vertex
                break

        web.add(spider)
        visited_vertex.add(spider)
        dfs(spider, graph, visited_vertex, web)
        count += 1

        with open(output, "a") as output_file:
            output_file.write("group_{num}\t{value}\n".format(num=count, value="\t".join(web)))


def Transcriptologs_output_parser(input_file, orthologs_dict):
    """
    The function parse input file, created by Transcriptologs.py, and append IDs of similar contigs to undirected graph
    @param input_file: input file with pairs of similar sequences
    @param orthologs_dict: python dictionary (class "undirected_graph")
    """
    head = input_file.readline()
    for line in input_file:
        pair_description = line.strip().split("\t")
        seq_1, seq_2 = pair_description[0], pair_description[1]
        orthologs_dict.add_edge(seq_1, seq_2)


if __name__ == "__main__":
    orthologs_dict = undirected_graph()
    print("***** Step 1: The script reads 'Psilo_cer_mar' table *****")
    Transcriptologs_output_parser(args.psilo_cer_mar, orthologs_dict)
    print("***** Step 2: The script reads 'Psilo_cer_red' table *****")
    Transcriptologs_output_parser(args.psilo_cer_red, orthologs_dict)
    print("***** Step 3: The script reads 'Psilo_mar_red' table *****")
    Transcriptologs_output_parser(args.psilo_mar_red, orthologs_dict)
    print("***** Step 4: The script reads 'Sphaer_cer_mar' table *****")
    Transcriptologs_output_parser(args.sphaer_cer_mar, orthologs_dict)
    print("***** Step 5: The script reads 'Sphaer_cer_red' table *****")
    Transcriptologs_output_parser(args.sphaer_cer_red, orthologs_dict)
    print("***** Step 6: The script reads 'Sphaer_red_mar' table *****")
    Transcriptologs_output_parser(args.sphaer_red_mar, orthologs_dict)
    print("***** Step 7: The script reads 'Psilo_Sphaer_cer' table *****")
    Transcriptologs_output_parser(args.psilo_sphaer_cer, orthologs_dict)
    print("***** Step 8: The script reads 'Psilo_Sphaer_red' table *****")
    Transcriptologs_output_parser(args.psilo_sphaer_red, orthologs_dict)
    print("***** Step 9: The script reads 'Psilo_Sphaer_mar' table *****")
    Transcriptologs_output_parser(args.psilo_sphaer_mar, orthologs_dict)
    print("***** Step 10: The script used depth-first searching for group reconstruction "
          "and writes this groups in {tag} output file".format(tag=args.tag))
    dfs_plus(orthologs_dict, "{tag}.tab".format(tag=args.tag))

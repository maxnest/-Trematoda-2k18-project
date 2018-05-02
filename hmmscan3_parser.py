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
    from Bio import SearchIO
except ImportError:
    print("Please check if module 'biopython' is installed")
    quit()

parser = argparse.ArgumentParser()
parser.add_argument('--tab', type=argparse.FileType('r'), required=True,
                    help="table with results of domains searching "
                         "[hmmer/hmmscan output file (--domtblout)]")
parser.add_argument('--tag', type=str, required=True, help="tag for output file")
args = parser.parse_args()


class directed_graph:

    def __init__(self):
        self.number_of_vertex = 0
        self.number_of_edges = 0
        self.nodes = {}

    def add_vertex(self, sequence, score, start, end):
        """
        The function adds domain and its properties (coordinates and hit-score) to the directed graph
        @param sequence: domain ID
        """
        self.sequence = sequence
        self.nodes[sequence] = [{"neighbourhood": []}, {"score": score}, {"start": start}, {"end": end}]
        self.number_of_vertex += 1

    def add_edge(self, seq, other_seq):
        """
        The function creates link between two vertex in directed graph
        @param seq: first domain ID
        @param other_seq: second domain ID
        """
        if other_seq not in self.nodes[seq][0]["neighbourhood"]:
            self.nodes[seq][0]["neighbourhood"].append(other_seq)
        self.number_of_edges += 1


def add_vertex_to_graph(dict_with_contigs, contig_id, potential_protein_graph):
    """
    The function reads table with domains and information in python dictionary (class "directed_graph")
    @param dict_with_contigs: python dictionary with contigs
    @param contig_id: ID of contig
    @param potential_protein_graph: python dictionary of a "directed_graph" class
    """
    for domain in dict_with_contigs[contig_id]:
        for key, value in domain.items():
            domain_name, start, end, score = key, value[-2][0], value[-2][1], value[-1]
            potential_protein_graph.add_vertex(domain_name, score, start, end)


def add_edges_between_domains(potential_protein_graph, vertex_pairs, best_path_for_contig):
    """
    If two domains doesn't overlapped,then function creates link between them.
    In contrast, if all domains in potential protein are overlapped, function selects domain with best score
    @param potential_protein_graph: python dictionary of a "directed_graph" class
    @param vertex_pairs: python list
    @param best_path_for_contig: python list ('best path' in this step means domain with best score)
    """
    for vertex in potential_protein_graph.nodes:
        for other_vertex in potential_protein_graph.nodes:
            if potential_protein_graph.nodes[vertex][3]["end"] \
                    < potential_protein_graph.nodes[other_vertex][2]["start"]:
                potential_protein_graph.add_edge(vertex, other_vertex)
                vertex_pairs.append([vertex, other_vertex])
    if len(vertex_pairs) == 0:  # if all domains in potential protein are overlapped
        the_best_score = 0
        for vertex in potential_protein.nodes:
            domain_score = potential_protein.nodes[vertex][1]["score"]
            if domain_score >= the_best_score:
                the_best_score = domain_score
                best_path_for_contig.clear()
                best_path_for_contig.append(vertex)


def dfs_paths(graph, start, goal):
    """
    dfs = depth first search
    with changes from http://eddmann.com/posts/depth-first-search-and-breadth-first-search-in-python/
    the function use dfs for reconstruction of path between start and goal
    @param graph: python dictionary of "directed graph" class
    @param start: one of the domains ("first" domain)
    @param goal: another domain ("final" domain)
    """
    stack = [(start, [start])]
    while stack:
        (vertex, path) = stack.pop()
        for next in set(graph.nodes[vertex][0]['neighbourhood']) - set(path):
            if next == goal:
                yield path + [next]
            else:
                stack.append((next, path + [next]))


def find_all_possible_paths(pairs, all_possible_paths):
    """
    At this step 'possible path' means combination of domains.
    The function founds all possible paths between pairs of domains
    @param pairs: python list with pairs of connected domains (vertices in directed graph)
    @param all_possible_paths: python list
    """
    for pair in pairs:
        all_possible_paths.extend(list(dfs_paths(potential_protein, pair[0], pair[1])))


def find_score_for_paths(all_possible_paths, paths_with_scores):
    """
    The function summarizes domains scores in all their potential combinations
    @param all_possible_paths: python list with combinations of domains
    @param paths_with_scores: python dictionary
    """
    for path in all_possible_paths:
        path_score = 0
        for domain in path:
            path_score += potential_protein.nodes[domain][1]["score"]
        paths_with_scores[path_score] = path


def select_best_path(paths_with_scores, best_path):
    """
    The function select only one ("best") path between domains
    @param paths_with_scores: previously created python dictionary with scores of paths
    @param best_path: python list
    """
    best_score = max(paths_with_scores.keys())
    for score, path in paths_with_scores.items():
        if score == best_score:
            for domain in path:
                best_path.append(domain)


def write_output(tag, contig, potential_protein, best_path_for_contig):
    domains_with_coord = []
    for domain in best_path_for_contig:
        start, end = potential_protein.nodes[domain][2]["start"], potential_protein.nodes[domain][3]["end"]
        domains_with_coord.append("{domain}: from {start} to {end}".format(domain=domain, start=start, end=end))
    with open("{tag}.tab".format(tag=tag), 'a') as output:
        output.write("{contig_ID}\t{num_of_domains}\t{all_domains}\t{coordinates}\n".format(
                contig_ID=contig, num_of_domains=len(best_path_for_contig),
                all_domains=" ".join(domain for domain in best_path_for_contig),
                coordinates=" ".join(domain for domain in domains_with_coord)
            ))

if __name__ == "__main__":
    dict_with_contigs = {}
    print("***** Step 1: The script parses the hmmscan output file *****")
    hmmscan3_parser = SearchIO.parse(args.tab, 'hmmscan3-domtab')
    for qresult in hmmscan3_parser:
        for hit in qresult.hits:
            for hsp in hit.hsps:
                contig_id, domain_id, domain_desc, coord, score = qresult.id, hit.id, hit.description, \
                                                                  hsp.query_range, hsp.bitscore
                if contig_id not in dict_with_contigs.keys():
                    dict_with_contigs[contig_id] = [{domain_id: [domain_desc, coord, score]}]
                else:
                    dict_with_contigs[contig_id].append({domain_id: [domain_desc, coord, score]})

    print("***** Step 2: The script checks domain architecture of potential protein "
          "and writes result in {tag} output file *****".format(tag=args.tag))

    with open("{tag}.tab".format(tag=args.tag), 'a') as output:
        output.write("Contig_ID\tNum_of_domains\tName_of_all_domains\tCoordinates\n")

    for contig in dict_with_contigs.keys():

        potential_protein, pairs_of_domains, all_possible_paths, \
        paths_with_scores, best_path_for_contig = directed_graph(), [], [], {}, []

        add_vertex_to_graph(dict_with_contigs, contig, potential_protein)
        if len(potential_protein.nodes) >= 2:
            add_edges_between_domains(potential_protein, pairs_of_domains, best_path_for_contig)
            if len(best_path_for_contig) != 0:
                best_domain = best_path_for_contig[0]
            else:
                find_all_possible_paths(pairs_of_domains, all_possible_paths)
                find_score_for_paths(all_possible_paths, paths_with_scores)
                select_best_path(paths_with_scores, best_path_for_contig)
        else:
                for key in potential_protein.nodes.keys():
                    best_path_for_contig.append(key)
        write_output(args.tag, contig, potential_protein, best_path_for_contig)



import logging
import time
import copy
import os
import math
import itertools
import networkx as nx
import json
import csv
from sets import Set
from scipy.stats import hypergeom
from sklearn.cluster import DBSCAN
from utilities import Utilities

MIN_NUM_OF_GENOMES_IN_CLUSTER = 4
WINDOW_SLIDE_SIZE = 16
MIN_NUM_OF_GENES_IN_INTERVAL = 2
MIN_RANKING_SCORE = 20
SIZE_OF_GENOME = 0
NUMBER_OF_GENOMES = 0
GENES_GAP = 2000
GENES_OVERLAP = 100


# The pipeline: CreateGraph: create a bipartite graph G=(A+B,E) where A={gi | gi gene at the centroid genome) A is order by the their gi.start}. B = { geneInterval(i)| each geneInterval(i) is one sliding window of d genes
# from the reference genome where its id is been represented by the color of geneInterval(i) }
def compute_length(root_dir):
    results = Utilities.files_in_dir(root_dir)
    return Utilities.calculate_lengths(root_dir, results)


def calculate_subsets(array_length, edges_group):
    total_subsets = []
    if WINDOW_SLIDE_SIZE < array_length:
        for i in range(WINDOW_SLIDE_SIZE, array_length + 1):
            tmp_arr = []
            for edge in edges_group:
                if i - WINDOW_SLIDE_SIZE <= edge < i:
                    tmp_arr.append(edge)
            tmp_subsets = list(itertools.chain(
                *(list(itertools.combinations(tmp_arr, n)) for n in
                  range(MIN_NUM_OF_GENES_IN_INTERVAL, len(tmp_arr) + 1))))
            total_subsets.append(tmp_subsets)
    else:
        tmp_arr = []
        for edge in edges_group:
            tmp_arr.append(edge)
        tmp_subsets = list(itertools.chain(
            *(list(itertools.combinations(tmp_arr, n)) for n in
              range(MIN_NUM_OF_GENES_IN_INTERVAL, len(tmp_arr) + 1))))
        total_subsets.append(tmp_subsets)
    logging.info("Finish subset")
    return total_subsets


# Calculate ranking score according to the following equation ((n-k)C(d-k)/(nCd)*(n-d))^c
# n: size of genome
# d: length of the interval/operon in the query
# k: number of genes in A (in the operon)
# c: number of genomes in the bicluster
def calculate_log_ranking_score(m, n, d, c, k):
    if d > n:
        d = n
    up = Utilities.ncr(n - k, d - k)
    down = Utilities.ncr(n, d)
    p = float(float(up) / float(down))
    p = p * (n - d)
    i_prime = math.floor(float((m + 1) * float(p)))
    if i_prime > m:
        i_prime = m
    elif i_prime < c:
        i_prime = c
    if p > 1:
        return -1
    ranking_score = -float(
        float(math.log10(m - c + 1)) + float(math.log10(Utilities.ncr(m, i_prime))) + float(i_prime * math.log10(p)) + float(
            (m - i_prime) * math.log10(1 - p)))
    # logging.info( "m=" + str(m) + " n=" + str(n) + " d=" + str(d) + " c=" + str(c) +  " k=" + str(k) + " m=" + str(m) + " up= " + str(up) + "  down= " + str(down) + " p=" + str(p) + " i_prime=" + str(i_prime) + " ranking_score=" + str(ranking_score)
    return ranking_score


# Implementation for edit distance which takes into consideration the strand of the genes. Deletion and mismatch are equal 1 and match with different strands equals 0.5.
# In this function we measure the edit distance of only the substring of genes from the query_set.
def edit_distance(s, t, direction_s, direction_t, query_set):
    s_prime = []
    t_prime = []
    directions_s_prime = []
    directions_t_prime = []
    i = 0
    for gene in s:
        if gene in query_set:
            s_prime.append(gene)
            directions_s_prime.append(direction_s[i])
        i += 1
    i = 0
    for gene in t:
        if gene in query_set:
            t_prime.append(gene)
            directions_t_prime.append(direction_t[i])
        i += 1
    s_prime[1:len(s_prime) + 1] = s_prime[0:len(s_prime)]
    s_prime[0] = ''
    directions_s_prime[1:len(directions_s_prime) + 1] = directions_s_prime[0:len(directions_s_prime)]
    directions_s_prime[0] = ''

    t_prime[1:len(t_prime) + 1] = t_prime[0:len(t_prime)]
    t_prime[0] = ''
    directions_t_prime[1:len(directions_t_prime) + 1] = directions_t_prime[0:len(directions_t_prime)]
    directions_t_prime[0] = ''
    d = {}
    s = len(s_prime)
    t = len(t_prime)
    for i in range(s):
        d[i, 0] = i
    for j in range(t):
        d[0, j] = j
    for j in range(1, t):
        for i in range(1, s):
            if s_prime[i] == t_prime[j]:
                if directions_s_prime[i] == directions_t_prime[j]:
                    d[i, j] = d[i - 1, j - 1]
                else:
                    d[i, j] = min(d[i - 1, j - 1] + 0.5, min(d[i - 1, j], d[i, j - 1]) + 1)
            else:
                d[i, j] = min(d[i - 1, j], d[i, j - 1], d[i - 1, j - 1]) + 1
    return d[s - 1, t - 1]


# Edit distance with all the genes, not only the query
def edit_distance2(s, t, direction_s, direction_t):
    s_prime = []
    t_prime = []
    directions_s_prime = []
    directions_t_prime = []
    i = 0
    for gene in s:
        s_prime.append(gene)
        directions_s_prime.append(direction_s[i])
        i += 1
    i = 0
    for gene in t:
        t_prime.append(gene)
        directions_t_prime.append(direction_t[i])
        i += 1
    s_prime[1:len(s_prime) + 1] = s_prime[0:len(s_prime)]
    s_prime[0] = ''
    directions_s_prime[1:len(directions_s_prime) + 1] = directions_s_prime[0:len(directions_s_prime)]
    directions_s_prime[0] = ''

    t_prime[1:len(t_prime) + 1] = t_prime[0:len(t_prime)]
    t_prime[0] = ''
    directions_t_prime[1:len(directions_t_prime) + 1] = directions_t_prime[0:len(directions_t_prime)]
    directions_t_prime[0] = ''
    d = {}
    s = len(s_prime)
    t = len(t_prime)
    for i in range(s):
        d[i, 0] = i
    for j in range(t):
        d[0, j] = j
    for j in range(1, t):
        for i in range(1, s):
            if s_prime[i] == t_prime[j]:
                if directions_s_prime[i] == directions_t_prime[j]:
                    d[i, j] = d[i - 1, j - 1]
                else:
                    d[i, j] = min(d[i - 1, j - 1] + 0.5, min(d[i - 1, j], d[i, j - 1]) + 1)
            else:
                d[i, j] = min(d[i - 1, j], d[i, j - 1], d[i - 1, j - 1]) + 1
    return d[s - 1, t - 1]


# Gene class,
class Gene(object):
    def __init__(self, id, name, gene_id, start, stop, strand, attribute, target_gene_name, target_gene_attribute, eval):
        self.id = id
        self.name = name
        self.gene_id = gene_id
        self.start = start
        self.stop = stop
        self.strand = strand
        self.attribute = attribute
        self.target_gene_name = target_gene_name
        self.target_gene_attribute = target_gene_attribute
        self.eval = eval

    def to_json(self):
        ans = {'id': self.id, 'name': self.name, 'gene_id': self.gene_id, 'start': self.start, 'stop': self.stop,
               'strand': self.strand,
               'attribute': self.attribute, 'target_gene_name': self.target_gene_name,
               'target_gene_attribute': self.target_gene_attribute, 'eval': self.eval}
        return ans


'''
GeneInterval class, Vertex in group B in the bipartite graph. Each vertex represents a group of genes which are close to each other from one genome.
Each genome is been represented by a color.
'''
class GeneInterval(object):
    def __init__(self, strain_file, color, interval_id):
        self.genes = []
        self.color = color
        self.interval_id = interval_id
        self.numOfGenes = 0
        self.organism = strain_file

    def add_gene(self, gene):
        self.genes.append(gene)
        self.genes.sort(key=lambda x: x.start, reverse=False)
        self.numOfGenes += 1

    def get_gene(self, i):
        return self.genes[i]

    def get_genes(self):
        return self.genes

    def tostring(self):
        ans = ""
        ans = ans + str(self.organism) + "(" + str(self.genes[0].start) + "|" + str(self.genes[-1].stop) + ") "
        for gene in self.genes:
            ans = ans + " " + gene.name + "[" + str(gene.start) + "-" + str(gene.stop) + "|" + str(gene.strand) + "]"
        return ans

    def to_json(self, organism, index, label):
        ans = {'color': self.color, 'id': self.interval_id, 'numOfGenes': self.numOfGenes, 'strain': self.organism,
               'organism': organism, 'genes': [], 'index': index, 'label': label}
        for gene in self.genes:
            ans['genes'].append(gene.to_json())
        return ans


# This class represents a full biclusters where A holds the genes name, B holds the geneIntervals.
class FullBicluster(object):
    def __init__(self):
        self.labels = []
        self.A = []
        self.B = []
        self.ranking_score = 0

    def set_ranking_scoreue(self, ranking_score):
        self.ranking_score = ranking_score

    def get_ranking_scoreue(self):
        return self.ranking_score

    def set_a(self, a):
        self.A = a

    def set_b(self, b):
        self.B = b

    def get_b(self):
        return self.B

    def get_a(self):
        return self.A

    def tostring(self):
        return "A=" + str(self.A) + " B=" + str(sorted(self.B))

    def to_json(self, bip_graph, org_to_index, hyper_geo_array):
        ans = {'hyper_geo_array': hyper_geo_array, 'A': [], 'B': [], 'ranking_scoreue': self.ranking_score}
        for a in self.A:
            ans['A'].append(bip_graph.A[a].to_json())
        i = 0
        for b in self.B:
            organism = bip_graph.strain_to_organism[bip_graph.B[b].organism]
            index = org_to_index[organism]
            ans['B'].append(bip_graph.B[b].to_json(organism, index, self.labels[i]))
            i += 1

        return ans


# BiPartiteGrpah class stores all the information about the bipartite graph
class BiPartiteGrpah(object):
    def __init__(self):
        self.bic_results = []
        self.ranking_score_results = []
        self.A = []
        self.B = []
        self.E = {}
        self.d = 0
        self.n = 0

    def set_d(self, d):
        self.d = d

    def get_d(self):
        return self.d

    def add_gene(self, gene):
        self.A.append(gene)

    def add_gene_interval(self, gene_interval):
        self.B.append(gene_interval)

    def get_edges(self):
        return self.E

    def get_a(self):
        return self.A

    def get_b(self):
        return self.B

    def set_n(self, n):
        self.n = n

    def get_n(self):
        return self.n

    '''
    Add new edge first check if Ai (node from A) is exist than verify that the color
     exist and than check if Bj is not already exist in Ai in color color
    '''

    def add_edge(self, ai, bj):
        if self.E.get(ai) is None:
            self.E[ai] = []
            (self.E[ai]).append(bj)
        else:
            if bj not in self.E[ai]:
                (self.E[ai]).append(bj)

    def count_colors(self, b_prime):
        if len(b_prime) > 0:
            b_prime = sorted(b_prime, key=lambda x: self.B[x].color, reverse=True)
            last_color = self.B[b_prime[0]].color
            ans = 1
            for b in b_prime[1:]:
                if self.B[b].color != last_color:
                    last_color = self.B[b].color
                    ans += 1
            return ans
        return -1

    '''
    Algorithm for computing all d-valid bicliques.
    1. For each node in A
    we compute its Na.
    2. Compute with dynamic programing K: K[i,l,s], all the d-valid bicliques which start in node i in A
       finish at i+l and are at the size of s.
    3. For each group we create an bitarrray size of |B| where each bit represents a node in B.
    '''
    def calculate_biclusters(self):
        a_groups = []
        logging.info("Calculate Bicluster() start")
        total_subsets = calculate_subsets(len(self.A) + 1, self.E)
        for indexes_subsets in total_subsets:
            round = 1
            # logging.info( "Round " + str(round) + " Number of subsets: " + str(len(indexes_subsets))
            round += 1
            for subset in indexes_subsets:
                # Check if the subset is a biclique
                if subset != () and subset[0] in self.get_edges():
                    b_prime = copy.deepcopy(self.get_edges()[subset[0]])
                    for i in subset:
                        if i in self.get_edges():
                            if len(b_prime) > 0:
                                b_prime = list(set(b_prime) & set(self.get_edges()[i]))
                        else:
                            b_prime = []
                    number_of_colors_in_bicluster = 0
                    if len(b_prime) > 0:
                        number_of_colors_in_bicluster = self.count_colors(b_prime)
                    if number_of_colors_in_bicluster > MIN_NUM_OF_GENOMES_IN_CLUSTER:
                        biclique = FullBicluster()
                        biclique.set_a(subset)
                        biclique.set_b(b_prime)
                        ranking_score = calculate_log_ranking_score(NUMBER_OF_GENOMES, SIZE_OF_GENOME, WINDOW_SLIDE_SIZE, number_of_colors_in_bicluster,len(subset))
                        biclique.set_ranking_scoreue(ranking_score)
                        if self.check_bic(biclique) == 1 and ranking_score > MIN_RANKING_SCORE:
                            self.bic_results.append(biclique)
                            a_groups.append(biclique.get_a())

        logging.info("Finish to calculate bicliques")

    '''
    Check for each biclique in the results array:
    If it contains the new biclique (smaller or equal): return 0 else return 1
    If the new biclique contains the biclique delete the biclique from the results array
    '''
    def check_bic(self, bic):
        bic_a = set(bic.get_a())
        delete_array = []
        for bic2 in self.bic_results:
            tmp_a = set(bic2.get_a())

            # the new biclique is contained in some other biclique in the results array
            if bic_a.issubset(tmp_a):
                if bic.get_b() == bic2.get_b():
                    return 0
            if tmp_a.issubset(bic_a):
                if bic.get_b() == bic2.get_b():
                    delete_array.append(bic2)
        for bic3 in delete_array:
            self.bic_results.remove(bic3)
        return 1

    def calculate_distance_matrices(self):
        for bic in self.bic_results:
            instances = []
            taxonomies = []
            direction_instances = []
            for b in bic.B:
                interval = []
                direction_interval = []
                # logging.info( self.B[b]
                taxonomies.append(self.B[b].taxonomy)
                for gene in self.B[b].genes:
                    interval.append(gene.name)
                    direction_interval.append(gene.strand)
                direction_instances.append(direction_interval)
                instances.append(interval)
            bic.distanceMatrix = []
            query_genes = []
            for a in bic.A:
                query_genes.append(self.A[a].name)
            i = 0
            for instance in instances:
                bic.distanceMatrix.append([])
                j = 0
                for instance2 in instances:
                    tmp = float(edit_distance2(copy.deepcopy(instance), copy.deepcopy(instance2), copy.deepcopy(direction_instances[i]),copy.deepcopy(direction_instances[j])))
                    tmp_direction = direction_instances[i]
                    inverse_instance = list(reversed(instance))
                    inverse_direction_instance = list(reversed(tmp_direction))
                    inverse_direction_instance[:] = [int(x) * -1 for x in inverse_direction_instance]
                    inverse_direction_instance[:] = [str(x) for x in inverse_direction_instance]
                    tmp2 = float(edit_distance2(copy.deepcopy(inverse_instance), copy.deepcopy(instance2),
                                     copy.deepcopy(inverse_direction_instance), copy.deepcopy(direction_instances[j])))
                    bic.distanceMatrix[i].append(float(min(tmp, tmp2)))
                    j += 1
                i += 1
            db = DBSCAN(metric='precomputed', eps=0.01, min_samples=1).fit(bic.distanceMatrix)
            bic.labels = db.labels_
            tmp_taxa_set = set(taxonomies)
            bic.taxa_to_index_dict = {}
            bic.index_to_taxa_Dict = {}

            i = 0
            for taxa in tmp_taxa_set:
                bic.taxa_to_index_dict[taxa] = i
                bic.index_to_taxa_Dict[i] = taxa
                i += 1

            bic.cluster_to_taxa = [[bic.labels[i], taxonomies[i]] for i in range(len(bic.B))]
            for tup in bic.cluster_to_taxa:
                tup[1] = bic.taxa_to_index_dict[tup[1]]
            bic.hyperDist = calc_hypergeom(bic.cluster_to_taxa, len(set(bic.labels)), len(bic.taxa_to_index_dict))


# This is not relevant for the 1.0.0 version.
def calc_hypergeom(cluster_taxa_tuples_array, num_of_clusters, num_of_taxa):
    """
    The method re
    Parameters:
      cluster_taxa_tuples_array - array of tuples where each tuple holds for (cluster id ,taxa)
      num_of_clusters - number of clusters

    Output:
        stats: array of num_of_clusters objects {best_dist,best_taxa}}. For each cluster there is an object that holds the taxa that got the best dist and the dist itself.
    """
    m = len(cluster_taxa_tuples_array)
    n_big_array = [0] * num_of_clusters
    n_array = [0] * num_of_taxa
    x_matrix = []
    tmp = [0 for i in range(num_of_taxa)]
    for j in range(num_of_clusters):
        x_matrix.append(copy.deepcopy(tmp))

    for tup in cluster_taxa_tuples_array:
        n_big_array[tup[0]] += 1
        n_array[tup[1]] += 1
        x_matrix[tup[0]][tup[1]] += 1

    stats = [{} for j in range(num_of_clusters)]
    for i in range(num_of_clusters):
        best_dist = 1.1
        best_taxa = -1
        N = n_big_array[i]
        for j in range(num_of_taxa):
            n = n_array[j]
            x = x_matrix[i][j]
            rv = hypergeom(m, n, N)
            tmp_dist = rv.pmf(x)
            if tmp_dist < best_dist:
                best_dist = tmp_dist
                best_taxa = j
        stats[i]["best_dist"] = best_dist
        stats[i]["best_taxa"] = best_taxa

    return stats


# This class stores the blocks graph that were founded in the bipartite graph.
# In this graph each vertex is corresponding to gene block. We connect two vertices in this graph only if they share a few genes from A in common
# Eventually we seek for cliques in the graph
def max_ranking_scoreue(clique, bip_graph):
    ranking_score = 0
    for block in clique:
        if bip_graph.bic_results[block].ranking_score > ranking_score:
            ranking_score = bip_graph.bic_results[block].ranking_score
    return ranking_score


def check_double_context_switch(block, bip_graph):
    tmp_strain_dict = {}
    check_context_set = set()
    for b in block.B:
        strain = bip_graph.B[b].organism
        if strain in tmp_strain_dict:
            check_context_set.add(strain)
        tmp_strain_dict[strain] = True
    return len(check_context_set) > 1


def check_context_switch(block, bip_graph):
    tmp_strain_dict = {}
    for b in block.B:
        strain = bip_graph.B[b].organism
        if strain in tmp_strain_dict:
            return True
        tmp_strain_dict[strain] = True
    return False


class BlockGraph(object):
    def __init__(self):
        self.best_cliques = []
        self.blocksGraph = nx.Graph()
        self.cliques = []
        self.super_cliques = []

    def init_block_graph(self, bip_graph, filter_context_switch):
        logging.info("Init block graph ")
        i = 0
        for bic in bip_graph.bic_results:
            add_block = True
            if filter_context_switch:
                add_block = check_context_switch(bic, bip_graph)
            # logging.info( add_block
            if add_block:
                self.blocksGraph.add_node(i)
                self.blocksGraph[i]['A'] = bic.get_a()
                self.blocksGraph[i]['B'] = bic.get_b()
                self.blocksGraph[i]['ranking_score'] = bic.get_ranking_scoreue()
                self.blocksGraph[i]['hyperGeoArray'] = []
                for d in bic.hyperDist:
                    dist = {'bestTaxa': bic.index_to_taxa_Dict[d['best_taxa']], 'bestDist': d['best_dist']}
                    self.blocksGraph[i]['hyperGeoArray'].append(dist)
                i += 1
        subsets = (itertools.combinations(self.blocksGraph.nodes(), 2))
        for sub in subsets:
            tmp_c = list(set(self.blocksGraph[sub[0]]['A']) & set(self.blocksGraph[sub[1]]['A']))
            if len(tmp_c) > 1:
                self.blocksGraph.add_edge(sub[0], sub[1], weight=len(tmp_c))

        self.cliques = list(nx.find_cliques(self.blocksGraph))

        # Create the cliques graph, Here each node is a clique of blocks from the blocks graph.
        # We connect two nodes in this graph only if they share more than 75% of their query subsequence genes
        # Eventually we seek for cliques of cliques of blocks and from each such clique we will pick in the greedy way
        # the best clique (for now we define the best clique to be the one that contains the most blocks.
        self.best_cliques = []
        cliques_genes_set = []
        cliques_graph = nx.Graph()
        i = 0
        # Create the nodes of the new graph, reminder. Each node is corresponding to a clqiue in the blocks graph
        for clique in self.cliques:
            cliques_graph.add_node(i)
            genes_set = Set([])
            for block in clique:
                for gene in self.blocksGraph[block]['A']:
                    genes_set.add(gene)
            cliques_genes_set.append(genes_set)
            i += 1
        # Now we add the edges in the cliques graph.
        for i in range(0, len(self.cliques)):
            for j in range(i + 1, len(self.cliques)):
                genes_set_a = cliques_genes_set[i]
                genes_set_b = cliques_genes_set[j]
                intersection = genes_set_a & genes_set_b
                if len(intersection) > len(genes_set_a) * 0.75 or len(intersection) > len(genes_set_b) * 0.75:
                    cliques_graph.add_edge(i, j)
        # Finding cliques of cliques
        self.super_cliques = list(nx.find_cliques(cliques_graph))
        # logging.info( self.super_cliques

        # From each clique of cliques we pick the best clique
        for superClique in self.super_cliques:
            best_clique = copy.deepcopy(superClique[0])
            for clique in superClique:
                if len(self.cliques[clique]) > len(self.cliques[best_clique]):
                    best_clique = copy.deepcopy(clique)
            if best_clique not in self.best_cliques:
                self.best_cliques.append(best_clique)
        tmp_cliques = []

        for i in self.best_cliques:
            # In this method we calculate for each calculate for eqach block its Hypergeometric distribution against all the taxonomies
            tmp_cliques.append(self.cliques[i])
        self.cliques = copy.deepcopy(tmp_cliques)

    def cliques_to_json(self, outfile, bip_graph):
        data_target = {'Organisms': []}
        org_to_index = {}
        for i in range(len(bip_graph.color_to_organims)):
            data_target['Organisms'].append(bip_graph.color_to_organims[i])
            org_to_index[bip_graph.color_to_organims[i]] = i
        with open(outfile + '_target.json', 'w') as outfile1:
            json.dump(data_target, outfile1)
        data_query = {'A': []}
        gene_name_to_number = {}
        for a in bip_graph.A:
            gene_name_to_number[a.to_json()['name']] = a.to_json()['id']
            data_query['A'].append(a.to_json())
        with open(outfile + '_query.json', 'w') as outfile2:
            json.dump(data_query, outfile2)
        with open(outfile + '_query_info_file.csv', 'w') as query_csv_file:
            query_csv_writer = csv.writer(query_csv_file)
            query_csv_writer.writerow(['Number', 'Start', 'End', 'Gene Id', 'Gene Name', 'Attribute', 'Strand'])
            for query_gene in data_query['A']:
                query_csv_writer.writerow(
                    [query_gene['id'] + 1, query_gene['start'], query_gene['stop'], query_gene['gene_id'],
                     query_gene['name'], query_gene['attribute'], query_gene['strand']])
        data_blocks = {'Minimum_Number_Genomes_In_Bicluster': str(MIN_NUM_OF_GENOMES_IN_CLUSTER + 1),
                       'Minimum_Number_Genes_Interval': str(MIN_NUM_OF_GENES_IN_INTERVAL + 1),
                       'Window_Size': str(WINDOW_SLIDE_SIZE), 'Genes_Gap': str(GENES_GAP),
                       'Number_Of_Cliques': str(len(self.cliques)), 'Cliques': []}
        i = 0
        self.cliques = sorted(self.cliques, key=lambda k: max_ranking_scoreue(k, bip_graph), reverse=True)
        for clique in self.cliques:
            clique = sorted(clique, key=lambda j: bip_graph.bic_results[j].ranking_score, reverse=True)
            bol = 0
            for block in clique:
                if bip_graph.bic_results[block].ranking_score > MIN_RANKING_SCORE:
                    bol = 1
            if bol == 1:
                data_blocks['Cliques'].append({'id': i})
                data_blocks['Cliques'][i]['Metadata'] = clique
                data_blocks['Cliques'][i]['Blocks'] = []
                for node in clique:
                    # logging.info(  self.blocksGraph[node]['hyperGeoArray']
                    data_blocks['Cliques'][i]['Blocks'].append(
                        bip_graph.bic_results[node].to_json(bip_graph, org_to_index,
                                                            self.blocksGraph[node]['hyperGeoArray']))
                i += 1
        with open(outfile + '.json', 'w') as outfile3:
            json.dump(data_blocks, outfile3)

        with open(outfile + '_results_file.csv', 'w') as results_csv_file:
            results_csv_writer = csv.writer(results_csv_file)
            for clique in data_blocks['Cliques']:
                results_csv_writer.writerow(['Clique Number', str(clique['id'] + 1)])
                block_index = 1
                for block in clique['Blocks']:
                    results_csv_writer.writerow([' '])
                    results_csv_writer.writerow(['Block Number', block_index])
                    results_csv_writer.writerow(['Ranking Score', block['ranking_scoreue']])
                    results_csv_writer.writerow(['Number of Genes from the Query genome', str(len(block['A']))])
                    results_csv_writer.writerow(['Number of Intervals from the Target genomes', str(len(block['B']))])
                    block_index += 1
                    results_csv_writer.writerow([' '])
                    results_csv_writer.writerow(['Group A (genes)'])
                    results_csv_writer.writerow(
                        ['Gene Number', 'Start', 'End', 'Gene Id', 'Gene Name', 'Attribute', 'Strand'])
                    for gene in block['A']:
                        results_csv_writer.writerow(
                            [gene['id'] + 1, gene['start'], gene['stop'], gene['gene_id'], gene['name'],
                             gene['attribute'], gene['strand']])
                    results_csv_writer.writerow([' '])
                    results_csv_writer.writerow(['Group B (Genes Intervals)'])
                    for index, gene_interval in enumerate(block['B']):
                        results_csv_writer.writerow([' '])
                        results_csv_writer.writerow(['Interval ' + str(index + 1)])
                        results_csv_writer.writerow(['Organism', gene_interval['organism']])
                        results_csv_writer.writerow(['Strain', gene_interval['strain']])
                        results_csv_writer.writerow(['Number of Genes',gene_interval['numOfGenes']])
                        results_csv_writer.writerow(
                            ['Gene Number', 'Start', 'End', 'Centroid Genome Gene Id', 'Centroid Genome Gene Name', 'Centroid Genome Attribute', 'Strand',
                             'Target Gene Id', 'Target Gene Attribute', 'Blast E-value'])
                        for gene in gene_interval['genes']:
                            results_csv_writer.writerow(
                                [gene_name_to_number[gene['name']] + 1, gene['start'], gene['stop'], gene['gene_id'],
                                 gene['name'], gene['attribute'], gene['strand'], gene['target_gene_name'],
                                 gene['target_gene_attribute'], gene['eval']])


def create_bipartite_graph(query_file, fname_list, refernce_folder):
    # Map each strain to its oragnism
    """

    :rtype : object
    """
    global NUMBER_OF_GENOMES

    bipartite_graph = BiPartiteGrpah()
    # Map each strain to its correspond organims and map each organism to its correspond color
    bipartite_graph.strain_to_organism = {}
    bipartite_graph.organism_to_color = {}
    bipartite_graph.color_to_organims = {}
    bipartite_graph.node_to_color = {}
    bipartite_graph.metadata = {}
    bipartite_graph.strain_to_length = {}
    color_index = 0
    logging.info("DEBUG: refernece folder " + refernce_folder)
    organisms = next(os.walk(refernce_folder))[1]
    organisms.sort()
    logging.info("After organisms.sort()")
    # Read the taxonomy of all the ref files
    with open('./TMP/taxonomy.json') as data_file:
        bipartite_graph.strain_to_taxonomy = json.load(data_file)
    for organism in organisms:
        bipartite_graph.organism_to_color[organism] = color_index
        bipartite_graph.color_to_organims[color_index] = organism
        color_index += 1
        strainFiles = next(os.walk(refernce_folder + '/' + organism))[2]
        for strain in strainFiles:
            bipartite_graph.strain_to_organism[strain.split(".")[0]] = organism
    geneId = 0
    genes_dict = {}
    gene_name_to_id = {}
    # Construct A.
    logging.info("Construct A")
    for line in [i.strip() for i in open(query_file).readlines()]:
        query_gene_name, query_gene_id, gene_start, gene_stop, query_gene_attribute, gene_strand = line.split("\t")
        gene = Gene(geneId, query_gene_name, query_gene_id, int(gene_start), int(gene_stop), gene_strand,
                    query_gene_attribute, "", "",
                    0)
        bipartite_graph.add_gene(gene)
        genes_dict[line.split("\t")[0]] = geneId
        gene_name_to_id[geneId] = line.split("\t")[0]
        geneId += 1
    # Now we finished to prepare all the nodes in A. In adddition we created the genes_dict that holds for each gene its unique id.
    bipartite_graph.set_d(geneId)
    geneId = 0
    interval_id = 0
    NUMBER_OF_GENOMES = len(fname_list)

    # Construct bipartite.B and create edges.
    logging.info("Construct B")
    for fname in fname_list:
        strain_id = (fname.split("/")[-1]).split(".")[0]
        # logging.info( bipartite_graph.strain_to_organism
        organism = bipartite_graph.strain_to_organism[strain_id]
        color = bipartite_graph.organism_to_color[organism]
        geneInterval = GeneInterval(strain_id, color, interval_id)
        lastGene = Gene(-1, '', '', -1000, -1000, 1, '', '', '', 0)
        lastEval = 100000
        bipartite_graph.strain_to_length[strain_id] = len([i.strip() for i in open(fname).readlines()])
        # logging.info( "Num of genes in genome " + str(organism) + " in strain " + str(strain_id) + " is: " + str(
        #     len([i.strip() for i in open(fname).readlines()]))
        for line in [i.strip() for i in open(fname).readlines()]:
            if line != "":
                gene_start, gene_stop, gene_strand, query_gene_name, query_gene_id, query_gene_attribute, e_val, a, b, gene_name, gene_attribute = line.split(
                    "|")
                gene = Gene(geneId, query_gene_name, query_gene_id, int(gene_start), int(gene_stop), gene_strand,
                            query_gene_attribute,
                            gene_name, gene_attribute, e_val)
                if abs(gene.start - lastGene.stop) < GENES_GAP or abs(
                                gene.start - lastGene.start) < GENES_OVERLAP or geneId == 0:
                    # check if the same gene in the reference was mapped to the same gene in operon, if so choose the one with the lowest evalue
                    if abs(gene.start - lastGene.start) < GENES_OVERLAP and geneId > 0:
                        if float(lastEval) > float(e_val):
                            genesLen = len(geneInterval.genes)
                            geneInterval.genes[genesLen - 1] = gene
                    else:
                        geneInterval.add_gene(gene)
                else:
                    # Quorom for the number of genes in Interval
                    if geneInterval.numOfGenes > MIN_NUM_OF_GENES_IN_INTERVAL:
                        # if geneInterval.numOfGenes > MIN_NUM_OF_GENES_IN_INTERVAL and geneInterval.numOfGenes < 0.90 * len(bipartite_graph.getA()):
                        geneInterval.taxonomy = bipartite_graph.strain_to_taxonomy[strain_id]['taxonomy']
                        bipartite_graph.add_gene_interval(geneInterval)
                        for geneI in geneInterval.genes:
                            bipartite_graph.add_edge(genes_dict[geneI.name], interval_id)
                        interval_id += 1
                    geneInterval = GeneInterval(strain_id, color, interval_id)
                    geneInterval.add_gene(gene)
                lastGene = gene
                lastEval = e_val
                geneId += 1
        # if MIN_NUM_OF_GENES_IN_INTERVAL < geneInterval.numOfGenes < 0.9 * len(bipartite_graph.getA()):
        if MIN_NUM_OF_GENES_IN_INTERVAL < geneInterval.numOfGenes:
            geneInterval.taxonomy = bipartite_graph.strain_to_taxonomy[strain_id]['taxonomy']
            bipartite_graph.add_gene_interval(geneInterval)
            for geneI in geneInterval.genes:
                bipartite_graph.add_edge(genes_dict[geneI.name], interval_id)
            interval_id += 1
        color += 1
    logging.info("Created bipartite graph with A= " + str(len(bipartite_graph.get_a())) + " B=" + str(
        len(bipartite_graph.get_b())) + " E=" + str(len(bipartite_graph.get_edges())))
    return bipartite_graph


# TODO here for each clique if it has at least one block with context switch, if yes add
def create_run_stats(bipartite_graph, block_graph):
    num_of_block = 0
    num_of_cliques = 0
    max_ranking_score = 10
    bol2 = 0
    check_context_switch_bol = False
    check_double_context_switch_bol = False
    for clique in block_graph.cliques:
        bol = 0
        for block in clique:
            check_double_context_switch_bol |= check_double_context_switch(bipartite_graph.bic_results[block],
                                                                           bipartite_graph)
            check_context_switch_bol |= check_context_switch(bipartite_graph.bic_results[block], bipartite_graph)
            # logging.info( 'real\ testy'
            # logging.info( check_context_switch_bol
            if bipartite_graph.bic_results[block].ranking_score > max_ranking_score:
                max_ranking_score = bipartite_graph.bic_results[block].ranking_score
                bol2 = 1
            if bipartite_graph.bic_results[block].ranking_score > MIN_RANKING_SCORE:
                bol = 1
        if bol == 1:
            num_of_cliques += 1
            num_of_block += len(clique)
    ans = {
        'context_switch': check_context_switch_bol,
        'double_context_switch': check_double_context_switch_bol,
        'numOfBlocks': num_of_block,
        'num_of_cliques': num_of_cliques,
        'avgBlockPerClique': 0,
        'accession': "",
        'length': -1,
        'number_of_genes': -1
    }
    if bol2 == 1:
        ans['max_ranking_score'] = max_ranking_score
    if ans['num_of_cliques'] != 0:
        ans['avgBlockPerClique'] = float(float(ans['numOfBlocks']) / float(ans['num_of_cliques']))
    return ans


def parallel_compute_biclusters_dict(query_file, blast_parse_folder, out_folder, reference_folder):
    bilcuster_out_folder = out_folder
    if not os.path.isdir(bilcuster_out_folder):
        os.makedirs(bilcuster_out_folder)
    fname_list = Utilities.return_recursive_dir_files(blast_parse_folder)

    bipartite_graph = create_bipartite_graph(query_file, fname_list, reference_folder)
    # In the case of on operon D is equal to the size of the operon, in the case where we divide the query into mulyiple operon regions we will have to set D for each one of them.
    bipartite_graph.calculate_biclusters()
    bipartite_graph.calculate_distance_matrices()
    logging.info("Finish create graph where |A|=" + str(len(bipartite_graph.A)) + " and |B|= " + str(
        len(bipartite_graph.A)) + " and |E|=" + str(len(bipartite_graph.E)))

    block_graph = BlockGraph()
    bol = False
    block_graph.init_block_graph(bipartite_graph, bol)
    block_graph.cliques_to_json(out_folder + query_file.split("/")[-1].split(".")[0], bipartite_graph)

    # bipartite_graph.parse_results(reference_folder, out_folder, query_file)

    return create_run_stats(bipartite_graph, block_graph)


def compute_bicluster(query_file, blast_parse_folder, out_folder, reference_folder, genome_size, min_genes_per_interval,
                      min_genomes_per_block, window_size, e_val, min_rank):
    global MIN_NUM_OF_GENOMES_IN_CLUSTER
    global WINDOW_SLIDE_SIZE
    global MIN_NUM_OF_GENES_IN_INTERVAL
    global MIN_RANKING_SCORE
    MIN_NUM_OF_GENOMES_IN_CLUSTER = min_genomes_per_block
    WINDOW_SLIDE_SIZE = window_size
    MIN_NUM_OF_GENES_IN_INTERVAL = min_genes_per_interval
    MIN_RANKING_SCORE = min_rank
    logging.info('Compute Biclusters ' + query_file + " " + blast_parse_folder + " " + out_folder + " " + reference_folder)
    global SIZE_OF_GENOME
    SIZE_OF_GENOME = genome_size
    logging.info("Genome size " + str(SIZE_OF_GENOME))

    start = time.time()
    stats_json = parallel_compute_biclusters_dict(query_file, blast_parse_folder, out_folder, reference_folder)

    logging.info(str(time.time() - start))
    return stats_json

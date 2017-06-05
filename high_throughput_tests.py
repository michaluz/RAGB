import logging
import time
import copy
import os
import math
import itertools
import networkx as nx
import json
import numpy
from scipy.stats import hypergeom
from sklearn.cluster import DBSCAN
from sets import Set


MIN_NUM_OF_GENOMES_IN_CLUSTER = 5
MIN_NUM_OF_GENES_IN_INTERVAL = 3
SIZE_OF_GENOME = 3200
GENES_GAP = 2000
GENES_OVERLAP = 100
MINIMUM_P_VALUE = 20


# the pipeline: CreateGraph: create a bipartite graph G=(A+B,E) where A={gi | gi gene at the centroid genome) A is order by the their gi.start. B = { geneInterval(i)| each geneInterval(i) is one sliding window of d genes
# from the reference genome where its id is been represented by the color of geneInterval(i)
def files_in_dir(root_dir):
    results = []
    for file in os.listdir(root_dir):
        if file.endswith(".ffc"):
            results.append(file)
    return results


def calculate_lengths(root_dir, results):
    genomes_length = {}
    for result in results:

        path = root_dir + "/" + result
        with open(path, 'r') as f:
            count = 0
            for line in f:
                if line[0] == ">":
                    count += 1
            genomes_length[result.split(".")[0]] = count
    return genomes_length


def compute_length(root_dir):
    results = files_in_dir(root_dir)
    return calculate_lengths(root_dir, results)


def calculate_subsets(array_length, edges_group,window):
    total_subsets = []
    if window < array_length:
        for i in range(window, array_length + 1):
            # logging.info("Calculate set for the " + str(i) + " sliding window")
            tmp_arr = []
            for edge in edges_group:
                if i - window <= edge < i:
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
            *(list(itertools.combinations(tmp_arr, n)) for n in range(MIN_NUM_OF_GENES_IN_INTERVAL, len(tmp_arr) + 1))))
        total_subsets.append(tmp_subsets)
    return total_subsets


def ncr(n, r):
    f = math.factorial
    return f(n) / f(r) / f(n - r)


# Calculate pvalue accordibng to the following equation ((n-k)C(d-k)/(nCd)*(n-d))^c
# n: size of genome
# d: length of the interval/operon in the query
# k: number of genes in A (in the operon)
# c: number of genomes in the bicluster
def calculate_log_pval(n, d, c, k):
    if d > n:
        d = n
    up = ncr(n - k, d - k) + 1
    down = ncr(n, d)
    tmp = float(float(up) / float(down))
    tmp = tmp * (n - d)
    pval = -(float(float(c) * math.log10(tmp)))
    return pval


def ld(s, t, direction_s, direction_t, query_set):
    s_prime = []
    t_prime = []
    direction_sprime = []
    direction_tprime = []
    i = 0
    for gene in s:
        if gene in query_set:
            s_prime.append(gene)
            direction_sprime.append(direction_s[i])
        i += 1
    i = 0
    for gene in t:
        if gene in query_set:
            t_prime.append(gene)
            direction_tprime.append(direction_t[i])
        i += 1
    s_prime[1:len(s_prime) + 1] = s_prime[0:len(s_prime)]
    s_prime[0] = ''
    direction_sprime[1:len(direction_sprime) + 1] = direction_sprime[0:len(direction_sprime)]
    direction_sprime[0] = ''

    t_prime[1:len(t_prime) + 1] = t_prime[0:len(t_prime)]
    t_prime[0] = ''
    direction_tprime[1:len(direction_tprime) + 1] = direction_tprime[0:len(direction_tprime)]
    direction_tprime[0] = ''
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
                if direction_sprime[i] == direction_tprime[j]:
                    d[i, j] = d[i - 1, j - 1]
                else:
                    d[i, j] = min(d[i - 1, j - 1] + 0.5, min(d[i - 1, j], d[i, j - 1]) + 1)
            else:
                d[i, j] = min(d[i - 1, j], d[i, j - 1], d[i - 1, j - 1]) + 1
    return d[s - 1, t - 1]


# Edit distance with all the genes, not only the query

def ld2(s, t, direction_s, direction_t):
    s_prime = []
    t_prime = []
    direction_sprime = []
    direction_tprime = []
    i = 0
    for gene in s:
        s_prime.append(gene)
        direction_sprime.append(direction_s[i])
        i += 1
    i = 0
    for gene in t:
        t_prime.append(gene)
        direction_tprime.append(direction_t[i])
        i += 1
    s_prime[1:len(s_prime) + 1] = s_prime[0:len(s_prime)]
    s_prime[0] = ''
    direction_sprime[1:len(direction_sprime) + 1] = direction_sprime[0:len(direction_sprime)]
    direction_sprime[0] = ''

    t_prime[1:len(t_prime) + 1] = t_prime[0:len(t_prime)]
    t_prime[0] = ''
    direction_tprime[1:len(direction_tprime) + 1] = direction_tprime[0:len(direction_tprime)]
    direction_tprime[0] = ''
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
                if direction_sprime[i] == direction_tprime[j]:
                    d[i, j] = d[i - 1, j - 1]
                else:
                    d[i, j] = min(d[i - 1, j - 1] + 0.5, min(d[i - 1, j], d[i, j - 1]) + 1)
            else:
                d[i, j] = min(d[i - 1, j], d[i, j - 1], d[i - 1, j - 1]) + 1
    return d[s - 1, t - 1]


# type of A vertex. Each vertex correspond to a gene in the centroid genome. Each node contains: id, name, start postion and end postion.
class Gene(object):
    def __init__(self, gene_id, name, start, stop, strand, attribute, target_gene_name, target_gene_attribute, eval):
        """
        :rtype : object
        """
        self.id = gene_id
        self.name = name
        self.start = start
        self.stop = stop
        self.strand = strand
        self.attribute = attribute
        self.target_gene_name = target_gene_name
        self.target_gene_attribute = target_gene_attribute
        self.eval = eval

    def to_json(self):
        ans = {'id': self.id, 'name': self.name, 'start': self.start, 'stop': self.stop, 'strand': self.strand,
               'attribute': self.attribute, 'target_gene_name': self.target_gene_name,
               'target_gene_attribute': self.target_gene_attribute, 'eval': self.eval}
        return ans


# type of B vertex. Each vertex represents a group of genes which are close to each other from one genome.
# Each genome is been represeneted by a color.
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

    def play(self):
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


# This class wil represents a full biclusters where A will holds the genes name, B will holds the geneInterval
class full_bicluster(object):
    def __init__(self):
        self.labels = []
        self.A = []
        self.B = []
        self.pval = 0

    def set_pvalue(self, pval):
        self.pval = pval

    def get_pvalue(self):
        return self.pval

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
        ans = {'hyper_geo_array': hyper_geo_array, 'A': [], 'B': [], 'pvalue': self.pval}
        for a in self.A:
            ans['A'].append(bip_graph.A[a].to_json())
        i = 0
        for b in self.B:
            organism = bip_graph.strain_to_organism[bip_graph.B[b].organism]
            index = org_to_index[organism]
            ans['B'].append(bip_graph.B[b].to_json(organism, index, self.labels[i]))
            i += 1

        return ans


class BiPartiteGrpah(object):
    def __init__(self):
        self.bic_results = []
        self.pval_results = []
        self.A = []
        self.B = []
        self.E = {}
        self.d = 0

    def set_d(self, d):
        self.d = d

    def ge_d(self):
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

    # noinspection PyPep8Naming
    def calculate_biclutsers(self,window):
        aGroups = []
        total_subsets = calculate_subsets(len(self.A) + 1, self.E,window)
        for indexes_subsets in total_subsets:
            round = 1
            round += 1
            for subset in indexes_subsets:
                if subset != () and subset[0] in self.get_edges():
                    b_prime = copy.deepcopy(self.get_edges()[subset[0]])
                    for i in subset:
                        if i in self.get_edges():
                            if len(b_prime) > 0:
                                b_prime = list(set(b_prime) & set(self.get_edges()[i]))
                        else:
                            b_prime = []
                    bic_color = 0
                    if len(b_prime) > 0:
                        bic_color = self.count_colors(b_prime)
                    if bic_color > MIN_NUM_OF_GENOMES_IN_CLUSTER:
                        biclique = full_bicluster()
                        biclique.set_a(subset)
                        biclique.set_b(b_prime)
                        pval = calculate_log_pval(SIZE_OF_GENOME, window, bic_color, len(subset))
                        biclique.set_pvalue(pval)
                        if pval > MINIMUM_P_VALUE and self.check_bic(biclique) == 1:
                            self.bic_results.append(biclique)
                            aGroups.append(biclique.get_a())
    # Check for each biclique in the results array:
    # If it contains the new biclique (smaller or equal): return 0 else return 1
    # If the new biclique contain the biclique delete the biclique from the results array
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
                # print self.B[b]
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
                    tmp = float(
                        ld2(copy.deepcopy(instance), copy.deepcopy(instance2), copy.deepcopy(direction_instances[i]),
                            copy.deepcopy(direction_instances[j])))
                    tmp_direction = direction_instances[i]
                    inverse_instance = list(reversed(instance))
                    inverse_direction_instance = list(reversed(tmp_direction))
                    inverse_direction_instance[:] = [int(x) * -1 for x in inverse_direction_instance]
                    inverse_direction_instance[:] = [str(x) for x in inverse_direction_instance]
                    tmp2 = float(ld2(copy.deepcopy(inverse_instance), copy.deepcopy(instance2),
                                     copy.deepcopy(inverse_direction_instance), copy.deepcopy(direction_instances[j])))
                    bic.distanceMatrix[i].append(float(min(tmp, tmp2)))
                    j += 1
                i += 1
            db = DBSCAN(metric='precomputed', eps=1, min_samples=1).fit(bic.distanceMatrix)
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


# This class represents the graph of all the blocks that were founded in the bipartite graph.
# In this graph each vertex is corresponding to gene block that was founded in the bipartite graph.
# We connect two vertices in this graph only if they share a few genes from A in common
# Eventaully we seek for cliques in the graph
def max_pvalue(clique, bip_graph):
    pval = 0
    for block in clique:
        if bip_graph.bic_results[block].pval > pval:
            pval = bip_graph.bic_results[block].pval
    return pval


class BlockGraph(object):
    def __init__(self):
        self.best_cliques = []
        self.blocksGraph = nx.Graph()
        self.cliques = []
        self.super_cliques = []

    def init_block_graph(self, bip_graph):
        # print "Init block graph"
        i = 0
        for bic in bip_graph.bic_results:

            self.blocksGraph.add_node(i)
            self.blocksGraph[i]['A'] = bic.get_a()
            self.blocksGraph[i]['B'] = bic.get_b()
            self.blocksGraph[i]['pval'] = bic.get_pvalue()
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
        # print self.super_cliques

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
        # with open(outfile + '_target.json', 'w') as outfile1:
        #     json.dump(data_target, outfile1)

        # print "query"
        data_query = {'A': []}
        for a in bip_graph.A:
            data_query['A'].append(a.to_json())
        # with open(outfile + '_query.json', 'w') as outfile2:
        #     json.dump(data_query, outfile2)

        # print "blocks"
        data_blocks = {'Minimum_Number_Genomes_In_Bicluster': str(MIN_NUM_OF_GENOMES_IN_CLUSTER + 1),
                       'Minimum_Number_Genes_Interval': str(MIN_NUM_OF_GENES_IN_INTERVAL + 1),
                       'Window_Size': str(WINDOW_SLIDE_SIZE), 'Genes_Gap': str(GENES_GAP),
                       'Number_Of_Cliques': str(len(self.cliques)), 'Cliques': []}
        i = 0
        self.cliques = sorted(self.cliques, key=lambda k: max_pvalue(k, bip_graph), reverse=True)
        for clique in self.cliques:
            clique = sorted(clique, key=lambda j: bip_graph.bic_results[j].pval, reverse=True)
            bol = 0
            for block in clique:
                if bip_graph.bic_results[block].pval > 40:
                    bol = 1
            if bol == 1:
                data_blocks['Cliques'].append({'id': i})
                data_blocks['Cliques'][i]['Metadata'] = clique
                data_blocks['Cliques'][i]['Blocks'] = []
                for node in clique:
                    # print  self.blocksGraph[node]['hyperGeoArray']
                    data_blocks['Cliques'][i]['Blocks'].append(
                        bip_graph.bic_results[node].to_json(bip_graph, org_to_index,
                                                            self.blocksGraph[node]['hyperGeoArray']))
                i += 1
        # with open(outfile + '.json', 'w') as outfile3:
        #     json.dump(data_blocks, outfile3)


# noinspection PyPep8Naming
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result


# noinspection PyPep8Naming
def create_bipartite_graph(query_file, fname_list, refernce_folder,gap):
    # Map each strain to its oragnism
    bipartite_graph = BiPartiteGrpah()
    # Map each strain to its correspond organims and map each organism to its correspond color
    bipartite_graph.strain_to_organism = {}
    bipartite_graph.organism_to_color = {}
    bipartite_graph.color_to_organims = {}
    bipartite_graph.node_to_color = {}
    bipartite_graph.metadata = {}
    bipartite_graph.strain_to_length = {}
    color_index = 0
    organisms = next(os.walk('./' + refernce_folder))[1]
    organisms.sort()
    # Read the taxonomy of all the ref files
    with open('./TMP/taxonomy.json') as data_file:
        bipartite_graph.strain_to_taxonomy = json.load(data_file)
    for organism in organisms:
        bipartite_graph.organism_to_color[organism] = color_index
        bipartite_graph.color_to_organims[color_index] = organism
        color_index += 1
        strainFiles = next(os.walk('./' + refernce_folder + '/' + organism))[2]
        for strain in strainFiles:
            bipartite_graph.strain_to_organism[strain.split(".")[0]] = organism
    geneId = 0
    genes_dict = {}
    gene_name_to_id = {}
    # Construct A.
    # logging.info("Construct A")
    for line in [i.strip() for i in open(query_file).readlines()]:
        # print query_file + " " + line
        query_gene_name, gene_start, gene_stop, query_gene_attribute, gene_strand = line.split("\t")
        gene = Gene(geneId, query_gene_name, int(gene_start), int(gene_stop), gene_strand, query_gene_attribute, "", "",0)
        bipartite_graph.add_gene(gene)
        genes_dict[line.split("\t")[0]] = geneId
        gene_name_to_id[geneId] = line.split("\t")[0]
        geneId += 1
    # Now we finished to prepare all the nodes in A. In adddition we created the genes_dict that holds for each gene its unique id.
    bipartite_graph.set_d(geneId)
    geneId = 0
    interval_id = 0
    # Construct bipartite.B and create edges.
    # logging.info("Construct B")
    for fname in fname_list:
        strain_id = (fname.split("/")[-1]).split(".")[0]
        organism = bipartite_graph.strain_to_organism[strain_id]
        color = bipartite_graph.organism_to_color[organism]
        geneInterval = GeneInterval(strain_id, color, interval_id)
        lastGene = Gene(-1, "", -1000, -1000, 1, "", "", "", 0)
        lastEval = 1000000000
        bipartite_graph.strain_to_length[strain_id] = len([i.strip() for i in open(fname).readlines()])
        # print "Num of genes in genome " + str(organism) + " in strain " + str(strain_id) + " is: " + str(
        #     len([i.strip() for i in open(fname).readlines()]))
        for line in [i.strip() for i in open(fname).readlines()]:
            if line != "":
                gene_start, gene_stop, gene_strand, query_gene_name, query_gene_attribute, e_val, a, b, gene_name, gene_attribute = line.split(
                    "|")
                gene = Gene(geneId, query_gene_name, int(gene_start), int(gene_stop), gene_strand, query_gene_attribute,
                            gene_name, gene_attribute, e_val)
                if abs(gene.start - lastGene.stop) < gap or abs(
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
    # logging.info("Created bipartite graph with A= " + str(len(bipartite_graph.get_a())) + " B=" + str(
    #     len(bipartite_graph.get_b())) + " E=" + str(len(bipartite_graph.get_edges())))
    # print query_file + " 4.5 "
    return bipartite_graph


def create_run_stats(bipartite_graph, block_graph):
    num_of_block = 0
    num_of_cliques = 0
    max_pval = 40
    bol2 = 0
    for clique in block_graph.cliques:
        bol = 0
        for block in clique:
            if bipartite_graph.bic_results[block].pval > max_pval:
                max_pval = bipartite_graph.bic_results[block].pval
                bol2 = 1
            if bipartite_graph.bic_results[block].pval > 40:
                bol = 1
        if bol == 1:
            num_of_cliques += 1
            num_of_block += len(clique)
    ans = {
        'numOfBlocks': num_of_block,
        'num_of_cliques': num_of_cliques,
        'avgBlockPerClique': 0,
        'accession': "",
        'length': -1,
        'number_of_genes': -1
    }
    if bol2 == 1:
        ans['max_pval'] = max_pval
    if ans['num_of_cliques'] != 0:
        ans['avgBlockPerClique'] = float(float(ans['numOfBlocks']) / float(ans['num_of_cliques']))
    return ans


def parallel_compute_biclusters_dict(query_file, blast_parse_folder, out_folder, reference_folder,gap,window):
    # print query_file + " 3 "
    bilcuster_out_folder = out_folder
    if not os.path.isdir(bilcuster_out_folder):
        os.makedirs(bilcuster_out_folder)
    fname_list = returnRecursiveDirFiles(blast_parse_folder)
    # print query_file + " 4 "
    bipartite_graph = create_bipartite_graph(query_file, fname_list, reference_folder,gap)
    # In the case of on operon D is equal to the size of the operon, in the case where we divide the query into mulyiple operon regions we will have to set D for each one of them.
    # print query_file + " 5 "
    bipartite_graph.calculate_biclutsers(window)
    bipartite_graph.calculate_distance_matrices()

    block_graph = BlockGraph()
    block_graph.init_block_graph(bipartite_graph)
    # block_graph.cliques_to_json(out_folder + query_file.split("/")[-1].split(".")[0], bipartite_graph)

    # bipartite_graph.parse_results(reference_folder, out_folder, query_file)
    return create_run_stats(bipartite_graph, block_graph)


def compute_bicluster(tupleParams):
    query_file, blast_parse_folder, out_folder, reference_folder = tupleParams
    # print "start compute_bicluster on " + query_file
    gapTimesArray = []
    windowsTimesArray = []
    # print query_file + " 1 "
    windowSize = list(range(5,10))

    gap = 2000
    for window in windowSize:
        start = time.time()
        parallel_compute_biclusters_dict(query_file, blast_parse_folder, out_folder, reference_folder,gap,window)
        windowsTimesArray.append(time.time() - start)
        print time.time() - start
    testAns = {
        'min_log_pval' : MINIMUM_P_VALUE,
        'min_num_genes_in_interval' : MIN_NUM_OF_GENES_IN_INTERVAL,
        'min_num_genomes_in_cluster' : MIN_NUM_OF_GENOMES_IN_CLUSTER,
        'constant_gap': gap,
        'constant_window' : window,
        'windowsTimesArray' : [],
    }
    j = 0
    for i in windowSize:
        testAns['windowsTimesArray'].append({'window':i,'time':windowsTimesArray[j]})
        j += 1
    with open('./highThroughput/highThroughput-' +os.path.splitext(query_file)[0].split('/')[-1] +  '.json', 'w') as outfile2:
        json.dump(testAns, outfile2)
    return 0

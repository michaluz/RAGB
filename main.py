from multiprocessing import Pool
import time
import os
import sys
import argparse
import shutil
import parse_input
import blast_script
import blast_parse
import compute_bicliques
import high_throughput_tests
import json
import logging
import random
import csv
from utilities import Utilities

__author__ = 'Arnon Benshahar'


def parser_code():
    parser = argparse.ArgumentParser(
        description='The purpose of this script is to run the full software suite that we have developed to study gene clusters.')

    parser.add_argument("-q", "--qfolder", dest="qfolder", metavar="DIRECTORY", default='./data/testing/iv/',
                        help="Folder containing the gbk or IslandViewer format files of the centroid query (or queries in case of multiple run).")

    parser.add_argument("-g", "--dbfolder", dest="dbfolder", metavar="DIRECTORY", default='./data/res/genomes/',
                        help="Folder containing all genbank files for use by the program as the reference genomes.")

    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="DIRECTORY", default='./results' + '/',
                        help="Folder where the results will be stored.")

    parser.add_argument("-d", "--window", dest="window_size", metavar="INT", default=12,
                        help="Size of our Biclustering algorithm's window")

    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default=os.sysconf("SC_NPROCESSORS_CONF"),
                        type=int,
                        help="Number of processors that you want this script to run on. The default is every CPU that the system has.")

    parser.add_argument("-iv", "--island_viewer_format", dest="island_viewer_format", metavar="STRING", default='F',
                        help="IslandViewer queries format, T for islandviewer format and F for normal gbk file.")

    parser.add_argument("-min_genomes", "--min_genomes_per_block", dest="min_genomes_per_block", metavar="INT",
                        default=4,
                        help="Minimum genome in a gene-block.")

    parser.add_argument("-min_genes", "--min_genes_per_interval", dest="min_genes_per_interval", metavar="INT",
                        default=3,
                        help="Minimum genes in a gene interval.")

    parser.add_argument("-rank", "--min_rank", dest="min_rank", metavar="INT", default=20,
                        help="Ranking score threshold")

    parser.add_argument("-e", "--e-val", dest="e_value", metavar="FLOAT", default='0.01',
                        help="E-value threshold for the BLAST search.")
    return parser.parse_args()


def check_options(parsed_args):
    if os.path.isdir(parsed_args.dbfolder):
        db_folder = parsed_args.dbfolder
    else:
        logging.info("The folder %s does not exist." % parsed_args.dbfolder)
        sys.exit()

    if os.path.isdir(parsed_args.qfolder):
        query_folder = parsed_args.qfolder
    else:
        logging.info("The folder %s does not exist." % parsed_args.qfolder)
        sys.exit()

    # if the directory that the user specifies does not exist, then the program makes it for them.
    if not os.path.isdir(parsed_args.outfolder):
        os.makedirs(parsed_args.outfolder)
    if parsed_args.outfolder[-1] != '/':
        out_folder = parsed_args.outfolder + '/'
    else:
        out_folder = parsed_args.outfolder

    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = int(parsed_args.num_proc)

    # validate the input for the window size
    try:
        window_size = int(parsed_args.window_size)
        if window_size <= 0:
            logging.info("The window that you entered %s is a negative number, please enter a positive integer." % parsed_args.max_gap)
            sys.exit()
        else:
            pass
    except Exception:
        logging.info("The window that you entered %s is not an integer, please enter a positive integer." % parsed_args.max_gap)
        sys.exit()

    # validate the query input format (isalndviewer or gbk)
    if parsed_args.island_viewer_format == 'F' or parsed_args.island_viewer_format == 'T':
        island_viewer_format = (parsed_args.island_viewer_format == 'T')
    else:
        logging.info("T for isalndviewer format and F for normal gbk format")
        sys.exit()

    # validate the input for the min_genomes_per_block
    try:
        min_genomes_per_block = int(parsed_args.min_genomes_per_block)
        if min_genomes_per_block <= 1:
            logging.info("The minimum genomes per block that you entered %s is less than 2, please enter a positive integer greater than 2." % parsed_args.max_gap)
            sys.exit()
        else:
            pass
    except:
        logging.info("The minimum genomes per block you entered %s is not an integer, please enter a positive integer." % parsed_args.max_gap)
        sys.exit()

    # validate the input for the min_genomes_per_block
    try:
        min_genes_per_interval = int(parsed_args.min_genes_per_interval)
        if min_genes_per_interval <= 1:
            logging.info( "The min genes per interval you entered %s is less than 2, please enter a positive integer greater than 2." % parsed_args.min_genes_per_interval)
            sys.exit()
        else:
            pass
    except:
        logging.info( "The min genes per interval you entered %s is not an integer, please enter a positive integer." % parsed_args.min_genes_per_interval)
        sys.exit()

    # validate the input for the min_genomes_per_block
    try:
        min_rank = int(parsed_args.min_rank)
        if min_rank <= 0:
            logging.info( "The min rank you entered %s is not an integer, please enter a positive integer." % parsed_args.min_rank)
            sys.exit()
        else:
            pass
    except:
        logging.info( "The min rank you entered %s is not an integer, please enter a positive integer." % parsed_args.min_rank)
        sys.exit()

    e_val = parsed_args.e_value
    return db_folder, query_folder, out_folder, num_proc, window_size, island_viewer_format, min_genes_per_interval, min_genomes_per_block, min_rank, e_val


def biclustering(tuple_list):
    logging.info( str(tuple_list))
    query_file, refernce_folder, ref_fasta, query_fasta, blast_results, blast_parse_dir, query_gene_list_dir, bicluster_results, max_genome_size, min_genes_per_interval, min_genomes_per_block,window_size, e_val, min_rank = tuple_list
    logging.info( "Parse blast results")
    list_file_name = query_gene_list_dir + query_file.split("/")[-1].split(".")[0] + ".txt"
    blast_parse.parse_blast(blast_results, blast_parse_dir, "", 10, list_file_name, query_file)

    logging.info( "Preforming our Biclustering Algorithm")
    return compute_bicliques.compute_bicluster(list_file_name, blast_parse_dir, bicluster_results, refernce_folder, max_genome_size, min_genes_per_interval, min_genomes_per_block,window_size, e_val, min_rank)


def main():
    start = time.time()
    logging.basicConfig(filename="./info.log",format='%(asctime)s %(message)s', level=logging.DEBUG)
    logging.info( "Start RAGBI program")
    # Parse all the user's arguments
    parsed_args = parser_code()
    db_folder, q_folder, outfolder, num_proc, window_size, island_viewer_format, min_genes_per_interval, min_genomes_per_block, min_rank, e_val = check_options(
        parsed_args)
    logging.info( "Command line arguments=" + str(check_options(parsed_args)))

    # Only for high_throughput
    high_throughput = False

    # This json will contain all the information about this run.
    stats_json = {'query_list': []}

    tmp_dir = './TMP'
    ref_fasta_dir = tmp_dir + '/ffc/ref/'
    query_fasta_dir = tmp_dir + '/ffc/qry/'

    if os.path.exists(outfolder):
        shutil.rmtree(outfolder)
    logging.info("Create output folder, " + outfolder )
    os.makedirs(outfolder)

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    logging.info("Create tmp folder, " + tmp_dir)
    os.makedirs(tmp_dir)

    if os.path.exists(ref_fasta_dir):
        shutil.rmtree(ref_fasta_dir)
    logging.info("Create ref fasta folder, " + ref_fasta_dir)
    os.makedirs(ref_fasta_dir)

    if os.path.exists(query_fasta_dir):
        shutil.rmtree(query_fasta_dir)
    logging.info("Create query fasta folder, " + query_fasta_dir)
    os.makedirs(query_fasta_dir)

    '''
    Convert the gbk files of the query and the reference to ffc format.
    In case the query is in normal gbk format (one file, one island) we convert the gbk file into one ffc file which contains the island.
    In case the query is in Islandviewer format we split each genebank file into multiple ffc files where each file contains one island.
    For more details about the input formats please go to the README.md file.
    '''
    if island_viewer_format:
        logging.info("Parse query in IslandViewer format," + query_fasta_dir)
        query_json = parse_input.parse_islandviewer(q_folder, query_fasta_dir)
    else:
        logging.info("Parse query in noraml format," + query_fasta_dir)
        query_json = parse_input.parse_gbk(q_folder, query_fasta_dir, "NONE", True, False)
    logging.info("Parse target genomes," + ref_fasta_dir)

    with open(outfolder + 'queries.json', 'w') as outfile1:
        json.dump(query_json, outfile1)

    with open(outfolder + 'centroid_genomes_summary.csv','w') as queries_csv_file:
        queries_csv_writer = csv.writer(queries_csv_file)
        if island_viewer_format:
            for q in query_json:
                queries_csv_writer.writerow(['Accession Number','Description','Number of islands','Length'])
                queries_csv_writer.writerow([q['accession'],q['description'],q['num_of_islands'],q['length']])
                queries_csv_writer.writerow(['Islands'])
                queries_csv_writer.writerow(['Start','End','Length', 'Number of Genes'])
                for query in q['islands']:
                    queries_csv_writer.writerow([query['start'],query['end'],query['length'],query['num_of_genes']])
                queries_csv_writer.writerow([])
        else :
            queries_csv_writer.writerow(['Organism','Accession Number','Description','Number of genes','Length'])
            for q in query_json:
                queries_csv_writer.writerow([q['organism'],q['accession'],q['description'],q['number_of_genes'],q['length']])

    target_json = parse_input.parse_gbk(db_folder, ref_fasta_dir, "NONE", True, True)

    with open(outfolder + 'targets.json', 'w') as outfile1:
        json.dump(target_json, outfile1)

    with open(outfolder + 'targets_summary.csv','w') as queries_csv_file:
        queries_csv_writer = csv.writer(queries_csv_file)
        queries_csv_writer.writerow(['Accession Number','Specie','Description','Length','Number of Genes'])
        for target in target_json:
            queries_csv_writer.writerow([target['accession'],target['organism'],target['description'],target['length'],target['number_of_genes']])

    # create the queries file
    blast_results_dir = tmp_dir + '/blast_results/'
    if os.path.exists(blast_results_dir):
        shutil.rmtree(blast_results_dir)
    logging.info("Create blast results folder, " + blast_results_dir)
    os.makedirs(blast_results_dir)

    blast_parse_tmp = tmp_dir + '/blast_parse/'
    if os.path.exists(blast_parse_tmp):
        shutil.rmtree(blast_parse_tmp)
    logging.info("Create blast parse folder, " + blast_parse_tmp)
    os.makedirs(blast_parse_tmp)

    query_gene_list_dir = tmp_dir + '/query_gene_list_dir/'
    if os.path.exists(query_gene_list_dir):
        shutil.rmtree(query_gene_list_dir)
    logging.info("Create query gene list folder, " + query_gene_list_dir)
    os.makedirs(query_gene_list_dir)

    query_fasta_list = []
    logging.info( 'Create blast parse folders' )
    for content in os.listdir(query_fasta_dir):  # "." means current directory
        if content.split(".")[-1] == "ffc":
            query_fasta_list.append(query_fasta_dir + content)
            if not os.path.exists(blast_results_dir + "" + content.split(".")[-2]):
                os.makedirs(blast_results_dir + "" + content.split(".")[-2])
            if not os.path.exists(blast_parse_tmp + "" + content.split(".")[-2]):
                os.makedirs(blast_parse_tmp + "" + content.split(".")[-2])

    logging.info( 'Create targets.json file' )
    with open(outfolder + 'targets.json') as data_file:
        targets_json = json.load(data_file)

    s = 0
    for target in targets_json:
        s += target['number_of_genes']
    genome_size = s / len(targets_json)

    logging.info( "Avg genome size " + str(genome_size) )

    general_stats = []
    file_num = 1
    tuple_list_array = []
    for file in query_fasta_list:
        blast_parse_tmp = tmp_dir + '/blast_parse/'
        blast_results_dir = tmp_dir + '/blast_results/'
        query_gene_list_dir = tmp_dir + '/query_gene_list_dir/'

        logging.info( 'File Number ' + str(file_num) )
        file_num += 1

        # Run blast with the query fasta vs the ref fastas
        query_file_name = file.split("/")[-1].split(".")[-2]
        blast_output = blast_results_dir + query_file_name + "/"
        logging.info("Run Blast Script on " + file)
        blast_script.blast(ref_fasta_dir, blast_output, "", os.sysconf("SC_NPROCESSORS_CONF"), file, e_val)

        blast_results_tmp = blast_results_dir + file.split("/")[-1].split(".")[-2] + "/"
        blast_parse_tmp = blast_parse_tmp + file.split("/")[-1].split(".")[-2] + "/"
        bicluster_results_tmp = outfolder + file.split("/")[-1].split(".")[-2] + "/"

        logging.info( str(file) + "," + str(db_folder) + "," + str(ref_fasta_dir) + "," + str(
            query_fasta_dir) + "," + str(blast_results_tmp) + "," + str(blast_parse) + "," + str(
            query_gene_list_dir) + "," + outfolder )

        tuple_list = (file, db_folder, ref_fasta_dir, query_fasta_dir, blast_results_tmp, blast_parse_tmp,
                      query_gene_list_dir, bicluster_results_tmp, genome_size, min_genes_per_interval, min_genomes_per_block, window_size, e_val, min_rank)

        logging.info( str(tuple_list))
        logging.info("Parse Blast results")
        list_file_name = query_gene_list_dir + file.split("/")[-1].split(".")[0] + ".txt"
        blast_parse.parse_blast(blast_results_tmp, blast_parse_tmp, "", 10, list_file_name, file)

        logging.info("Preforming our Biclustering Algorithm")
        file_stats = compute_bicliques.compute_bicluster(list_file_name, blast_parse_tmp, bicluster_results_tmp, db_folder, genome_size, min_genes_per_interval, min_genomes_per_block,window_size, e_val, min_rank)

        if not high_throughput:
            file_stats['accession'] = query_file_name
            with open(outfolder + 'queries.json') as data_file:
                query_json = json.load(data_file)
            for query in query_json:
                if query['accession'] == query_file_name:
                    file_stats['length'] = query['length']
                    file_stats['number_of_genes'] = query['number_of_genes']
            if file_stats['num_of_cliques'] > 0:
                general_stats.append(file_stats)
                with open(outfolder + 'resultStats.json', 'w') as outfile1:
                    json.dump(general_stats, outfile1)
        else:
            tuple_list_array.append((query_gene_list_dir + file.split("/")[-1].split(".")[0] + ".txt", blast_parse_tmp, './', db_folder))

    with open(outfolder + 'general_results.csv','w') as general_results_csv_file:
        general_results_csv_writer = csv.writer(general_results_csv_file)
        general_results_csv_writer.writerow(['Arguments'])
        general_results_csv_writer.writerow(["Window's Size",window_size,"Minimum Number of Genes per Interval",min_genes_per_interval,"Minimum Number of Genomes per Block", min_genomes_per_block, "Ranking Score Threshold", min_rank, "E-value Threshold", e_val])
        general_results_csv_writer.writerow([' '])
        if island_viewer_format:
            general_results_csv_writer.writerow(['Accession Number- Start of Island - End of Island','Max Ranking Score','Number of Cliques','Number of Blocks','Avg Blocks per Clique','Context Switch'])
        else:
            general_results_csv_writer.writerow(['Accession Number','Max Ranking Score','Number of Cliques','Number of Blocks','Avg Blocks per Clique','Context Switch'])
        for result in general_stats:
            general_results_csv_writer.writerow([result['accession'],result['max_ranking_score'],result['num_of_cliques'],result['numOfBlocks'],result['avgBlockPerClique'],result['context_switch']])

    for root, dirs, files in os.walk(outfolder):
        for file in files:
            if file.endswith('.json'):
                os.remove(root + '/' + file)

    shutil.rmtree(tmp_dir)

    logging.info("Run time: " + str(time.time() - start))


if __name__ == '__main__':
    main()

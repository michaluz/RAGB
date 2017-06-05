import shutil
import json
import time
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.SeqRecord import SeqRecord
from utilities import  Utilities
import uniprot as uni
import csv
import logging

__author__ = 'Arnon Benshahar'

# This function will check if the gene is locate in one of the islands. If it is, it returns an array with all the occurrences of the gene. Each occurrence is represented by the start and the end of the island, else it returns -1
def find_island(start_end_array, island_list, start, end, gene_name, gene_locus):
    locations_array = []
    for i_start, i_end in start_end_array:
        if int(i_start) - 2000 <= int(start) and int(i_end) + 2000 >= int(end):
            genes = island_list[str(i_start + '-' + i_end)]['genes']
            loci = island_list[str(i_start + '-' + i_end)]['loci']
            for i in range(len(genes)):
                if genes[i] == gene_name:
                    locations_array.append([i_start, i_end])
                elif loci[i] == gene_locus:
                    locations_array.append([i_start, i_end])
    if len(locations_array) > 0:
        return locations_array
    else:
        return -1


# This function will convert input file of IslandViewer format into multiple ffc file (for more information about the input format please go to README.md).
def convert_islandviewer(genbank_path, csv_path, db_directory):
    island_list = {}
    start_end_array = []
    with open(csv_path, 'rb') as f:
        # parse the csv file
        reader = csv.reader(f)
        for row in reader:
            if str(row[0]) != 'Island start':
                start, end, length, method, gene_name, gene_id, locus, gene_start, gene_end, strand, product, extranl_annotations = row
                if str(start + '-' + end) in island_list:
                    if gene_name != "":
                        island_list[str(start + '-' + end)]['genes'].append(gene_name.replace('_', ''))
                        island_list[str(start + '-' + end)]['loci'].append(locus.replace('_', ''))
                    else:
                        island_list[str(start + '-' + end)]['genes'].append(gene_id.replace('_', ''))
                        island_list[str(start + '-' + end)]['loci'].append(locus.replace('_', ''))
                else:
                    start_end_array.append([start, end])
                    island_list[str(start + '-' + end)] = {}
                    island_list[str(start + '-' + end)]['genes'] = []
                    island_list[str(start + '-' + end)]['loci'] = []
                    if gene_name != "":
                        island_list[str(start + '-' + end)]['genes'].append(gene_name.replace('_', ''))
                        island_list[str(start + '-' + end)]['loci'].append(locus.replace('_', ''))
                    else:
                        island_list[str(start + '-' + end)]['genes'].append(gene_id.replace('_', ''))
                        island_list[str(start + '-' + end)]['loci'].append(locus.replace('_', ''))

        # parse the gbk file
        record_list = {}
        islands_from_query = {}
        seq_record = SeqIO.parse(open(genbank_path), "genbank").next()
        accession = seq_record.annotations['accessions'][0]
        organism = seq_record.annotations['organism'].replace(' ', '_')
        islands_from_query['accession'] = seq_record.annotations['accessions'][0]
        islands_from_query['organism'] = seq_record.annotations["organism"]
        islands_from_query['length'] = len(seq_record)
        islands_from_query['description'] = seq_record.description
        num_of_cds = 0
        for feature_number, feature in enumerate(seq_record.features):
            err_flag = False
            error_in_field = False
            if feature.type == 'CDS' and ('protein_id' in feature.qualifiers or 'gene' in feature.qualifiers):
                start = feature.location.start
                end = feature.location.end
                try:
                    locus = feature.qualifiers['locus_tag'][0].replace('_', '')
                except:
                    try:
                        locus = feature.qualifiers['old_locus_tag'][0].replace('_', '')
                    except:
                        locus = 'unknown'
                try:
                    old_locus = feature.qualifiers['old_locus_tag'][0].replace('_', '')
                except:
                    try:
                        old_locus = feature.qualifiers['old_locus_tag'][0].replace('_', '')
                    except:
                        old_locus = 'unknown'
                try:
                    gene = feature.qualifiers['protein_id'][0]
                    locations = find_island(start_end_array, island_list, start, end,
                                            feature.qualifiers['protein_id'][0].replace('_', ''), old_locus)
                except:
                    gene = feature.qualifiers['gene'][0]
                    locations = find_island(start_end_array, island_list, start, end,
                                            feature.qualifiers['gene'][0].replace('_', ''), old_locus)
                if locations != -1:
                    for location in locations:
                        i_start, i_end = location[0], location[1]

                        # if 'locus_tag' in feature.qualifiers and find_island(start_end_array, island_list, start, end, feature.qualifiers['protein_id'][0].replace('_','')) != -1:
                        #     location = find_island(start_end_array, island_list, start, end,
                        #                            feature.qualifiers['locus_tag'][0].replace('_', ''))
                        #     locus = feature.qualifiers['locus_tag'][0].replace('_', '')
                        #     if location != -1:
                        #         i_start, i_end = location[0], location[1]
                        # elif 'old_locus_tag' in feature.qualifiers and find_island(start_end_array, island_list, start, end, feature.qualifiers['protein_id'][0].replace('_', '')) != -1:
                        #     location = find_island(start_end_array, island_list, start, end,
                        #                            feature.qualifiers['old_locus_tag'][0].replace('_', ''))
                        #     locus = feature.qualifiers['old_locus_tag'][0].replace('_', '')
                        #     if location != -1:
                        #         i_start, i_end = location[0], location[1]

                        try:
                            start = int(feature.location.start)
                            stop = int(feature.location.end)
                        except:
                            error_in_field = True

                        strand = feature.strand
                        try:
                            product = feature.qualifiers['product'][0]
                        except:
                            product = "unknown"

                        try:
                            note = feature.qualifiers['note'][0]
                        except:
                            note = 'unknown'

                        try:
                            gene_real_name = feature.qualifiers['gene'][0]
                        except:
                            gene_real_name = 'unknown'

                        if 'translation' in feature.qualifiers.keys():
                            prot_seq = Seq(''.join(feature.qualifiers['translation']), IUPAC.protein)
                            # print "prot_seq", type(prot_seq), prot_seq
                            dna_seq = seq_record.seq[start:stop]
                            # print "dna_seq", type(dna_seq), dna_seq
                            gc = GC(dna_seq)
                            gc = "%3.2f" % gc

                            seq_rec_to_store = SeqRecord(prot_seq, id='|'.join(
                                [accession, organism, locus, gene, gene_real_name, product, str(start), str(stop),
                                 str(strand),
                                 gc, old_locus]).replace(' ', '_'), description='')
                            if str(str(i_start) + '-' + str(i_end)) in record_list:
                                record_list[str(str(i_start) + '-' + str(i_end))].append(seq_rec_to_store)
                            else:
                                record_list[str(str(i_start) + '-' + str(i_end))] = []
                                record_list[str(str(i_start) + '-' + str(i_end))].append(seq_rec_to_store)

                                # elif 'protein_id' in feature.qualifiers:
                                #     gene = feature.qualifiers['protein_id'][0]
                                #     # record_list.append(SeqRecord(prot_seq, id = '|'.join([accession, organism, locus, gene, str(start), str(stop), str(strand), gc]).replace(' ', ''), description = ''))
                                #     seq_rec_to_store = SeqRecord(prot_seq, id='|'.join(
                                #         [accession, organism, locus, gene, description, str(start), str(stop), str(strand),
                                #          gc]).replace(' ', '_'), description='')
                                #     if str(str(i_start) + '-' + str(i_end)) in record_list:
                                #         record_list[str(str(i_start) + '-' + str(i_end))].append(seq_rec_to_store)
                                #     else:
                                #         record_list[str(str(i_start) + '-' + str(i_end))] = []
                                #         record_list[str(str(i_start) + '-' + str(i_end))].append(seq_rec_to_store)
                                # else:
                                #     print "No name for gene"
        islands_from_query['num_of_islands'] = len(record_list)
        islands_from_query['islands'] = []
        for key in record_list:
            if len(record_list[key]) == len(island_list[key]['genes']):
                islands_from_query['islands'].append({
                    'start': key.split('-')[0],
                    'end': key.split('-')[1],
                    'length': int(key.split('-')[1]) - int(key.split('-')[0]),
                    'num_of_genes': len(record_list[key])
                })
                outpath = db_directory + islands_from_query['accession'] + '-' + key + '.ffc'
                print outpath
                logging.info("Create new ffc file " + islands_from_query['accession'] + '-' + key + '.ffc')
                out_handle = open(outpath, "w")
                SeqIO.write(record_list[key], out_handle, "fasta")
                out_handle.close()
        return islands_from_query


# take the genbank file specified by genbank path, and save the customized result file in the db_directory folder
def convert_genbank(genbank_tuple, create_tax):
    global seq_rec_to_store, err_flag
    query_json = []
    accession_to_taxonomy = {}
    for tup in genbank_tuple:
        query = {}
        logging.info("Parsing tuple " + str(tup))
        print "Parsing tuple " + str(tup)
        genbank_path, db_directory, error_fname, do_protein = tup
        record_list = []
        seq_record = SeqIO.parse(open(genbank_path), "genbank").next()
        accession = seq_record.annotations['accessions'][0]
        organism = seq_record.annotations['organism'].replace(' ', '_')
        query['accession'] = seq_record.annotations['accessions'][0]
        query['organism'] = seq_record.annotations["organism"]
        query['length'] = len(seq_record)
        query['description'] = seq_record.description

        tmp = {'taxonomy': seq_record.annotations["taxonomy"][2], 'organism': genbank_path.split('/')[:-1][-1]}
        accession_to_taxonomy[query['accession']] = tmp

        err_log = []
        gc_list = []
        num_of_cds = 0
        # loop over the genbank file
        for fnum, feature in enumerate(seq_record.features):
            err_flag = False
            error_in_field = False
            if feature.type == 'CDS':
                num_of_cds += 1
                try:
                    start = int(feature.location.start)
                    stop = int(feature.location.end)
                except:
                    error_in_field = True

                strand = feature.strand
                dna_seq = seq_record.seq[start:stop]
                # print "dna_seq", type(dna_seq), dna_seq
                gc = GC(dna_seq)
                gc_list.append(gc)
                gc = "%3.2f" % gc
                try:
                    locus = feature.qualifiers['locus_tag'][0]
                except:
                    try:
                        locus = feature.qualifiers['gene'][0]
                    except:
                        locus = 'unknown'

                try:
                    description = feature.qualifiers['product'][0]
                except:
                    try:
                        description = feature.qualifiers['note'][0]
                    except:
                        description = 'unknown'

                if do_protein:
                    # seq = seq.translate()
                    # print type(seq)
                    # print feature.qualifiers.keys()
                    # seq = dir(feature)
                    try:
                        if 'translation' in feature.qualifiers.keys():
                            prot_seq = Seq(''.join(feature.qualifiers['translation']), IUPAC.protein)
                            # print "prot_seq", type(prot_seq), prot_seq

                            if 'gene' in feature.qualifiers:
                                gene = feature.qualifiers['gene'][0]
                                # record_list.append(SeqRecord(prot_seq, id = '|'.join([accession, organism, locus, gene, str(start), str(stop), str(strand), gc]).replace(' ', ''), description = ''))
                                seq_rec_to_store = SeqRecord(prot_seq, id='|'.join(
                                    [accession, organism, locus, gene, description, str(start), str(stop), str(strand),
                                     gc]).replace(' ', '_'), description='')
                            elif 'protein_id' in feature.qualifiers:
                                gene = feature.qualifiers['protein_id'][0]
                                # record_list.append(SeqRecord(prot_seq, id = '|'.join([accession, organism, locus, gene, str(start), str(stop), str(strand), gc]).replace(' ', ''), description = ''))
                                seq_rec_to_store = SeqRecord(prot_seq, id='|'.join(
                                    [accession, organism, locus, gene, description, str(start), str(stop), str(strand),
                                     gc]).replace(' ', '_'), description='')
                            else:
                                print "No name for gene"
                    except Exception:
                        err_log.append(Exception)
                    pass
                if not error_in_field:
                    record_list.append(seq_rec_to_store)
                else:
                    logging.info("A record was omitted")
        query['number_of_genes'] = num_of_cds
        query_json.append(query)
        for i in err_log:
            print str(i)
            print i
            logging.error('Error:' + i)
        if not err_flag:
            outpath = db_directory + os.path.splitext(os.path.basename(genbank_path))[0] + '.ffc'
            out_handle = open(outpath, "w")
            SeqIO.write(record_list, out_handle, "fasta")
            out_handle.close()
        if do_protein:
            cmd = "makeblastdb -in %s -dbtype prot" % outpath
        else:
            cmd = "makeblastdb -in %s -dbtype nucl" % outpath
        os.system(cmd)
    if create_tax:
        with open('./TMP/taxonomy.json', 'w') as outfile1:
            json.dump(accession_to_taxonomy, outfile1)
    return query_json


def parallel_convert_genbank(file_list, outfolder, do_protein, create_tax, error_fname="./error_log.txt"):
    if os.path.isfile(error_fname):
        os.remove(error_fname)
    # Package the variables for the convert_genbank function so everything can be run in parallel
    tuple_list = [(i, outfolder, error_fname, do_protein) for i in file_list]
    return convert_genbank(tuple_list, create_tax)


def parse_gbk(infolder, outfolder, filter_file, do_protein, create_tax):
    flist = Utilities.return_recursive_dir_files(infolder)

    if filter_file != 'NONE':
        filter_list = [i.strip() for i in open(filter_file).readlines()]
        file_list = [i for i in flist if i.split('/')[-1].split('.')[0] in filter_list]
    else:
        file_list = flist
    ans = parallel_convert_genbank(file_list, outfolder, do_protein, create_tax)
    return ans


def parse_islandviewer(in_path, ffc_path):
    queries = []
    for dir in [x[0] for x in os.walk(in_path)][1:]:
        accession = dir.split('/')[-1]
        csv_path = dir + '/' + accession + '.csv'
        genbank_path = dir + '/' + accession + '.gbk'
        queries.append(convert_islandviewer(genbank_path, csv_path, ffc_path))
    return queries

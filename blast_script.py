# Copyright(C) 2014 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

#!/usr/bin/python

from multiprocessing import Pool
import logging
import time
import os
import os, shutil
import sys
import argparse
from utilities import Utilities
from Bio.Blast.Applications import NcbiblastxCommandline

# This function will return all of the files that are in a directory. os.walk is recursive traversal.
# def return_recursive_dir_files(root_dir):
#     result = []
#     for path, dir_name, flist in os.walk(root_dir):
#         for f in flist:
#             fname = os.path.join(path, f)
#             if os.path.isfile(fname):
#                 result.append(fname)
#     return result


def do_parallel_blast(arg_tuple):
    db, query_file, blast_result_folder, num_processors, eval_threshold = arg_tuple
    if (db.split('/')[-1] != query_file.split('/')[-1]):
        out_file = "%s%s.txt" % (blast_result_folder, db.split('/')[-1].split('.')[0])
        print "Before blast query file: " + query_file + ", db:" + db
        logging.info("Before blast query file: " + query_file + ", db:" + db)
        cmd = "blastp -query %s -db %s -evalue %s -out %s -outfmt 6" % (query_file, db, float(eval_threshold), out_file)
        os.system( cmd )
        print "Finish  blast on query file: " + query_file + ", db:" + db
        logging.info("Finish  blast on query file: " + query_file + ", db:" + db)


def parallel_blast(database_folder, outfolder, filter_file, num_proc, query_file, e_val):
    #Delete all the current files in out folder
    import os, shutil
    logging.info('Parallerl blast')
    print 'Parallerl blast'
    print query_file
    folder = outfolder
    print folder
    for the_file in os.listdir(folder):
        file_path = os.path.join(folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            #elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            logging.info(Exception)
            print(e)
    
    unfiltered_db_list = [i for i in Utilities.return_recursive_dir_files(database_folder) if i.split('/')[-1].split('.')[-1] == 'ffc']

    if filter_file == '':
        db_list = unfiltered_db_list
    else:
        filter_list = [i.strip() for i in open(filter_file).readlines()]
        db_list = [i for i in unfiltered_db_list if i.split('/')[-1].split('.')[0] in filter_list]

    blast_arg_list = [(i, query_file, outfolder, 1, e_val) for i in db_list]
    pool = Pool(processes = num_proc)
    pool.map(do_parallel_blast, blast_arg_list)


def blast(database_folder, outfolder, filter_file, num_proc, query_file, e_val):
    
    start = time.time()

    print database_folder, outfolder, filter_file, num_proc, query_file, e_val
    
    parallel_blast(database_folder, outfolder, filter_file, num_proc, query_file, e_val)

    print time.time() - start
    


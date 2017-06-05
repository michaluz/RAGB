import os
import math
import logging
from Bio import SeqIO

__author__ = 'Arnon Benshahar'


class Utilities(object):

    @staticmethod
    def ffc_files_in_dir(root_dir):
        results = []
        for file in os.listdir(root_dir):
            if file.endswith(".ffc"):
                results.append(file)
        return results

    # This function will return all the files that are in a directory. os.walk is recursive traversal.
    @staticmethod
    def return_recursive_dir_files(root_dir):
        result = []
        for path, dir_name, flist in os.walk(root_dir):
            for f in flist:
                fname = os.path.join(path, f)
                if os.path.isfile(fname):
                    result.append(fname)
        return result

    @staticmethod
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

    @staticmethod
    def ncr(n, r):
        f = math.factorial
        return f(n) / f(r) / f(n - r)

    @staticmethod
    # Calculate ranking score according to the following equation ((n-k)C(d-k)/(nCd)*(n-d))^c
    # n: size of genome
    # d: length of the interval/operon in the query
    # k: number of genes in A (in the operon)
    # c: number of genomes in the bicluster
    def calculate_ranking_score(m, n, d, c, k):
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
        pval = -float(
            float(math.log10(m - c + 1)) + float(math.log10(Utilities.ncr(m, i_prime))) + float(i_prime * math.log10(p)) + float(
                (m - i_prime) * math.log10(1 - p)))
        # logging.info( "m=" + str(m) + " n=" + str(n) + " d=" + str(d) + " c=" + str(c) +  " k=" + str(k) + " m=" + str(m) + " up= " + str(up) + "  down= " + str(down) + " p=" + str(p) + " i_prime=" + str(i_prime) + " pval=" + str(pval)
        return pval


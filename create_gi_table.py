import json
import time
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.SeqRecord import SeqRecord
import uniprot as uni

__author__ = 'user'

# Maybe will be added in the future. Give each gene its uniprot id.
def add_gi(feature,gi_to_go):
    if feature.qualifiers['db_xref'][0].split(':')[0] == 'GI':
        gi = feature.qualifiers['db_xref'][0].split(':')[1]
        go = []
        print "Finding uniprot of " + str(gi)
        tmpUnirprot = uni.map(gi,f='P_GI', t='ACC')
        print "TMP unit prot " + str(tmpUnirprot)
        tmpUnirprot = tmpUnirprot[gi].pop()
        print "AFTER POP " + str(tmpUnirprot)
        full_gene = uni.retrieve(tmpUnirprot).split('\n')
        print "FULL GENE " + str(full_gene)
        for line in full_gene:
            if len(line.split("   ")) > 1:
                print line
                if line.split("   ")[1].split(";")[0]=='GO':
                    print line
                    go_att = {}
                    go_att['go_number'] = line.split("   ")[1].split(";")[1].split(':')[1]
                    go_att['attributes'] = line.split("   ")[1].split(";")[2:]
                    go.append(go_att)
        # print dir(feature.location)
        gi_to_go[gi] = go
        print 'gi to go ' + str(gi_to_go)
    else:
        print "NO GI"
        noGi += 1
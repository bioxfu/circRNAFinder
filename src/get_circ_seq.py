#! /usr/bin/env python

import sys
import re

id2seq = {}
circ = {}

for line in sys.stdin:
    ID, seq = line.strip().split('\t')
    seq = seq.upper()
    id2seq[ID] = seq
    ID = re.sub('\|\d+$', '', ID)
    circ[ID] = 1

for ID in sorted(circ):
    ID_list = ID.split('|')
    strand = ID_list[3]
    num = int(ID_list[-1]) + 1
    seq = ''
    
    for i in range(num):
        if strand == '+':
            seq = seq + id2seq['%s|%s' % (ID, i)]
        else:
            seq = id2seq['%s|%s' % (ID, i)] + seq
    
    print ">%s\n%s" % ('|'.join(ID_list[:4]), seq)

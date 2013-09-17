#! /usr/bin/env python

import sys

id2lib = {}
id2reads = {}
lib_dict = {}

hit_mat = open(sys.argv[1], 'w')
bed_file = open(sys.argv[2], 'w')

for line in sys.stdin:
    ID, lib , reads = line.strip().split('\t')
    lib_dict[lib] = 1
    id2lib[ID] = lib + ',' + id2lib.get(ID, '')
    id2reads[ID] = reads + ';' + id2reads.get(ID, '')
    
lib_list = lib_dict.keys()
hit_mat.write('id\t%s\n' % ('\t'.join(lib_list)))

for i in sorted(id2lib):
    
    chrom, start, end, strand = i.split('|')
    start = int(start)
    end = int(end)
    start -= 1
    genomic_len = end - start
    reads = id2reads[i].strip(';').split(';')
    count = len(set(reads))
    reads = ';'.join(set(reads))
    libs = id2lib[i].strip(',').split(',')

    if count >= 2:
        bed_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' 
                       % (chrom, start, end, i, genomic_len, strand, ';'.join(libs), count, reads))
        hit = []
        
        for j in lib_list:
            if j in libs:
                hit.append('1')
            else:
                hit.append('0')
        hit_mat.write('%s\t%s\n' % (i, '\t'.join(hit)))

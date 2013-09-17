#! /usr/bin/env python

import sys

input_file = open(sys.argv[1], 'r')
output_anno = open(sys.argv[2], 'w')
output_bed = open(sys.argv[3], 'w')

id2circ = {}
id2sense_starts = {}
id2sense_ends = {}
id2sense_genes = {}
id2sense_biotype = {}
id2antisen_starts = {}
id2antisen_ends = {}
id2antisen_genes = {}
id2antisen_biotype = {}

for line in input_file:
    line_list = line.strip().split('\t')
    chrom, circ_start, circ_end, ID, genomic_len, circ_str, lib, count, reads = line_list[:9]
    cds_start, cds_end, cds_str, exon_start, exon_end, ref, name, biotype = line_list[12:20]
    ## all start coordinates are 0-based, all end coordinates are 1-based ##
    id2circ[ID] = '\t'.join(line_list[:9])
    ## circRNA overlap with sense/antisense gene ##
    if circ_str == cds_str:
        id2sense_genes[ID] = name + ',' + id2sense_genes.get(ID,'')
        id2sense_biotype[ID] = biotype + ',' + id2sense_biotype.get(ID,'')
        id2sense_starts[ID] = exon_start + ',' + id2sense_starts.get(ID,'')
        id2sense_ends[ID] = exon_end + ',' + id2sense_ends.get(ID,'')
    else:
        id2antisen_genes[ID] = name + ',' + id2antisen_genes.get(ID,'')
        id2antisen_biotype[ID] = biotype + ',' + id2antisen_biotype.get(ID,'')
        id2antisen_starts[ID] = exon_start + ',' + id2antisen_starts.get(ID,'')
        id2antisen_ends[ID] = exon_end + ',' + id2antisen_ends.get(ID,'')
input_file.close()

for ID in sorted(id2circ):
    chrom, circ_start, circ_end, ID, genomic_len, circ_str, lib, count, reads = id2circ[ID].split('\t')
    circ_exons = []
    sense = 'sense'
    locat = 'exon'
    splice = 0
    genes, biotype, starts, ends = '','','',''
    if ID in id2sense_genes:
        genes = id2sense_genes[ID].strip(',')
        biotype = id2sense_biotype[ID].strip(',')
        starts = id2sense_starts[ID].strip(',')
        ends = id2sense_ends[ID].strip(',')
    else:
        genes = id2antisen_genes[ID].strip(',')
        biotype = id2antisen_biotype[ID].strip(',')
        starts = id2antisen_starts[ID].strip(',')
        ends = id2antisen_ends[ID].strip(',')
        sense = 'antisen'
    
    if sense == 'antisen':
        circ_exons.append('%s|%s' % (circ_start, circ_end))
    else:
        exon_starts = starts.split(',')
        exon_ends = ends.split(',')
        for i in range(len(exon_starts)):
            if exon_starts[i] >= circ_end or exon_ends[i] <= circ_start:
                continue
            circ_exons.append('%s|%s' % (exon_starts[i], exon_ends[i]))
    circ_exons.sort()
    
    if len(circ_exons) == 0:
        circ_exons.append('%s|%s' % (circ_start, circ_end))
        locat = 'intron'
    
    temp_list = circ_exons[0].split('|')
    if circ_start == temp_list[0]:
        splice += 1
    else:
        circ_exons[0] = '%s|%s' % (circ_start, temp_list[1])
    
    temp_list = circ_exons[-1].split('|')
    if circ_end == temp_list[1]:
        splice += 1
    else:
        circ_exons[-1] = '%s|%s' % (temp_list[0], circ_end)
    
    if locat == 'intron':
        splice = 0
    
    gene_num = len(genes.split(','))
    output_anno.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ID, genes, gene_num, biotype
                                                       , splice, sense, locat))
    exon_num = len(circ_exons)
    for i in range(exon_num):
        temp = circ_exons[i].split('|')
        output_bed.write('%s\t%s\t%s\t%s|%s|%s\t.\t%s\n' % (chrom, temp[0], temp[1], ID, exon_num-1, i, circ_str))   

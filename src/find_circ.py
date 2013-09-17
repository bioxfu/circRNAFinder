#! /usr/bin/env python

'find_circ.py -- find the circRNA'

import sys, os

if len(sys.argv) == 4:
    in_file = open(sys.argv[1],'r')
    bed_file = open(sys.argv[1]+'.bed','w')
    genome_fa = sys.argv[2]
    out_file = open(sys.argv[3],'w')
else:
    #print("usage: ./find_circ.py input_bwt output")
    #print("or python find_circ.py input_bwt output")
    sys.exit(1)

reads_dict = {'A': {'chrom': {}, 'strand': {}, 'start': {}, 'seq': {}}, 
              'B': {'chrom': {}, 'strand': {}, 'start': {}, 'seq': {}}}

for line in in_file:
    field = line.strip('\n').split('\t')
    ID, strand, chrom, start, subseq = field[0:5]
    read, part = '.'.join(ID.split('.')[:-1]), ID.split('.')[-1]
    # start: 0-based offset into the forward reference strand where leftmost character of the alignment occurs
    # subseq: Read sequence (reverse-complemented if orientation is `-`).
    
    if part == 'A':
        reads_dict['A']['chrom'][read] = chrom + ',' + reads_dict['A']['chrom'].get(read, '')
        reads_dict['A']['strand'][read] = strand + ',' + reads_dict['A']['strand'].get(read, '')
        reads_dict['A']['start'][read] = start + ',' + reads_dict['A']['start'].get(read, '')
        reads_dict['A']['seq'][read] = subseq + ',' + reads_dict['A']['seq'].get(read, '')
    if part == 'B':
        reads_dict['B']['chrom'][read] = chrom + ',' + reads_dict['B']['chrom'].get(read, '')
        reads_dict['B']['strand'][read] = strand + ',' + reads_dict['B']['strand'].get(read, '')
        reads_dict['B']['start'][read] = start + ',' + reads_dict['B']['start'].get(read, '')
        reads_dict['B']['seq'][read] = subseq + ',' + reads_dict['B']['seq'].get(read, '')
                
circ = {}

for x in reads_dict['A']['chrom']:
    if x in reads_dict['B']['chrom']:
        Achrom = reads_dict['A']['chrom'][x][:-1].split(',')
        Astrand = reads_dict['A']['strand'][x][:-1].split(',')
        Astart = [int(n) for n in reads_dict['A']['start'][x][:-1].split(',')]
        Aseq = reads_dict['A']['seq'][x][:-1].split(',')
        Bchrom = reads_dict['B']['chrom'][x][:-1].split(',')
        Bstrand = reads_dict['B']['strand'][x][:-1].split(',')
        Bstart = [int(n) for n in reads_dict['B']['start'][x][:-1].split(',')]
        Bseq = reads_dict['B']['seq'][x][:-1].split(',')
        
        if len(Achrom) == 1:
            for i in range(len(Bchrom)): 
                if Achrom[0] == Bchrom[i] and Astrand[0] == Bstrand[i]:
                    seqA_3p, seqB_5p = 0, 0
                    if Astrand[0] == '+':
                        seqA_3p = Astart[0] + len(Aseq[0])  # 1-based
                        seqB_5p = Bstart[i] + 1  #1-based
                        if seqA_3p > seqB_5p:
                            rec = '%s\t%s\t%s' % (Achrom[0], seqB_5p, seqA_3p)
                            circ[rec] = x + ';' + circ.get(rec, '')
                    else:
                        seqA_3p = Astart[0] + 1  # 1-based
                        seqB_5p = Bstart[i] + len(Bseq[i])  #1-based
                        if seqA_3p < seqB_5p:
                            rec = '%s\t%s\t%s' % (Achrom[0], seqA_3p, seqB_5p)
                            circ[rec] = x + ';' + circ.get(rec, '')
        elif len(Bchrom) == 1:
            for i in range(len(Achrom)): 
                if Achrom[i] == Bchrom[0] and Astrand[i] == Bstrand[0]:
                    seqA_3p, seqB_5p = 0, 0
                    if Bstrand[0] == '+':
                        seqA_3p = Astart[i] + len(Aseq[i])  # 1-based
                        seqB_5p = Bstart[0] + 1  #1-based
                        if seqA_3p > seqB_5p:
                            rec = '%s\t%s\t%s' % (Bchrom[0], seqB_5p, seqA_3p)
                            circ[rec] = x + ';' + circ.get(rec, '')
                    else:
                        seqA_3p = Astart[i] + 1  # 1-based
                        seqB_5p = Bstart[0] + len(Bseq[0])  #1-based
                        if seqA_3p < seqB_5p:
                            rec = '%s\t%s\t%s' % (Bchrom[0], seqA_3p, seqB_5p)
                            circ[rec] = x + ';' + circ.get(rec, '')

for i in circ:
    reads = circ[i][:-1]
    field = i.split('\t')
    ID = '|'.join(field)
    if int(field[2]) - int(field[1]) < 100000 and int(field[1])-3 > 0:
        bed_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (field[0], int(field[1])-3, int(field[2])+2, ID+'|+|'+reads, '.', '+'))
        bed_file.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (field[0], int(field[1])-3, int(field[2])+2, ID+'|-|'+reads, '.', '-'))

in_file.close()
bed_file.close()

bed_cmd = 'bedtools getfasta -name -tab -s -fi %s -bed %s -fo %s.seq' % (genome_fa, sys.argv[1]+'.bed', sys.argv[1]+'.bed') 
os.system(bed_cmd)

in_file = open(sys.argv[1]+'.bed.seq', 'r')
for line in in_file:
    field = line.strip().split('\t')
    seq = field[1].upper()
    if seq.startswith('AG') and seq.endswith('GT'):
        chrom, start, end, strand, reads = field[0].split('|')
        counts = len(reads.split(';'))
        ID = chrom + '|' + start + '|' + end
        out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom, start, end, ID, counts, strand, reads))
in_file.close()
out_file.close()

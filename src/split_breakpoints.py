#! /usr/bin/env python

'''split_breakpoints.py -- filter the anchors alignments and split the read at 
all sites between two anchors'''

import sys

if len(sys.argv) == 3:
    in_file = open(sys.argv[1],'r')
    out_file = open(sys.argv[2], 'w')
else:
    print("usage: ./split_breakpoints.py input_bwt output_fa")
    print("or python split_breakpoints.py input_bwt output_fa")
    sys.exit(1)

reads_dict = {'A': {'chrom': {}, 'strand': {}, 'start': {}}, 
              'B': {'chrom': {}, 'strand': {}, 'start': {}}}

for line in in_file:
    field = line.strip('\n').split('\t')
    ID, strand, chrom, start = field[:4]
    read, part = ID[:-2], ID[-1]
    # start: 0-based offset into the forward reference strand where leftmost character of the alignment occurs
    # subseq: Read sequence (reverse-complemented if orientation is `-`).
    
    if part == 'A':
        reads_dict['A']['chrom'][read] = chrom + ',' + reads_dict['A']['chrom'].get(read, '')
        reads_dict['A']['strand'][read] = strand + ',' + reads_dict['A']['strand'].get(read, '')
        reads_dict['A']['start'][read] = start + ',' + reads_dict['A']['start'].get(read, '')
    if part == 'B':
        reads_dict['B']['chrom'][read] = chrom + ',' + reads_dict['B']['chrom'].get(read, '')
        reads_dict['B']['strand'][read] = strand + ',' + reads_dict['B']['strand'].get(read, '')
        reads_dict['B']['start'][read] = start + ',' + reads_dict['B']['start'].get(read, '')

for x in reads_dict['A']['chrom']:
    if x in reads_dict['B']['chrom']:
        Achrom = reads_dict['A']['chrom'][x][:-1].split(',')
        Astrand = reads_dict['A']['strand'][x][:-1].split(',')
        Astart = [int(n) for n in reads_dict['A']['start'][x][:-1].split(',')]
        Bchrom = reads_dict['B']['chrom'][x][:-1].split(',')
        Bstrand = reads_dict['B']['strand'][x][:-1].split(',')
        Bstart = [int(n) for n in reads_dict['B']['start'][x][:-1].split(',')]
        
        if len(Achrom) == 1:
            temp1 = []
            temp2 = []
            for i in range(len(Bchrom)): 
                if Achrom[0] == Bchrom[i] and Astrand[0] == Bstrand[i]:
                    temp1.append(Bstart[i])
                    if (Astrand[0] == '+' and Astart[0] > Bstart[i]) or (Astrand[0] == '-' and Astart[0] < Bstart[i]):
                        temp2.append(Bstart[i])
            if len(temp1) > 0 and len(temp1) == len(temp2):
                for n in range(17, (len(x) - 17 + 1)): # n is the length of the left part
                    #left = s[:n] + 'GT'
                    #right = 'AG' + s[n:]
                    left = x[:n]
                    right = x[n:]
                    out_file.write('>%s.%d.A\n%s\n' % (x, n, left))
                    out_file.write('>%s.%d.B\n%s\n' % (x, n, right))
                
        elif len(Bchrom) == 1:
            temp1 = []
            temp2 = []
            for i in range(len(Achrom)): 
                if Achrom[i] == Bchrom[0] and Astrand[i] == Bstrand[0]:
                    temp1.append(Astart[i])
                    if (Bstrand[0] == '+' and Astart[i] > Bstart[0]) or (Bstrand[0] == '-' and Astart[i] < Bstart[0]):
                        temp2.append(Astart[i])
            if len(temp1) > 0 and len(temp1) == len(temp2):
                for n in range(17, (len(x) - 17 + 1)): # n is the length of the left part
                    #left = s[:n] + 'GT'
                    #right = 'AG' + s[n:]
                    left = x[:n]
                    right = x[n:]
                    out_file.write('>%s.%d.A\n%s\n' % (x, n, left))
                    out_file.write('>%s.%d.B\n%s\n' % (x, n, right))


in_file.close()
out_file.close()
    

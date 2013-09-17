#! /usr/bin/env python

import sys

intergen_file = open(sys.argv[1], 'r')
genic_file = open(sys.argv[2], 'r')
seq_file = open(sys.argv[3], 'r')
input_file = open(sys.argv[4], 'r')
output_file = open(sys.argv[5], 'w')
gff_file = open(sys.argv[6], 'w')
id_prefix = sys.argv[7]

id2intergen = {}
id2genic = {}
id2seq = {}

for line in intergen_file:
    line_list = line.strip().split('\t')
    id2intergen[line_list[3]] = line_list[-1]
intergen_file.close()

for line in genic_file:
    line_list = line.strip().split('\t')
    id2genic[line_list[0]] = '\t'.join(line_list[1:])
genic_file.close()

for line in seq_file:
    line_list = line.strip().split('\t')
    id2seq[line_list[0][1:]] = line_list[1]
seq_file.close()


output_file.write("id\tchr\tstart\tend\tstr\tcounts\tloci\tgenomic_len\trna_len\tlib_list\t \
lib_num\tgene_list\tgene_num\tbiotype\tsplice_site\tsense\treads\tseq\n")

index = 1;
for line in input_file:
    chrom, start, end, ID, genomic_len, strand, lib_list, counts, reads = line.strip().split('\t')
    lib_num = len(lib_list.split(';'))
    gene_list, gene_num, biotype, splice_site, sense, loci = ('.', '.', '.', '.', '.', 'intergenic')
    if ID in id2genic:
        gene_list, gene_num, biotype, splice_site, sense, loci = id2genic[ID].split('\t')
    ID = id_prefix + '_circRNA%06d' % index    
    start = int(start) + 1;
    seq = id2seq['%s|%s|%s|%s'% (chrom, start, end, strand)]
    rna_len = len(seq) - 1; ## remove the '='
    output_file.write('\t'.join([ID, chrom, str(start), end, strand, counts, loci, genomic_len, str(rna_len), 
                                 lib_list, str(lib_num), gene_list, gene_num, biotype, splice_site, sense,
                                 reads, seq]) + '\n')
    gff_file.write('%s\tcircRNAtools\tcircRNA\t%s\t%s\t.\t%s\t.\tID=%s\n' % (chrom, start, end, strand, ID))
    index += 1

#! /usr/bin/env python

import argparse, re, os
parser = argparse.ArgumentParser(description = 'Download genome/transcript sequences and gene annotations')
parser.add_argument('species', choices=['hg19','mm10','TAIR10'], help='choose a species (Human, Mouse, Arabidopsis)')
parser.add_argument('-d', '--download', action='store_true', help='download sequences or annotations')
parser.add_argument('-b', '--build', action='store_true', help='build sequences or annotations')
parser.add_argument('-g', '--genome', action='store_true', help='download or build genome sequences')
parser.add_argument('-t', '--transcriptome', action='store_true', help='download or build transcriptome sequences')
parser.add_argument('-a', '--annotation', action='store_true', help='download or build gene annotations')

args = parser.parse_args()

genome_dict = {'hg19': 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz', 
'mm10': 'ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz',
'TAIR10': 'ftp://ftp.ensemblgenomes.org/pub/current/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.*.dna.toplevel.fa.gz'}

trans_dict = {'hg19': ['ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.*.cdna.abinitio.fa.gz',
						'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.*.cdna.all.fa.gz',
						'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/ncrna/Homo_sapiens.*.ncrna.fa.gz'],
'mm10': ['ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/cdna/Mus_musculus.*.cdna.abinitio.fa.gz',
		'ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/cdna/Mus_musculus.*.cdna.all.fa.gz',
		'ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/ncrna/Mus_musculus.*.ncrna.fa.gz'],
'TAIR10': ['ftp://ftp.ensemblgenomes.org/pub/current/plants/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.*.cdna.abinitio.fa.gz',
			'ftp://ftp.ensemblgenomes.org/pub/current/plants/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.*.cdna.all.fa.gz',
			'ftp://ftp.ensemblgenomes.org/pub/current/plants/fasta/arabidopsis_thaliana/ncrna/Arabidopsis_thaliana.*.ncrna.fa.gz']}
anno_dict = {'hg19': 'ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/*.gtf.gz', 
'mm10': 'ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/*.gtf.gz',
'TAIR10': 'ftp://ftp.ensemblgenomes.org/pub/current/plants//gtf/arabidopsis_thaliana/*.gtf.gz'}

def gtf_build(gtf, build):
	input_file = open(gtf,'r')
	output_file = open(build,'w')
	tx2gene = {}
	tx2exon_starts = {}
	tx2exon_ends = {}
	tx2cds_starts = {}
	tx2cds_ends = {}
	
	for line in input_file:
		if line.startswith('#'):
			continue
		line_list = line.strip().split('\t')
		chrom, biotype, feature, start, end, strand, ID = (line_list[0],line_list[1],line_list[2],line_list[3],line_list[4],line_list[6],line_list[8])
		if gtf == 'hg19.gtf' or gtf == 'mm10.gtf':
			chrom = 'chr' + chrom		
		start = str(int(start) - 1) ## 0-based
		if re.search('gene_id \"(.+?)\".+transcript_id \"(.+?)\"', ID) is not None:
			gene_id, tx_id = re.search('gene_id \"(.+?)\".+transcript_id \"(.+?)\"', ID).groups()
			tx2gene[tx_id] = '%s|%s|%s|%s' % (chrom, strand, gene_id, biotype)
			if feature == 'exon':
				tx2exon_starts[tx_id] = start + ',' + tx2exon_starts.get(tx_id, '')
				tx2exon_ends[tx_id] = end + ',' + tx2exon_ends.get(tx_id, '')
			if feature == 'CDS':
				tx2cds_starts[tx_id] = start + ',' + tx2cds_starts.get(tx_id, '')
				tx2cds_ends[tx_id] = end + ',' + tx2cds_ends.get(tx_id, '')
	
	gene2repretx = {} ## representative transcript (repretx) is the longest transcript for each gene 
	trans2len = {}   
	for tx_id in tx2gene:
		chrom, strand, gene_id, biotype = tx2gene[tx_id].split('|')
		exon_starts = sorted([int(i) for i in tx2exon_starts[tx_id].strip(',').split(',')])
		exon_ends = sorted([int(i) for i in tx2exon_ends[tx_id].strip(',').split(',')])
		tx_len = 0
		for i in range(len(exon_starts)):
			tx_len += (exon_ends[i] - exon_starts[i])
		trans2len[tx_id] = tx_len
		if gene_id in gene2repretx:
			if tx_len > trans2len[gene2repretx[gene_id]]:
				gene2repretx[gene_id] = tx_id
		else:
			gene2repretx[gene_id] = tx_id
	
	for tx_id in sorted(tx2gene):
		chrom, strand, gene_id, biotype = tx2gene[tx_id].split('|')
		if tx_id == gene2repretx[gene_id]:
			exon_starts = [str(j) for j in sorted([int(i) for i in tx2exon_starts[tx_id].strip(',').split(',')])]
			exon_ends = [str(j) for j in sorted([int(i) for i in tx2exon_ends[tx_id].strip(',').split(',')])]
			tx_start = exon_starts[0]
			tx_end = exon_ends[-1]
		
			cds_start = '.'
			cds_end = '.'
		
			if tx_id in tx2cds_starts:
				cds_starts = [str(j) for j in sorted([int(i) for i in tx2cds_starts[tx_id].strip(',').split(',')])]
				cds_ends = [str(j) for j in sorted([int(i) for i in tx2cds_ends[tx_id].strip(',').split(',')])]
				cds_start = cds_starts[0]
				cds_end = cds_ends[-1]
		
			output_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom, tx_start, tx_end, cds_start, cds_end, strand, ','.join(exon_starts), ','.join(exon_ends), tx_id, gene_id, biotype))


if args.download:
	if args.genome:
		print '[download %s genome]' % args.species
		if args.species == 'hg19' or args.species == 'mm10':
			print 'wget -q %s -O %s.tar.gz' % (genome_dict[args.species], args.species)
			os.system('wget -q %s -O %s.tar.gz' % (genome_dict[args.species], args.species))
			print 'tar -zxf %s.tar.gz' % args.species
			os.system('tar -zxf %s.tar.gz' % args.species)
			print 'cat chr*.fa > %s_dna.fa' % args.species
			os.system('cat chr*.fa > %s_dna.fa' % args.species)
			print 'rm chr*.fa'
			os.system('rm chr*.fa')
		else:
			print 'wget -q %s -O %s.fa.gz' % (genome_dict[args.species], args.species)
			os.system('wget -q %s -O %s.fa.gz' % (genome_dict[args.species], args.species))
			print 'zcat %s.fa.gz > %s_dna.fa' % (args.species, args.species)
			os.system('zcat %s.fa.gz > %s_dna.fa' % (args.species, args.species))
			print 'rm %s.fa.gz' % args.species
			os.system('rm %s.fa.gz' % args.species)

	elif args.transcriptome:
		print '[download %s transcriptome]' % args.species
		for i in trans_dict[args.species]:
			print 'wget -q %s' % i
			os.system('wget -q %s' % i)
		print 'zcat *.fa.gz > %s_trans.fa' % args.species
		os.system('zcat *.fa.gz > %s_trans.fa' % args.species)
		print 'rm *.fa.gz'
		os.system('rm *.fa.gz')
	elif args.annotation:
		print '[download %s gene annotation]' % args.species
		print 'wget -q %s -O %s.gtf.gz' % (anno_dict[args.species], args.species)
		os.system('wget -q %s -O %s.gtf.gz' % (anno_dict[args.species], args.species))
		print 'gzip -d %s.gtf.gz' % args.species
		os.system('gzip -d %s.gtf.gz' % args.species)
	else:
		print 'please specify -g/--genome or -t/--transcriptome or -a/--annotation'
elif args.build:
	if args.genome:
		print '[build %s genome]' % args.species
		print 'bowtie-build %s_dna.fa %s_dna' % (args.species, args.species)
		os.system('bowtie-build %s_dna.fa %s_dna' % (args.species, args.species))
	elif args.transcriptome:
		print '[build %s transcriptome]' % args.species
		print 'bowtie-build %s_trans.fa %s_trans' % (args.species, args.species)
		os.system('bowtie-build %s_trans.fa %s_trans' % (args.species, args.species))
	elif args.annotation:
		print '[build %s gene annotation]' % args.species
		print 'gtf_build(%s.gtf, %s.gtf.build)' % (args.species, args.species)
		gtf_build(args.species+'.gtf', args.species+'.gtf.build')
	else:
		print 'please specify -g/--genome or -t/--transcriptome or -a/--annotation'
else:
	print 'please specify -d/--download or -b/--build'



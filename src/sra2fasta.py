#! /usr/bin/env python

import os, argparse
parser = argparse.ArgumentParser(description = 'Convert .sra to .fasta file')
group = parser.add_mutually_exclusive_group()
parser.add_argument('input', help='.sra file')
group.add_argument('-s', '--single-end', dest='single', action='store_true', help='')
group.add_argument('-p', '--paired-end', dest='paired', action='store_true', help='')
args = parser.parse_args()

if args.single:
	cmd = 'fastq-dump --split-spot --skip-technical --gzip -I --fasta 0 %s' % args.input
	print cmd
	os.system(cmd) 
elif args.paired:
	cmd = 'fastq-dump --split-spot --skip-technical --split-files --gzip -I --fasta 0 %s' % args.input
	print cmd
	os.system(cmd)
else:
	print 'please specify -s/--single-end or -p/--paired-end'


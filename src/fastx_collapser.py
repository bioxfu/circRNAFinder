#! /usr/bin/env python

import os, argparse, re
parser = argparse.ArgumentParser(description = '')
group = parser.add_mutually_exclusive_group()
parser.add_argument('input', help='file1.gz,file2.gz,...')
parser.add_argument('output', help='.fasta file')
group.add_argument('-f', '--fasta', action='store_true', help='input is fasta file')
group.add_argument('-q', '--fastq', action='store_true', help='input is fastq file')
args = parser.parse_args()

args.input = re.sub(',',' ',args.input)

if args.fasta:
	cmd = 'zcat %s | fastx_collapser -o %s' % (args.input, args.output)
	print cmd
	os.system(cmd) 
elif args.fastq:
	cmd = 'zcat %s | fastq_to_fasta | fastx_collapser -o %s' % (args.input, args.output)
	print cmd
	os.system(cmd)
else:
	print 'please specify -f/--fasta or -q/--fastq'


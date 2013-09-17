#! /usr/bin/env python

"circFind.py -- "

import sys, os, ConfigParser, argparse
parser = argparse.ArgumentParser(description='circFind is a pipeline to find circular RNA')
parser.add_argument('config_file', help='')
parser.add_argument('-s', '--step', choices=range(8), type=int, help='0:all steps, 1:ori_reads_genome_align, 2:ori_reads_trans_align, 3:split_anchors, 4:anchors_align, 5:split_breakpoints, 6:breakpoints_align, 7:find_circrna')

parser.add_argument('-r', '--remove-temp', dest='rm', action='store_true', help='')
args = parser.parse_args()

#=========================================================
def ori_reads_genome_align(output_dir, genome_idx, fa_file, run):
    print '=== ori_reads_genome_align ==='
    output = {'genome_unalign': output_dir + 'genome_unalign.fa', 'genome_align': output_dir + 'genome_align.bwt'}
    cmd = 'bowtie -f -v3 -p4 --un %s %s %s %s' % (output['genome_unalign'], genome_idx, fa_file, output['genome_align'])
    if run == True:
        print cmd
        os.system(cmd)
    return(output)

def ori_reads_trans_align(output_dir, trans_idx, fa_file, run):
    print '=== ori_reads_trans_align ==='
    output = {'trans_unalign': output_dir + 'trans_unalign.fa', 'trans_align': output_dir + 'trans_align.bwt'}
    cmd = 'bowtie -f -v3 -p4 --un %s %s %s %s' % (output['trans_unalign'], trans_idx, fa_file, output['trans_align'])
    if run == True:
        print cmd
        os.system(cmd)
    return(output)

def split_anchors(output_dir, fa_file, run):
    print '=== split anchors ==='
    output = {'anchors': output_dir + 'anchors.fa'}
    cmd = 'split_anchors.py %s %s' % (fa_file, output['anchors'])
    if run == True:
        print cmd
        os.system(cmd)
    return(output)

def anchors_align(output_dir, genome_idx, fa_file, run):
    print '=== anchors align ==='
    output = {'anchors_align': output_dir + 'anchors.align.bwt'}
    cmd = 'bowtie -f -v3 --best --strata -k11 -m10 -p4 %s %s %s' % (genome_idx, fa_file, output['anchors_align'])
    if run == True:
        print cmd
        os.system(cmd)
    return(output)

def split_breakpoints(output_dir, bwt, run):
    print '=== split reads at breakpoints ==='
    output = {'breakpoints': output_dir + 'breakpoints.fa'}
    cmd = 'split_breakpoints.py %s %s' % (bwt, output['breakpoints'])
    if run == True:
        print cmd
        os.system(cmd)
    return(output)

def breakpoints_align(output_dir, genome_idx, fa_file, run):
    print '=== breakpoints mapping ==='
    output = {'breakpoints_align': output_dir + 'breakpoints.align.bwt'}
    cmd = 'bowtie -f -v3 --best --strata -k11 -m10 -p4 %s %s %s' % (genome_idx, fa_file, output['breakpoints_align'])
    if run == True:
        print cmd 
        os.system(cmd)
    return(output)

def find_circrna(output_dir, sample_name, genome_fa, bwt, run):
    print "=== find circRNAs ==="    
    output = {'circ': output_dir + sample_name + '.circ.txt'}
    
    cmds = ['find_circ.py %s %s %s.circ' % (bwt, genome_fa, bwt),
            'cat %s.circ|sort -k1,1 -k2,2n|uniq > %s' % (bwt, output['circ']),
            ]
    
    for cmd in cmds:
        if run == True:
            print cmd
            os.system(cmd)

#================================================================
cf = ConfigParser.ConfigParser()  
cf.read(args.config_file)  

for sec in cf.sections():
    sample_name = cf.get(sec, 'sample_name')
    reads_file = cf.get(sec, 'reads_file')
    genome_index = cf.get(sec, 'genome_index')
    genome_seq = cf.get(sec, 'genome_seq')
    trans_index = cf.get(sec, 'trans_index')
    
    temp_dir = './' + sample_name + '/temp/'
    output_dir = './' + sample_name + '/output/'
    
    os.system('mkdir -p ' + temp_dir)
    os.system('mkdir -p ' + output_dir)
    
    tag = [False, False, False, False, False, False, False, False, False]
    if args.step > 0:
        tag[args.step] = True
    elif args.step == 0:
        tag = [True, True, True, True, True, True, True, True, True]
    else:
        print 'please specify the -s/--step'
        sys.exit()
    
    ori_reads_genome_align_out = ori_reads_genome_align(temp_dir, genome_index, reads_file, tag[1])

    ori_reads_trans_align_out = ori_reads_trans_align(temp_dir, trans_index, ori_reads_genome_align_out['genome_unalign'], tag[2])
        
    split_anchors_out = split_anchors(temp_dir, ori_reads_trans_align_out['trans_unalign'], tag[3])
     
    anchors_align_out = anchors_align(temp_dir, genome_index, split_anchors_out['anchors'], tag[4])
       
    split_breakpoints_out = split_breakpoints(temp_dir, anchors_align_out['anchors_align'], tag[5])
       
    breakpoints_align_out = breakpoints_align(temp_dir, genome_index, split_breakpoints_out['breakpoints'], tag[6])
   
    find_circrna_out = find_circrna(output_dir, sample_name, genome_seq, breakpoints_align_out['breakpoints_align'], tag[7])
   
    if args.rm == True:
        print "=== clean the temp directory ==="
        os.system('rm -r %s' % temp_dir)


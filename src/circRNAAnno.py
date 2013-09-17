#! /usr/bin/env python

"circAnno.py -- "

import sys, os, ConfigParser, argparse
parser = argparse.ArgumentParser(description='circAnno is a pipeline to annotate circular RNA')
parser.add_argument('config_file', help='')
parser.add_argument('-s', '--step', choices=range(5), type=int, help='0:all steps, 1:combine libraries, 2:intergenic, 3:genic, 4:maketable')

parser.add_argument('-r', '--remove-temp', dest='rm', action='store_true', help='')
args = parser.parse_args()

#================================================================
def combine_libs(output_dir, lib_name_list, lib_path_list, run):
    print '=== combine libraries==='
    output = {'id2lib': output_dir + '/temp/id2lib',
              'combLib_hit_mat': output_dir + '/temp/combLib.hit.mat',
              'combLib_bed': output_dir + '/temp/combLib.bed'}
    
    if os.path.exists(output['id2lib']):
        os.system('rm ' + output['id2lib'])
    
    for i in (range(len(lib_name_list))):
        cmd = '''awk '{print $4"|"$6"\\t%s\\t"$7}' %s >> %s''' % (lib_name_list[i], lib_path_list[i], output['id2lib'])
        if run == True:
            print cmd
            os.system(cmd)
    
    cmds = ['cat %s|combLibs.py %s %s' % (output['id2lib'], output['combLib_hit_mat'], output['combLib_bed']),
            'sort -k1,1 -k2,2n %s > %s.temp' % (output['combLib_bed'], output['combLib_bed']),
            'mv %s.temp %s' % (output['combLib_bed'], output['combLib_bed']),
            'compLib.r %s %s' % (output['combLib_hit_mat'], output_dir + '/figure/')]
    
    for cmd in cmds:
        if run == True:
            print cmd
            os.system(cmd)
    
    return(output)

def intergen_circ_anno(output_dir, combLib_bed, gene_bed, genome_seq, run):
    print '=== get intergenic circRNA anno and sequence ==='
    output = {'intergen_bed': output_dir + '/temp/combLib.intergen.bed',
              'intergen_fa': output_dir + '/temp/combLib.intergen.fa'}

    cmds = ['''bedtools intersect -a %s -b %s -wa -wb -v|sort -k1,1 -k2,2n|awk '{print $0"\\tintergen"}'> %s''' % (combLib_bed, gene_bed, output['intergen_bed']),
            '''bedtools getfasta -s -name -fi %s -bed %s -fo %s''' % (genome_seq, output['intergen_bed'], output['intergen_fa'])]
    
    for cmd in cmds:
        if run == True:
            print cmd
            os.system(cmd)
    
    return(output)

def genic_circ_anno(output_dir, combLib_bed, gene_bed, genome_seq, run):
    print '=== get genic circRNA anno and sequence ==='
    output = {'genic_bed': output_dir + '/temp/combLib.genic.bed',
              'genic_anno': output_dir + '/temp/combLib.genic.anno',
              'genic_exon_bed': output_dir + '/temp/combLib.genic.exon.bed',
              'genic_exon_bed_seq': output_dir + '/temp/combLib.genic.exon.bed.seq',
              'genic_fa': output_dir + '/temp/combLib.genic.fa'}
    
    cmds = ['''bedtools intersect -a %s -b %s -wa -wb|sort -k1,1 -k2,2n > %s''' % (combLib_bed, gene_bed, output['genic_bed']),
            '''find_circ_exons.py %s %s %s''' % (output['genic_bed'], output['genic_anno'], output['genic_exon_bed']),
            '''sort -t '|' -k1,1 -k2,2n %s > %s.temp''' % (output['genic_anno'], output['genic_anno']),
            '''mv %s.temp %s''' % (output['genic_anno'], output['genic_anno']),
            '''bedtools getfasta -tab -s -name -fi %s -bed %s -fo %s''' % (genome_seq, output['genic_exon_bed'], output['genic_exon_bed_seq']),
            '''cat %s|get_circ_seq.py > %s''' % (output['genic_exon_bed_seq'], output['genic_fa'])]
    
    for cmd in cmds:
        if run == True:
            print cmd
            os.system(cmd)
    
    return(output)


def make_annoTab(output_dir, intergen_fa, genic_fa, intergen_bed, genic_anno, combLib_bed, id_prefix, run):
    print '=== make annoTable ==='
    output = {'all_fa': output_dir + '/temp/combLib.all.fa',
              'all_cut': output_dir + '/temp/combLib.all.cut.txt',
              'annoTab': output_dir + '/table/combLib_annoTab.txt',
              'gff': output_dir + '/table/combLib_annoTab.gff'}

    cmds = ['cat %s %s > %s' % (intergen_fa, genic_fa, output['all_fa']),
            '''cat %s|split_circ_seq.py |sort -t '|' -k1,1 -k2,2n > %s''' % (output['all_fa'], output['all_cut']),
            'make_annoTable.py %s %s %s %s %s %s %s' % (intergen_bed, genic_anno, output['all_cut'], combLib_bed, output['annoTab'], output['gff'], id_prefix)]

    for cmd in cmds:
        if run == True:
            print cmd
            os.system(cmd)
    
    return(output)

# def mirna_target(dir, intergen_fa, genic_fa, mirna_fa, run):
#     print '=== search miRNA target ==='
#     output = {'all_fa': dir + '/temp/combLib.all.fa',
#               'mirna_bwt': dir + '/temp/combLib.all.mirna.bwt',
#               'mirna_hit': dir + '/temp/combLib.all.mirna.hit',
#               'genic_exon_bed_seq': dir + '/temp/combLib.genic.exon.bed.seq',
#               'genic_fa': dir + '/temp/combLib.genic.fa'}
#     
#     cmds = ['cat %s %s > %s' % (intergen_fa, genic_fa, output['all_fa']),
#             'bowtie-build %s %s' % (output['all_fa'], output['all_fa']),
#             'bowtie -f -v3 -a -p4 %s %s %s' % (output['all_fa'], mirna_fa, output['mirna_bwt']),
#             '''cat %s|mirna_match.py|sort -t '|' -k1,1 -k2,2n > %s''' % (output['mirna_bwt'], output['mirna_hit'])]
#     
#     for cmd in cmds:
#         if run == True:
#             print cmd
#             os.system(cmd)
#     
#     return(output)
# def circ_junc_quant(dir, id_prefix, sample_names, libs_size, id2lib, run):
#     print '=== quantify and normalize the circRNA junction reads ==='
#     output = {'circ_junc_quant': dir + '/table/circ_junc.quant',
#               'circ_junc_quant_rpm': dir + '/table/circ_junc.quant.rpm.norm'}
#     
#     cmds = ['''echo \#id,loci,%s|tr ',' '\t' > %s''' % (sample_names, output['circ_junc_quant']),
#             '''circ_junc_quant.py %s %s %s|sort -t '|' -k1,1 -k2,2n|add_ID.py >> %s''' % (id_prefix, sample_names, id2lib, output['circ_junc_quant']),
#             'RPM_normalize.py %s %s > %s' % (libs_size, output['circ_junc_quant'], output['circ_junc_quant_rpm'])]
# 
#     for cmd in cmds:
#         if run == True:
#             print cmd
#             os.system(cmd)
#     
#     return(output)

#================================================================
cf = ConfigParser.ConfigParser()  
cf.read(args.config_file)  

for sec in cf.sections():
    sample_names = cf.get(sec, 'sample_names')
    bowtie_index = cf.get(sec, 'genome_index')
    genome_seq = cf.get(sec, 'genome_seq')
    gene_bed = cf.get(sec, 'gene_bed')
    id_prefix = cf.get(sec, 'id_prefix')
    
    temp_dir = './' + sec + '/temp/'
    figure_dir = './' + sec + '/figure/'
    table_dir = './' + sec + '/table/'
    
    os.system('mkdir -p ' + temp_dir)
    os.system('mkdir -p ' + figure_dir)
    os.system('mkdir -p ' + table_dir)

    lib_name_list = [i.strip() for i in sample_names.split(',')]
    lib_path_list = [i+'/output/'+i+'.circ.txt' for i in lib_name_list]
    
    tag = [False, False, False, False, False, False]
    if args.step > 0:
        tag[args.step] = True
    elif args.step == 0:
        tag = [True, True, True, True, True, True]
    else:
        print 'please specify the -s/--step'
        sys.exit()
        
    combine_libs_out = combine_libs(sec, lib_name_list, lib_path_list, tag[1])

    intergen_circ_anno_out = intergen_circ_anno(sec, combine_libs_out['combLib_bed'], gene_bed, genome_seq, tag[2])
      
    genic_circ_anno_out = genic_circ_anno(sec, combine_libs_out['combLib_bed'], gene_bed, genome_seq, tag[3])
      
    make_annoTab_out = make_annoTab(sec, intergen_circ_anno_out['intergen_fa'], genic_circ_anno_out['genic_fa'], intergen_circ_anno_out['intergen_bed'], genic_circ_anno_out['genic_anno'], combine_libs_out['combLib_bed'], id_prefix, tag[4])
    
    if args.rm == True:
        print "=== clean the temp directory ==="
        os.system('rm -r %s' % temp_dir)

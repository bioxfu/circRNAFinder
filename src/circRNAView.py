#! /usr/bin/env python

"circView.py -- "

import sys
import ConfigParser
import igv
import re

if len(sys.argv) == 2:
    cfg_name =  sys.argv[1]
else:
    print("usage: circView.py config_file")
    sys.exit(1)

#================================================================
cf = ConfigParser.ConfigParser()  
cf.read(cfg_name)  

for sec in cf.sections():
    gff_file = cf.get(sec, 'gff_file')
    png_folder = cf.get(sec, 'png_folder')
    Genome = cf.get(sec, 'Genome')
    
    gff = open(gff_file, 'r')
    gff_track = re.sub('.+/', '', gff_file)
    print 'input: %s\noutput:%s\ntrack:%s\n' % (gff_file, png_folder, gff_track)

    igv = igv.IGV()
    print 'loading genome...'
    igv.genome(Genome)
     
    print 'loading track...'
    igv.load(gff_file)
     
    print 'expanding track...'
    igv.expand('Gene')
    igv.expand(gff_track)
    
    for g in gff:
        ID = re.sub('ID=', '', g.strip().split('\t')[-1])
        png = png_folder + ID + '.png'
        igv.go(ID)
        igv.save(png)
        print 'plot %s into %s' % (ID, png)

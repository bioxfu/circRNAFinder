#! /usr/bin/env python

'''split_anchors.py -- split the unaligned reads to get the 20nt anchors
   left part (20nt at 5' end) of the reads is 'A'
   right part (20nt at 3' end) of the reads is 'B' '''

import sys

if len(sys.argv) == 3:
    in_file = open(sys.argv[1], 'r')
    out_file = open(sys.argv[2], 'w')
else:
    print("usage: ./split_anchors.py input_fa output_fa")
    print("or python split_anchors.py input_fa output_fa")
    sys.exit(1)
################################################################################

for line in in_file:
    line = line.strip()
    if not line.startswith('>') and len(line) >= 40:
        anchorA = line[:20]
        anchorA_id = line + ".A"
        anchorB = line[-20:]
        anchorB_id = line + ".B"
        out_file.write('>%s\n%s\n>%s\n%s\n' % (anchorA_id, anchorA, anchorB_id, anchorB))

in_file.close()
out_file.close()

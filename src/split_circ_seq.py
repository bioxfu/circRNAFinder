#! /usr/bin/env python

import sys

ID_list = []
seq_list = []

for line in sys.stdin:
    if line.startswith('>'):
        ID_list.append(line.strip())
    else:
        seq_list.append(line.strip())

for i in range(len(ID_list)):
    ID = ID_list[i]
    seq = seq_list[i]
    seq_len = len(seq) / 2
    left = seq[:seq_len]
    right = seq[seq_len:]
    new_seq = right + '=' + left
    print('%s\t%s' % (ID, new_seq))

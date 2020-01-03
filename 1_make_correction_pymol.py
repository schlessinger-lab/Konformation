#!/usr/bin/env python3

# load all pdb that need to check for missing residues for correction
import sys

if len(sys.argv) != 2: sys.exit('\n\t> x.py\n\t\t[missing.file]\n')

with open(sys.argv[1], 'r') as fi:
  pdb = [[l.split('|')[0], l.split('|')[1]] for l in fi]

with open(sys.argv[1]+'.pml', 'w') as fo:
  fo.write('load 1atp.pdb\n')
  for p in pdb:
    name = p[0].split('.')[0]
    fo.write('load {0}, {1}.{2}\n'.format(p[0], name, p[1]))
  fo.write('hide everything\nshow cartoon\n')
  fo.write('color white\ncolor red, 1atp\n')
  fo.write('save correction.pse\n')


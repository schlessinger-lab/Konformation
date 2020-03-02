#!/usr/bin/env python3

##########################################################################
#
#  Peter M.U. Ung @ MSSM
#
#  Load all pdb that need to check for missing residues for correction

#  Compile a PyMOL session file with kinase structures that 
#  2_kinase_conf_classifier.py could not get the required structural
#  data (either because of misalignment, failed superposition, missing
#  residues, unrecognizable residue pattern).
#  From the PyMOL, user will inspect the structures, and manually select
#  the residues required for structural parameter generation and then
#  *rerun* 2_kinase_conf_classifier.py to incorporate the missing residues.
#
##########################################################################

import sys,os
msg = '''
  > {0}
    [ List of PDB with missing structural info ]\n
e.g.>  x.py   kinfo_pdb.example.missing.txt\n'''.format(sys.argv[0])
if len(sys.argv) != 2: sys.exit(msg)

with open(sys.argv[1], 'r') as fi:
  pdb = [[l.split('|')[0], l.split('|')[1]] for l in fi]

with open(sys.argv[1]+'.pml', 'w') as fo:
#  fo.write('load 1atp.pdb\n')
  for p in pdb:
    name = p[0].split('.')[0]
    fo.write('load {0}, {1}.{2}\n'.format(p[0], name, p[1]))

  fo.write('hide everything\nshow cartoon\n')
  fo.write('color white\ncolor red, 1atp\n')
  fo.write('save {0}.pse\n'.format(sys.argv[1]))

os.system('pymol -c {0}.pml'.format(sys.argv[1]))

#!/usr/bin/env python3

##########################################################################
##
##	Peter M.U. Ung @ MSSM
##
##	v1.0 - 14.02.10     4_PDB_3D_superimpose.py
##	v2.0 - 16.06.13     Revised the alignment sequence. Updated to only
##                          include the C-lobe catalytic beta-sheets (most
##                          consistent in all typical kinases), exclude
##                          most of the alpha-helices
##	v3.0 - 17.01.09     Rewrite the process into OOP format to allow
##			    check on 'super'. If failed, use 'align' instead
##  v4.0 - 18.03.12     renamed and rebuilt as a function
##  v4.1   20.02.07   search "Executive: RMS" to cover both PyMOL 1.x and 2.x
##
##	Purpose: Read in the template structure and the model structure(s)
##		 and use PyMOL's Superimposition function to align the 
##		 models to the template. Certain residues of the template
##		 can be selected for alignment by specifying them in the 
##		 input.
##
##	e.g:	> x.py templ.pdb '*.model.pdb' pdb ent '10-20+30-120'
##
##      Use the catalytic beta-sheets of C-lobe for typical kinase alignment
##      for both S/T- and Y-kinases. Need to check the results since sometimes
##      the superposition fails in PyMOL (too few atoms)
##          1ATP    'resi 122-138+162-183'
##
####      Use the C-lobe as base for typical kinase alignment (v1.0):
####      For serine/theorine-kinases:    For tyrosine-kinases:
####      1ATP    'resi 125-182+205-305'  2BDF    'resi 343-402+430-507'
##
##########################################################################

import sys

msg = """\n  > {0}
      [PyMOL executable]
      [Template Structure] [Model Structure(s)]
      [Input Model File Extension] [Output Model File Extension]

      -a=[Template residues for alignment: -a=(PyMOL format)]
      -r=[Remove ligands, ions, or water: -r=(PyMOL format)]\n
      #  Using PyMOL function "super" to superimpose all model structures
      #  to the supplied Template structure.\n
      #  Specifying residues for alignment in PyMOL format. 
        e.g. 'resi 10-20+30-120'
      #  Specifying residues for removal in PyMOL format. 
        e.g.  'resname HOH+CL+NA+SO4+PO4+UNX' or 
              'not poly and not org and not resname MG+MN'\n
      e.g. > {0}
        /usr/bin/pymol template.pdb 'model.*.pdb' pdb mod.pdb -a='resi 12'\n
      # For most kinases: 1ATP 'resi 122-183+162-183'\n""".format(sys.argv[0])
#if len(sys.argv) < 5 or len(sys.argv) > 8: sys.exit(msg)
#print("### Input variables ###\n{0}\n".format(sys.argv))

##########################################################################
import re,os

##########################################################################
def SuperposePDB( pymol_exec, templ_pdb, Targets, extIn, extOt, resid, outpref ):
  Targets.sort()
  print('\n\033[34m### Template structure file ###\033[0m\n{0}'.format(templ_pdb))
  print('\n\033[34m### Target structure file(s) ###\033[0m\n{0}\n'.format(len(Targets)))

  templ = templ_pdb.split('/')[-1].split('.pdb')[0]
  AlignStructures(pymol_exec, templ_pdb, resid, outpref, Targets, extIn, extOt, 'super')

  ## Check pymol alignment Log
  Mdls, Atoms = [], []
  # Extract information from pymol-log file
  pymol_pref = '{0}.{1}.{2}'.format(outpref, templ, 'super')
  with open('{0}.pymol-log'.format(pymol_pref), 'r') as fi:
    for line in fi:
      if re.search(r'PyMOL>load', line):
        curr_pdb = line.split('load ')[1].split(', ')
        Mdls.append(curr_pdb)
      if re.search(r'Executive: RMS', line):
        atom = int(line.split('=')[1].split('(')[1].split('to')[0])
        Atoms.append([line.split(':')[1], atom])
      if re.search(r'Executive: Error', line):
        # odd case where no atom left after rms refinement
        Atoms.append([line, 0])
        print('{0} | {1}'.format(curr_pdb, line))

  del Mdls[0]   # remove the first item, the reference structure

  # Write out extracted information and identify bad alignment
  Aligns = []
  with open('{0}.mod-super.log'.format(pymol_pref), 'w') as fo:
    # write to alignment, save models that fail to have 'enough atoms'
    for idx, Names in enumerate(Mdls):
      if Atoms[idx][1] < 60:
        print('  ** Insufficient Number of Atom for Structure Superposition **')
        print('     Atom: {0} -- fewer than threshold 60'.format(Atoms[idx][1]))
        fo.write('\n  #1# Superpose Warning: Check {0}: atom < 60\n'.format(Names[0]))
        Aligns.append(Names[0])
      fo.write('{0} : {1}'.format(Names[1].rstrip(), Atoms[idx][0]))

  ## Redo the alignment of models that failed the 'superimpose' with 'align'
  AlignStructures(pymol_exec, templ_pdb, resid, outpref, Aligns, extIn, extOt, 'align')


##########################################################################

def AlignStructures( pymol_exec, templ_pdb, resid, outpref, Models, 
                      extIn, extOt, align_mode ):

  templ = templ_pdb.split('/')[-1].split('.')[0]
  pymol_pref = '{0}.{1}.{2}'.format(outpref, templ, align_mode)
  m = open(pymol_pref+'.pml','w')

  ## Load the Template PDB, with object name "templ"
  m.write('load {0}, templ\n'.format(templ_pdb))
  if resid is not None:
    m.write('sele templ_resid, templ and poly and {0}\n'.format(resid))
  else:
    m.write('sele templ_resid, templ and poly\n')

  ## Load each Model PDB and superimpose it onto "templ"
  for model in Models:
    mod_name = model.split('/')[-1].split('.')[0]
    m.write('load {0}, {1}\n'.format(model, mod_name))

    m.write('{0} {1}, {2}\n'.format(align_mode, mod_name, "templ_resid"))
    m.write('save {0}.{1}, {2}\n'.format(mod_name, extOt, mod_name))
  
  m.write('hide everything\nshow ribbon\nshow lines, org\n')
  m.write('color white, all\ncolor red, templ\ncenter templ\n')
  m.write('save {0}.pse\n'.format(pymol_pref))
  m.close()

  os.system('{0} -c {1}.pml > {1}.pymol-log; wait'.format(pymol_exec, pymol_pref))
#  os.system('bzip2 -f {0}.pse *.{1}'.format(pymol_pref, extOt))


##########################################################################

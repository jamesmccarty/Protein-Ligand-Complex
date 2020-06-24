#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='combine gromacs drug and protein .gro and .top files')

parser.add_argument('-protein',help='name of protein .gro file',type=str,default='conf.gro')
parser.add_argument('-ligand',help='name of ligand .gro file',type=str,default='ligand.gro')
parser.add_argument('-proteintop',help='name of protein .top file',type=str,default='topol.top')
parser.add_argument('-ligandtop',help='name of ligand .top file',type=str,default='ligand.top')

args=parser.parse_args()

print 'Protein .gro file is ', args.protein
print 'Ligand .gro file is ', args.ligand
print 'Protein .top file is ', args.proteintop
print 'Ligand .top file is ', args.ligandtop

fprotein = open(args.protein,'r')
fligand = open(args.ligand,'r')

fgromacs = open('complex.gro','w')

## First generate combined gromacs structure file

line=fprotein.readline()
Nprotein=fprotein.readline()

n_protein_atoms = int(Nprotein)

print 'Reading ', n_protein_atoms, 'protein atoms'

line=fligand.readline()
Nligand=fligand.readline()

n_ligand_atoms = int(Nligand)

print 'Reading ', n_ligand_atoms, 'ligand atoms'

n_atoms = n_protein_atoms + n_ligand_atoms

fgromacs.write('protein ligand complex \n')
fgromacs.write(str(n_atoms)+'\n')

for i in range(n_protein_atoms):
    line=fprotein.readline()
    fgromacs.write(line)

for i in range(n_ligand_atoms):
    line=fligand.readline()
    fgromacs.write(line)

line=fprotein.readline()
fgromacs.write(line)

print 'generated combined structure file', 'complex.gro'

fprotein.close()
fligand.close()
fgromacs.close()

## Now create topology file

fproteintop = open(args.proteintop,'r')
fligandtop = open(args.ligandtop,'r')
ftopology = open('complex.top','w')
fliganditp = open('ligand.itp','w')
fproteinitp = open('protein.itp','w')

ftopology.write(';\n')
ftopology.write(';   topology file\n')
ftopology.write('\n')

found_force_field = False

for line in fproteintop:
    line=line.strip()
    if (line=='; Include forcefield parameters'):
        ftopology.write(line)
        ftopology.write('\n')
        ffparams=(next(fproteintop,'').strip())
        ftopology.write(ffparams)
        ffparams = ffparams.split()
        ftopology.write('\n\n')
        found_force_field = True
        break

if(found_force_field):
    print 'protein force field parameters found: ', ffparams[1]
else:
    print 'WARNING: protein force field parameters not found'

read_atom_types = False
atomcount = 0

for line in fligandtop:
    line=line.strip()
    if (line=='[ atomtypes ]'):
        read_atom_types = True
    if (read_atom_types):
        atomcount = atomcount + 1
        ftopology.write(line)
        ftopology.write('\n')
        if(len(line)) == 0:
            read_atom_types = False
            ftopology.write('\n')

print 'Writing ',atomcount, 'atomtypes'

fligandtop.seek(0)

ftopology.write('; Include chain topologies\n')
ftopology.write('#include "protein.itp"\n')
ftopology.write('#include "ligand.itp"\n\n')

fproteintop.seek(0)

found_water_ffield = False

for line in fproteintop:
    line=line.strip()
    if (line=='; Include water topology'):
        ftopology.write(line)
        ftopology.write('\n')
        ffparams=(next(fproteintop,'').strip())
        ftopology.write(ffparams)
        ffparams = ffparams.split()
        ftopology.write('\n\n')
        ftopology.write('#ifdef POSRES_WATER\n')
        ftopology.write('; Position restraint for each water oxygen\n')
        ftopology.write('[ position_restraints ]\n')
        ftopology.write(';  i funct       fcx        fcy        fcz\n')
        ftopology.write('   1    1       1000       1000       1000\n')
        ftopology.write('#endif\n\n')
        found_water_ffield = True
        break

if(found_water_ffield):
    print 'water force field parameters found: ', ffparams[1]
else:
    print 'WARNING: water force field parameters not found'

fproteintop.seek(0)

found_ion_ffield = False

for line in fproteintop:
    line=line.strip()
    if (line=='; Include topology for ions'):
        ftopology.write(line)
        ftopology.write('\n')
        ffparams=(next(fproteintop,'').strip())
        ftopology.write(ffparams)
        ffparams = ffparams.split()
        ftopology.write('\n\n')
        found_ion_ffield = True
        break

if(found_ion_ffield):
    print 'ion force field parameters found: ', ffparams[1]
else:
    print 'WARNING: ion force field parameters not found'

ftopology.write('[ system ]\n')
ftopology.write('Protein and ligand\n\n')
ftopology.write('[ molecules ]\n')
ftopology.write('; Compound        #mols\n')

fproteintop.seek(0)

found_protein = False

for line in fproteintop:
    line=line.strip()
    if (line=='; Compound        #mols'):
        ffparams=(next(fproteintop,'').strip())
        ftopology.write(ffparams)
        ffparams = ffparams.split()
        ftopology.write('\n')
        found_protein = True
        break

if(found_protein):
    print 'protein chain is: ', ffparams[0]
else:
    print 'WARNING: protein chain not found in toplogy'

fproteintop.seek(0)

read_molecule_types = False
molecule_count = 0

for line in fproteintop:
    line=line.strip()
    if (line=='[ moleculetype ]'):
        read_molecule_types = True
        molecule_count = molecule_count + 1
    if (read_molecule_types):
        fproteinitp.write(line)
        fproteinitp.write('\n')
        if(line == '#endif'):
            read_molecule_types = False
            fproteinitp.write('\n\n')

print 'Writing ',molecule_count, 'molecules in protein.itp file'

read_molecule_types = False
molecule_count = 0

for line in fligandtop:
    line=line.strip()
    if (line=='[ moleculetype ]'):
        read_molecule_types = True
        molecule_count = molecule_count + 1
    if (read_molecule_types):
        if(line == '[ system ]'):
            read_molecule_types = False
            fliganditp.write('\n\n')
        else:
            fliganditp.write(line)
            fliganditp.write('\n')

print 'Writing ',molecule_count, 'molecules in ligand.itp file'

fligandtop.seek(0)

found_ligand = False

for line in fligandtop:
    line=line.strip()
    if (line=='; Compound        nmols'):
        ffparams=(next(fligandtop,'').strip())
        ftopology.write(ffparams)
        ffparams = ffparams.split()
        ftopology.write('\n')
        found_ligand = True
        break

if(found_ligand):
    print 'ligand molecule is: ', ffparams[0]
else:
    print 'WARNING: ligand molecule not found in toplogy'

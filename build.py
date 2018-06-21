#!/usr/bin/python

import os
import numpy as np
import sys

class MOL:
    def __init__(self,paths, opts,typemol = 'reac'):
        """
        MOL object is reac, prod, or ts.
        """ 
        self.typemol    = typemol
        self.paths      = paths
        self.convert    = paths['x2z'] 
        global ob, pa, io, qc
        sys.path.insert(0, paths['qtc'])
        import obtools as ob
        import iotools as io
        import patools as pa
        import qctools as qc

        #OPTIONS##############################
        self.nsamps     = opts[0]    #number of MonteCarlo sampling points
        self.abcd       = opts[1]
        self.nodes      = opts[2]
        self.interval   = opts[3]    #Interval for tors scan (should be same geometry as 0 degree)
        self.nsteps     = opts[4]    #Number of points to take on PES
        self.XYZ        = opts[5]    #if QTC provides XYZ, 'true' (smiles.xyz), logfile 'anything.log', or just use smiles 'false'
        self.xyzstart   = opts[6]
        self.MDTAU      = opts[7]
        ######################################
        self.ijk        = [0, 0, 0]
        self.sort       = None
        
    
    def build_xyzfile(self, smiles):

        """
        Finds user-specified xyz or geo coordinates or generates openbabel coordinates and places
        them in data/<smiles>.xyz for x2z to use
 
        INPUT:
        smiles         -- SMILES name for molecule
        OUTPUT:
        smilesfilename -- name of cartesian coordinate file
        """
        smilesfilename = ob.get_smiles_filename(smiles) 
        cartlines = ''
        if '_m' in smiles:
            smiles, self.mult = smiles.split('_m')
        else:
            self.mult = ob.get_multiplicity(ob.get_mol(smiles))

        if self.XYZ.lower() == 'false':
            cartlines = build_obcart(smiles,self.mult)

        elif '.log' in self.XYZ.lower():
            cartlines = io.read_file('../' + self.XYZ)
            cartlines = pa.gaussian_xyz_foresk(cartlines)

        elif '.xyz' in self.XYZ.lower():
            if io.check_file(self.XYZ):
                cartlines = io.read_file(self.XYZ)
            else:
                print('ERROR: no .geo or .xyz provided')
                print('...Using openbabel instead')
                cartlines = build_obcart(smiles, self.mult)
        
        elif len(self.XYZ.split('/') ) < 2:
            cartlines, ijk = read_cart(smiles, self.mult)
            if ijk[0] != 0:
                self.ijk = ijk

        elif len(self.XYZ.split('/')) > 2:  
            self.XYZ = self.XYZ.replace('g09','gaussian')

            if io.check_file(io.db_opt_path(self.XYZ.split('/')[0], self.XYZ.split('/')[1], self.XYZ.split('/')[2], None, smiles) + '/' + smilesfilename + '.geo'):
                cartlines = io.db_get_opt_prop(smiles, 'geo', None, self.XYZ.split('/')[0], self.XYZ.split('/')[1], self.XYZ.split('/')[2])
                cartlines = cartlines.replace('\n\n','')
                cartlines = str(len(cartlines.split('\n'))) + ' \n\n' + cartlines
            elif io.check_file(io.db_opt_path(self.XYZ.split('/')[0], self.XYZ.split('/')[1], self.XYZ.split('/')[2], None, smiles) + '/' + smilesfilename + '.xyz'):
                cartlines = io.db_get_opt_prop(smiles, 'xyz', None, self.XYZ.split('/')[0], self.XYZ.split('/')[1], self.XYZ.split('/')[2]).split('\n\n')[1]
                cartlines = cartlines.split('\n')[0] + ' \n\n' + cartlines

            else:
                print ('\nERROR: No geometry found at ' + io.db_opt_path(self.XYZ.split('/')[0], self.XYZ.split('/')[1], self.XYZ.split('/')[2], None, smiles)
                               + smilesfilename + '.xyz' + '\n...building using OpenBabel instead')
                cartlines = build_obcart(smiles,self.mult)
        else:
            print('\nERROR: You have not specified a valid way to get the coordinates.  Use false, true, smiles.log, smiles.geo, smiles.xyz, or prog/method/basis')
        io.write_file(cartlines,smilesfilename + '.xyz')
        return smilesfilename

    def cart2zmat(self,smiles): 
        """
        Runs x2z by Yuri Georgievski on a file of Cartesian coordinates and collects
        the internal coordinate and torsional angle information that it outputs
        """
        import re

        ######   initialize  ###############
        ####################################
        atoms   = []
        measure = []
        angles  = []
        consts  = []

        smilesfilename = self.build_xyzfile(smiles)
        mol         = ob.get_mol(smiles,make3D=True)
        self.charge = ob.get_charge(mol)
        self.stoich = ob.get_formula(mol
)
        #######   RUN X2Z   #################
        #####################################
        tempfile = 'temp'
        paths = self.paths
        gcc     = paths['gcc']
        intel   = paths['intel']

        if io.check_file(smilesfilename + '.xyz'):
            os.system('{0}; {1}; {2} {3}.xyz > {4}'.format(gcc, intel, self.convert, smilesfilename,tempfile))
        else:
            os.system('{0}; {1}; {2} {3}_m{4}.xyz >  {5}'.format(gcc, intel, self.convert, smilesfilename, str(self.mult),tempfile))

        if os.stat(tempfile).st_size < 80:
            print('Failed')
            print('Please check that directory name and cartesian coordinate file name are equivalent')
            print('Please check that test_chem is in location: ' +  self.convert)
            print('Using OpenBabel zmat: no hindered rotors will be specified')
            atoms, measure = build_obzmat(smiles)                     #If test_chem fails use openbabel to get zmat
            self.symnum = ' 1'
            self.ilin   = ' 0'
            if len(atoms) < 3:                                        #Linear if dihedral (NEEDS TO ACTUALLY BE COMPUTER)
                self.ilin =' 1' 
            return atoms, measure, angles

        ####Get relevant data from Test_Chem output file########
        ########################################################
        props, lines = io.read_file(tempfile).split('Z-Matrix:\n')
        io.rm(tempfile)
        
        #linearity
        self.ilin =' 1'
        result = re.search("molecule is (\w+)", props).groups()[0]
        if ('non' in result or 'plan' in result):
            self.ilin = ' 0'

        #enantiomers
        ent = False
        result = re.search("has enantiomer?(\w+)", props)
        if result: 
            result = result.groups()[0]
            if ('yes' in result):
                ent = True

        #Symmetry number
        self.symnum = re.search("symmetry number = (\w+)", props).groups()[0]#Symmetry factor

        #Rearrangement
        props, order = props.split('Z-Matrix atom order:')
        self.sort = {}
        for index in order.split('\n')[1:-2]:
            self.sort[int(index.split('>')[1])] = str(index.split('-->')[0])
        #rotational groups
        groups = ''
        if len(lines.split("Rotational groups:")) > 1:
            lines, groups = lines.split("Rotational groups:") #if x2z prints rot groups

        #zmatrix connectivity
        lines = lines.split('\n')
        for i,line in enumerate(lines):
            if line == '':
                break
            atoms.extend([line.rstrip('\n').replace(' ','').split(',')])     #Main part of ZMAT

        #zmatrix parameters
        for j in range(i+1,len(lines)):
            if lines[j]:
                if 'const' in lines[j].lower() or 'Rot' in lines[j] or 'Beta' in lines[j]:
                    break
            measure.extend(lines[j].replace('=',' ').rstrip('\n').split())   #Gets parameters

        #Hindered rotor dihedral angles and constant angles
        for n in range(len(lines[j:])):
            if "Const" in lines[j+n] and lines[j].split(':')[1].rstrip('\n').strip() != '':
                consts  = lines[j+n].split(':')[1].strip().split()   
                if "Rot" in lines[j+n] and lines[j+n].split(':')[1].rstrip('\n').strip() != '':
                    angles  = lines[j+n].replace(" ","").upper().split(':')[1].strip().split(',')   
            elif "Rot" in lines[j+n] and lines[j+n].split(':')[1].rstrip('\n').strip() != '':
                    angles  = lines[j+n].replace(" ","").upper().split(':')[1].strip().split(',')  
        
        #Reformat zmatrix parameters
        measure = np.array(measure)                     
        if  (len(measure)%2 != 0):
            measure = measure[:-1]
        measure = measure.reshape( len(measure)/2, 2)      #Puts measurements into two columns
        for angle in measure:
            if 'R' in angle[0]:
                angle[1] = str(float(angle[1]) * 0.529177) #bohr to angstrom

        #Put dummy atoms paramters inside zmat
        for i, row in enumerate(atoms):
            if 'X' in row[0]:
                j=0
                for meas in measure:
                    if   meas[0] == row[2]:
                        atoms[i][2] = '{:.6f}'.format(float(meas[1]))
                        measure = np.delete(measure, j, axis=0)
                        j-=1
                    elif meas[0] == row[4]:
                        atoms[i][4] = '{:.2f}'.format(float(meas[1]))
                        measure = np.delete(measure, j, axis=0)
                        j-=1
                    elif len(row) > 5:
                         if meas[0] == row[6]:
                            atoms[i][6] = '{:.2f}'.format(float(meas[1]))
                            measure = np.delete(measure, j, axis=0)
                            j-=1
                    j+=1

        #Put constant angles inside zmat 
        for angle in consts:
            for i, row in enumerate(atoms):
                if len(row) > 5:
                    if angle in row[6]:
                        j=0
                        for meas in measure:
                            if  meas[0] == angle:
                                atoms[i][6] = '{:.6f}'.format(float(meas[1]))
                                measure = np.delete(measure, j, axis=0)
                                j-=1
                            j+=1
                elif len(row) > 3:
                    if angle in row[4]:
                        j=0
                        for meas in measure:
                            if  meas[0] == angle:
                                atoms[i][4] = '{:.6f}'.format(float(meas[1]))
                                measure = np.delete(measure, j, axis=0)
                                j-=1
                            j+=1
        #Count rotors
        nmethylgroups = 0 
        for rotor in groups.split('\n'):
            grouplist = rotor[3:].lower().split()
            if 'c1h3' in grouplist:
               nmethylgroups += 1
        self.nrotors = len(angles)
        self.nrotors = self.nrotors - nmethylgroups

        return atoms, measure, angles
         

    def build(self, n, smiles, angles, atoms = [], measure = []):
        """ 
        Builds reacn.dat or prodn.dat for EStokTP withh user defined nosmps (Monte Carlo sampling points
        for geometry search) and nhindsteps (number of points on the PES) 
        """
        smilesfilename = ob.get_smiles_filename(smiles)
        
        if self.typemol == 'reac' or self.typemol == 'prod':
            atoms, measure, angles  = update_interns(n,atoms,measure,angles)
            if not self.nsamps:
                if len(self.abcd.split(',')) >3:
                    a, b, c, d = self.abcd.split(',')
                    a = int(a)
                    b = int(b)
                    c = int(c)
                    d = int(d)
                    self.nsamps = min(a + b * c**self.nrotors, d)
            self.nsamps = str(int(np.ceil(float(self.nsamps) / self.nodes)))
            #MC Parameters#############
            zmatstring  = 'nosmp dthresh ethresh\n'     
            zmatstring += self.nsamps + '  1.0  0.00001\n'
            #Torsional Scan Parameters#############
            zmatstring += tau_hind_str(atoms, angles, self.interval, self.nsteps, self.MDTAU)
            #Size and linearity of molecule###########
            zmatstring += '\nnatom natomt ilin\n'
            ndummy = count_dummy(atoms)
            zmatstring += str(len(atoms)-ndummy) + ' ' + str(len(atoms)) + self.ilin + '\n'
        else:
            if not self.nsamps:
                if len(self.abcd.split(',')) >3:
                    a, b, c, d = self.abcd.split(',')
                    a = int(a)
                    b = int(b)
                    c = int(c)
                    d = int(d)
                    self.nsamps = min(a + b * c**1, d)
            self.nsamps = str(np.ceil(float(self.nsamps) / self.nodes))
            #MC Parameters#############
            zmatstring  = 'nosmp dthresh ethresh\n'     
            zmatstring += self.nsamps + '  1.0  0.00001\n'
            #Torsional Scan Parameters#############
            zmatstring += tau_hind_str(atoms, angles, self.interval, self.nsteps, self.MDTAU)
            if n == 'ts':
                import get_sites
                #i,j,k sites###########################
                zmatstring += '\nisite jsite ksite\n'
                if self.sort:
                    for i in range(3):
                        self.ijk[i] = self.sort[int(self.ijk[i])-1]
                if self.ijk[0] != 0:
                    zmatstring += ' '.join(self.ijk)
                else:
                    lines = io.read_file('../rmg.dat')
                    zmatstring += ' '.join(get_sites.sites(lines))
                zmatstring += '\n\nrmin rmax nr\n 1.0 2.5 8\n  -->aabs1,babs1,aabs2,babs2,babs3\n 90., 180., 90., 175., 90.\n'

        #Typical Z-Matrix########################
        zmatstring += '\ncharge  spin  atomlabel\n'
        zmatstring += str(self.charge) + ' ' + str(self.mult) + '\n'

        if self.typemol == 'reac' or self.typemol == 'prod':
            for row in atoms:
                for j in range(len(row)):
                    zmatstring += row[j] + ' '
                zmatstring +='\n'
            zmatstring += '\nintcoor'

            deletedangles=[]
            for hin in angles:                                       #Deletes hindered angles from
                for i,meas in enumerate(measure):                    #The internal coordinate list
                    if hin.lower().strip() == meas[0].lower().strip():
                        deletedangles.append(measure[i][1]) 
                        measure = np.array([np.delete(measure.T[0],i),np.delete(measure.T[1],i)]).T
            for meas in measure:
                zmatstring += '\n' + meas[0] + ' ' + meas[1]
            
            zmatstring += '\n'
        #Sym factor and no. of electronic states#
        zmatstring += '\nSymmetryFactor\n' + self.symnum + '\n'
        zmatstring += '\nnelec\n1\n 0.  ' + str(self.mult) + '\n\nend\n'

        #Build Reac/Prodnum_opt.out for starting after level0 or level1
        if self.xyzstart == '0' or self.xyzstart == '1':
            optim = build_optout(self.xyzstart, self.XYZ, measure, angles, deletedangles, smiles, smilesfilename)
            io.write_file(optim, '../output/' + self.typemol + str(n) + '_opt.out')


        return zmatstring

def build_molpro(meth,freqcalc,opt):
    """
    Builds molpro theory file 
    """
    method, basis = meth[2].split('/')
    nondft = ['hf', 'ccsd', 'cisd','mp']
    dft = True
    for key in nondft:
        if key in method.lower():
            dft = False
    molstr  = 'nosym\nEnd1\n\n'
    molstr += '!closed shell input\n\nbasis=' + basis +'\n'
    if dft:
        molstr += 'dft=['+method+']\n'
        molstr += 'hf\ndft'
    else:
        molstr += 'hf\n' + method.lower()
    if opt:
        molstr += '\noptg'
    if freqcalc:
        molstr += '\nfrequencies,symm=auto,numerical'
    molstr += '\nENERGY=energy\n\n'
    molstr += 'CBSen=energy\n\n'
    molstr += '! these lines must be always included in molpro input\n'
    molstr += '! CBSen should be defined as desired\n'
    molstr += '! the molden line should be left as it is\n\nput,molden,molpro.molden\n\n---\n\n'
    molstr += 'End2\n\n!open shell input\n\n'
    molstr += 'basis=' + basis +'\n'
    if dft:
        molstr += 'dft=['+method+']\n'
        molstr += 'uhf\ndft'
    elif 'ccsd' in method:
        molstr += 'rhf\n' + 'u' + method.lower()
    else:
        molstr += 'uhf\n' + method.lower()
    if opt:
        molstr += '\noptg'
    if freqcalc:
        molstr += '\nfrequencies,symm=auto,numerical'
    molstr += '\nENERGY=energy\n\n'
    molstr += 'CBSen=energy\n\n'
    molstr += '! these lines must be always included in molpro input\n'
    molstr += '! CBSen should be defined as desired\n'
    molstr += '! the molden line should be left as it is\n\nput,molden,molpro.molden\n\n---\n\n'
    molstr += 'End3\n\n\n\n'

    molfile = meth[0]
    if 'hind' in molfile:
        molfile = 'onedtau'
    if 'hlevel' in molfile:
        molfile = 'hl'
    if 'symm' in molfile:
        molfile = 'symm'
    molfile += '_molpro.dat'
    return molstr, molfile


def build_theory(meths,nTS, optoptions):
    """
    Builds theory.dat 
    meth[0] is module, meth[1] is program, and meth[2] is theory/basis
    """
    theory    = ''
    tsopt  = ' opt=(ts,calcfc,noeig,intern,maxcyc=50)\n '
    rpopt  = ' opt=(' + optoptions +  ')\n '
    vwopt  = ' opt=(internal,calcall) scf=qc\n '
    allint = 'int=ultrafine nosym '

    for meth in meths:
        if 'molpro' in meth[1]:

            theory += meth[0] + ' ' + meth[1] + '\n\n'
            if meth[0] == 'level1':
                freqcalc = True
            else:
                freqcalc = False
            if meth[0] == 'hlevel':
                opt = False
            else:
                opt = True
            molpro  = build_molpro(meth,freqcalc, opt)
            io.write_file(molpro[0], molpro[1])

        elif  'g09' in meth[1]:   
            theory += meth[0] + ' ' + meth[1] + '\n '
            theory += meth[2] + rpopt + allint
            if meth[0] == 'level1':
                theory += ' freq'
            theory += '\n\n'
            if nTS > 2 and meth[0] == 'level1':
                theory += meth[0] + '_61 ' + meth[1] + '\n '
                theory += meth[2] + vwopt + allint
                theory += ' freq\n\n'
            if nTS > 1 and meth[0] == 'level1':
                theory += meth[0] + '_51 ' + meth[1] + '\n '
                theory += meth[2] + vwopt + allint
                theory += ' freq\n\n'
            if nTS > 0 and meth[0] != 'hlevel' and meth[0] != 'irc':
                theory += meth[0] + '_ts ' + meth[1] + '\n '
                theory += meth[2] + tsopt + allint
                if meth[0] == 'level1':
                    theory += ' freq'
                if meth[0] == 'hind_rotor':
                    theory += '\n ' + meth[2] + rpopt + allint
                theory += '\n\n'
        elif meth[0].lower() != 'ktp':
            print meth[0] + ' is not a recognized program.\n'
            
    theory += 'End'

    return theory 

def build_estoktp(params, jobs, nreacs, nprods, nTS):
    """
    Builds estoktp.dat
    """
    stoichs    = params[0]
    reactype  = params[1]
    coresh    = params[2]
    coresl    = params[3]
    mem       = params[4]

    fullstoich = {}
    if len(stoichs) > 2:
        stoichs = stoichs[:2]
    for stoich in stoichs:
        stoich = qc.get_atom_stoich(stoich)
        if 'C' in stoich:
           if 'C' in fullstoich:
               fullstoich['C'] += stoich['C']
           else:
               fullstoich['C']  = stoich['C']
        if 'H' in stoich:
           if 'H' in fullstoich:
               fullstoich['H'] += stoich['H']
           else:
               fullstoich['H']  = stoich['H']
        if 'O' in stoich:
           if 'O' in fullstoich:
               fullstoich['O'] += stoich['O']
           else:
               fullstoich['O']  = stoich['O']
        if 'N' in stoich:
           if 'N' in fullstoich:
               fullstoich['N'] += stoich['N']
           else:
               fullstoich['N']  = stoich['N']
        if 'S' in stoich:
           if 'S' in fullstoich:
               fullstoich['S'] += stoich['S']
           else:
               fullstoich['S']  = stoich['S']
    stoich = ''
    if 'C' in fullstoich:
       if fullstoich['C'] > 1:
           stoich += 'C{:g}'.format(fullstoich['C'])
       else:
           stoich += 'C'
    if 'H' in fullstoich:
       if fullstoich['H'] > 1:
           stoich += 'H{:g}'.format(fullstoich['H'])
       else:
           stoich += 'H'
    if 'O' in fullstoich:
       if fullstoich['O'] > 1:
           stoich += 'O{:g}'.format(fullstoich['O'])
       else:
           stoich += 'O'
    if 'N' in fullstoich:
       if fullstoich['N'] > 1:
           stoich += 'N{:g}'.format(fullstoich['N'])
       else:
           stoich += 'N'
    if 'S' in fullstoich:
       if fullstoich['S'] > 1:
           stoich += 'S{:g}'.format(fullstoich['S'])
       else:
           stoich += 'S'
    eststring = ' Stoichiometry\t' + stoich.upper()
    
    PossibleRxns = ['addition','abstraction','isomerization','betascission','well','']
    Tstype       = ['TS','wellr','wellp']
 
    if reactype.lower() not in PossibleRxns:
        print('ReactionType ' + reactype +' is unrecognized, please use: Addition, Abstraction, Isomerization, or Betascission')
    elif reactype != '' and reactype.lower() != 'well':
        eststring  += '\n ReactionType   ' + reactype
        if nTS != 0:
            eststring += '  ' + str(nTS) + 'TS'
            if nTS > 1:
                eststring += '\n WellR findgeom level1'
                if nTS > 2:
                    eststring += '\n WellP findgeom level1'
    if 'Irc' in jobs:
        eststring += '\n Variational'
    if nprods > 0:
        eststring += '\n Prods'
    eststring +='\n Debug  2\n'
 
    for job in jobs:
        if 'kTP' in job or 'irc' in job.lower():
            eststring += job
        else:
            for n in range(nreacs):
                if 'Opt_1' in job:
                    eststring += '\n ' + job.rstrip('_1') + '_Reac' + str(n+1) + '_1'
                else:
                    eststring += '\n ' + job + '_Reac' + str(n+1)
 
            for n in range(nprods):
                if 'Opt_1' in job:
                    eststring += '\n ' + job.rstrip('_1') + '_Prod' + str(n+1) + '_1'
                else:
                    eststring += '\n ' +job + '_Prod' + str(n+1)
 
            for n in range(nTS):
                if 'Opt_1' in job:
                    eststring += '\n ' + job.rstrip('_1') + '_' + Tstype[n]  + '_1'
                elif 'Tau' in job:
                    if n < 1:
                        eststring += '\n ' +job + '_' + Tstype[n] 
                elif 'Opt' in job and n <1:
                    eststring += '\n Grid_' + job + '_' + Tstype[n]
                    eststring += '\n ' + job + '_' + Tstype[n] + '_0'
                    if 'nOpt' in job:
                        eststring += '\n nTauo_' + Tstype[n]
                    else:
                        eststring += '\n Tauo_' + Tstype[n]
                elif 'Opt' not in job:
                    eststring += '\n ' + job + '_' + Tstype[n]
                    
        eststring += '\n'
 
    eststring += '\nEnd'
    eststring += '\n ' + coresh + ',' + coresl + '\n numprocll,numprochl\n'
    eststring += ' ' + mem + 'MW  ' + mem + 'MW\n gmemll gmemhl\n'
 
    return eststring

def is_dummy(atom):
   if 'xe' in atom.lower():
      return False
   elif 'x' in atom.lower():
      return True
   return False

def count_dummy(atoms):
  dummy = 0
  for atom in atoms:
      if is_dummy(atom[0]):
          dummy += 1
  return dummy

def update_interns(n, atoms, measure ,angles):
    """
    Converts internal coordinate information from x2z to the form for
    a Z-Matrix required to run EStokTP
    """
    if len(atoms) == 0:
        return atoms, measure, angles
    for index,atom in enumerate(atoms):
        if n > 1:
            atoms[index][0] = atom[0].lower() + str(index+99)
        else:
            atoms[index][0] = atom[0].lower() + str( index+1)
        if len(atom) > 1:
            atoms[index][1] = atoms[int(atom[1])-1][0]
            if n > 1:
                atoms[index][2] = atoms[index][2].upper().replace('R','R10')
            if len(atom) > 3:
                atoms[index][3] = atoms[int(atom[3])-1][0]
                if n > 1:
                    atoms[index][4] = atoms[index][4].upper().replace('A','A10')
                if len(atom) > 5:
                    atoms[index][5] = atoms[int(atom[5])-1][0]
                    if n > 1:
                        atoms[index][6] = atoms[index][6].upper().replace('D','D10')
    if n > 1:
        for i in range(len(measure)):
            measure[i][0] = measure[i][0].upper().replace('R','R10')
            measure[i][0] = measure[i][0].upper().replace('A','A10')
            measure[i][0] = measure[i][0].upper().replace('D','D10')
        for i in range(len(angles)):
            angles[i] = angles[i].upper().replace('D','D10')
    return atoms, measure, angles

def find_period(zmat,hin):
    """
    Rough way of determining internal symmetry number (Hydrogen counting)
    """
    sym1 = sym2 = k = 0
    for i, row in enumerate(zmat[3:]):
        if hin.upper() in row[6]:
            k = i+3
    atom1 = zmat[k][1]
    atom2 = zmat[k][3]
    for row in zmat[1:]:
        if atom1 == row[1]:
            if 'h' in row[0]:
                sym1 += 1
            elif atom2 != row[0]:
                sym1 -= 100
        elif atom1 == row[0]:
            if 'h' in row[1]:
                sym1 += 1
            elif atom2 != row[1]:
                sym1 -= 100
        elif atom2 == row[1]:
            if 'h' in row[0]:
                sym2 += 1
            elif atom1 != row[0]:
                sym2 -= 100
        elif atom2 == row[0]:
            if 'h' in row[1]:
                sym2 += 1
            elif atom1 != row[1]:
                sym2 -= 100
    if sym1 <= 1 or sym2 <=1:
        period =  max(sym1, sym2)
    elif sym1 == sym2:
        period =  sym1 #+ sym2
    else:
        period = sym1 * sym2
    if period < 1:
        period = 1
    return period

def tau_hind_str(atoms, angles, interval, nsteps, mdtau):
    #TAU
    string  = '\nntau number of sampled coordinates\n'
    string += str(len(angles)) + '\n'
    string += ' -->nametau, taumin, taumax\n'
    for angle in angles:
        periodicity = find_period(atoms, angle)
        string += angle + ' 0 ' + str(interval) + '\n'

    #1 and 2D HIND
    string += '\nnhind\n'
    string += str(len(angles)) + '\n'
    string += ' -->namehind,hindmin,hindmax,nhindsteps,period\n'
    for hin in angles:
        periodicity = find_period(atoms, hin)
        string += hin + ' 0 ' + str(float(interval)/periodicity)  + ' ' + str(int(round(float(nsteps)/periodicity))) + ' ' + str(periodicity)  + '\n'   
    if mdtau and len(angles) > 1:
        mdtau   = mdtau.strip('D').strip('d')
        string += '\nnhind' + mdtau + 'D\n'
        string += '1\n'
        string += ' -->namehind,hindmin,hindmax,nhindsteps,period\n'
        for i in range(int(mdtau)):
            periodicity = find_period(atoms, angles[i])
            string += angles[i] + ' 0 ' + str(float(interval)/periodicity) + ' ' + str(int(round(float(nsteps)/periodicity))) + ' ' + str(periodicity)  + '\n'   
    return string

def build_obzmat(smiles):
    """
    Uses QTC interface by Murat Keceli to Openbabel to generate zmat file based on smile string
    
    INPUT:
    smiles -- SMILES string for molecule
    OUTPUT:
    atoms   -- The atoms and their connection to previous atoms (in EStokTP format)
    measure -- The distance and angle parameters for the zmat   (in EStokTP format)
    """
    
    mol   = ob.get_mol(smiles,make3D=True)
    zmat  = ob.get_zmat(  mol)
    atoms, measure = zmat.split('\nVariables:\n')
    atoms = atoms.split('\n')
    for i in range(len(atoms)):
        atoms[i] = atoms[i].split()
    measure = measure.split('\n')
    del measure[-1]
    for i in range(len(measure)):
        measure[i] = measure[i].upper().split('= ')
    return atoms, measure

def build_obcart(smiles, mult):
    """
    Uses QTC interface by Murat Keceli to Openbabel to generate cartesian coordinates based 
    on SMILES string

    INPUT:
    smiles -- SMILES string for molecule
    mult   -- mult for molecule
    OUTPUT:
    lines  -- cartesian coordinates, formatted for a xyz file
    """
    filename    =  ob.get_smiles_filename(smiles) + '_m' + str(mult) + '.xyz'
    mol         =  ob.get_mol(smiles,make3D=True)
    lines       =  ob.get_xyz(mol)#.split('\n')
    return lines
 
def build_optout(xyzstart, XYZ, measure, angles, deletedangles, smiles, smilesfilename):
    """
    Creates the output/reac1_opt.out input for when an initial xyz is specified to start 
    EStokTP with at either level0 or level1
    """
    optim = ''
    if '0' in xyzstart:
        XYZ = XYZ.replace('g09','gaussian')
        optim = 'opt geom           1'
        for meas in measure:
            optim += '\n\t' + meas[1]              
        for i in range(len(angles)):
            optim += '\n\t{}'.format(deletedangles[i])             

        if '.log' in XYZ:
            E= str(pa.gaussian_energy(io.read_file('../' + XYZ)))
        elif '.xyz' in XYZ:
            E = io.read_file(XYZ).split('\n')[1]
        elif io.check_file('../' + smilesfilename + '.ene'):
            E = io.read_file('../' + smilesfilename + '.ene')
        elif len(XYZ.split('/')) > 2 :
            split = XYZ>split('/')
            E = io.db_get_sp_prop(smilesfilename, 'ene', None, split[0], split[1], split[2],
                                                         split[0], split[1], split[2])
            if E == None:
                print('\nNo energy found at ' + io.db_sp_path(split[0], split[1], split[2], 
                        None, smiles, split[0], split[1], split[2]) + smilesfilename + '.ene')
        else:
            E = '0'

        if E == None:
           if io.check_file(XYZ.replace('xyz','ene')):
               E = io.read_file(XYZ.replace('xyz','ene'))
        optim += '\n\t' + E 

    elif '1' in xyzstart:
        optim = 'opt level1 0'
        for meas in measure:
            optim += '\n\t' + meas[1]                
        for i in range(len(angles)):
            optim += '\n\t{}'.format(deletedangles[i])             
 
    return optim

def read_cart(smiles, mult):
    """
    Finds and Reads an initial XYZ file (and set i,j,k sites for a transition state)
    """

    smilesfilename = ob.get_smiles_filename(smiles)
    cartlines = ''
    ijk = [0, 0, 0]
    import os
    if io.check_file('../' + smilesfilename + '_m' + str(mult) + '.xyz'):               
      #  cartlines = io.read_file('../' + smilesfilename + '_m' + str(mult) + '.xyz').split('\n\n')[1]
      #  cartlines =  str(len(cartlines.split('\n'))-1) + ' \n\n' + cartlines
        cartlines = io.read_file('../' + smilesfilename + '_m' + str(mult) + '.xyz')
    elif io.check_file('../' + smilesfilename + '.xyz'):                
        #cartlines = io.read_file('../' + smilesfilename + '.xyz').split('\n\n')[1]
        #cartlines =  str(len(cartlines.split('\n'))-1) + ' \n\n' + cartlines
        cartlines = io.read_file('../' + smilesfilename +  '.xyz')
    elif io.check_file('../' + smilesfilename + '_m' + str(mult) + '.geo'):
        cartlines = io.read_file('../' + smilesfilename + '_m' + str(mult) + '.geo')
        cartlines = str(len(cartlines.split('\n'))-1) + ' \n\n' + cartlines
    elif io.check_file('../' + smilesfilename + '.geo'):
        cartlines = io.read_file('../' + smilesfilename + '.geo')
        cartlines = str(len(cartlines.split('\n'))-1) + ' \n\n' + cartlines
    else:
        print('ERROR: no .geo or .xyz provided')
        print('...Using openbabel instead')
        cartlines = build_obcart(smiles, mult)
        return cartlines, ijk

    #Find if i,j,k site is specified:
    cartlines = cartlines.splitlines()
    for i,line in enumerate(cartlines[1:], start=1):
        if len(line.split()) > 4:
            cartlines[i] = '   '.join(line.split()[1:])
            if line.split()[0] == '1':
                ijk[1] = str(i)
            elif line.split()[0] == '2':
                ijk[0] = str(i)
            elif line.split()[0] == '3':
                ijk[2] = str(i)
            elif line.split()[0] == '4':
                temp = cartlines[i]
                del cartlines[i]
                cartlines.insert(2,temp)
    cartlines = '\n'.join(cartlines)
    return cartlines, ijk 


#!/usr/bin/python

import os
import numpy as np
import sys
sys.path.insert(0, '/home/elliott/Packages/QTC/')
import patools as pa
import iotools as io
import obtools as ob

class MOL:
    def __init__(self,opts,typemol = 'reac'):
         
        self.convert    = io.get_path('/home/elliott/Packages/TorsScan/test_chem') 
        self.typemol    = typemol
        #OPTIONS##############################
        self.nsamps     = opts[0]    #number of MonteCarlo sampling points
        self.interval   = opts[1]    #Interval for tors scan (should be same geometry as 0 degree)
        self.nsteps     = opts[2]    #Number of points to take on PES
        self.XYZ        = opts[3]    #if QTC provides XYZ, 'true' (smiles.xyz), logfile 'anything.log', or just use smiles 'false'
        self.xyzstart   = opts[4]
        self.MDTAU      = opts[5]
        ######################################
        self.ijk        = [0, 0, 0]
        self.sort       = None

    def build_cart(self,smiles):
        """
        Uses QTC interface by Murat Keceli to Openbabel to generate cartesian coorinate file based 
        on SMILES string
        """
        mol         = ob.get_mol(smiles)
        lines       =  ob.get_xyz(mol).split('\n')
        lines[0] = 'Geometry ' + lines[0] + ' Angstrom'
        del lines[1]
        io.write_file('\n'.join(lines), smiles + '.xyz')
        return 
    
    def ob_zmat(self,smiles):
        """
        Uses QTC interface by Murat Keceli to Openbabel to generate zmat file based on smile string
        """
        import obtools as ob
        
        mol   = ob.get_mol(smiles)
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

    def read_cart(self, smiles):
        if io.check_file('../' + smiles + '.xyz'):                
            cartlines = io.read_file('../' + smiles + '.xyz').split('\n\n')[1]
            cartlines = 'Geometry ' + str(len(cartlines.split('\n'))-1) + ' Angstrom\n' + cartlines
            io.write_file(cartlines,smiles + '.xyz')
        elif io.check_file('../' + smiles + '.geo'):
            cartlines = io.read_file('../' + smiles + '.geo')
            cartlines = 'Geometry ' + str(len(cartlines.split('\n'))-1) + ' Angstrom\n' + cartlines
            io.write_file(cartlines,smiles + '.xyz')
        else:
            print('ERROR: no .geo or .xyz provided')
            print('...Using openbabel instead')
            self.build_cart(smiles)
            return 
        #Find if i,j,k site is specified:
        cartlines = cartlines.split('\n')
        for i,line in enumerate(cartlines[1:], start=1):
            if len(line.split()) > 4:
                cartlines[i] = '   '.join(line.split()[1:])
                if line.split()[0] == '1':
                    self.ijk[1] = str(i)
                elif line.split()[0] == '2':
                    self.ijk[0] = str(i)
                elif line.split()[0] == '3':
                    self.ijk[2] = str(i)
                elif line.split()[0] == '4':
                    temp = cartlines[0]
                    cartlines[0] = cartlines[i]
                    del cartlines[i]
                    cartlines.insert(0,temp)
        cartlines = '\n'.join(cartlines)
        io.write_file(cartlines,smiles + '.xyz')
        return 

    def cart2zmat(self,smiles): 
        """
        Runs test_chem by Yuri Georgievski on a file of Cartesian coordinates and collects
        the internal coordinate and torsional angle information that it outputs
        """
        import re

        #initialize########################33
        atoms   = []
        measure = []
        angles  = []

        if self.XYZ.lower() == 'false':
            self.build_cart(smiles)
        elif '.log' in self.XYZ.lower():
            cartlines = io.read_file('../' + self.XYZ)
            io.write_file(pa.gaussian_xyz_foresk(cartlines),smiles + '.xyz')
        elif len(self.XYZ.split('/') ) < 2:
            cartlines = self.read_cart(smiles)
        elif len(self.XYZ.split('/')) > 2:  
            if io.check_file(io.db_opt_path(self.XYZ.split('/')[0], self.XYZ.split('/')[1], self.XYZ.split('/')[2], None, smiles) + '/' + smiles + '.geo'):
                cartlines = io.db_get_opt_prop(smiles, 'geo', None, self.XYZ.split('/')[0], self.XYZ.split('/')[1], self.XYZ.split('/')[2])
                cartlines = 'Geometry ' + str(len(cartlines.split('\n'))-1) + ' Angstrom\n' + cartlines
                io.write_file(cartlines,smiles + '.xyz')
            elif io.check_file(io.db_opt_path(self.XYZ.split('/')[0], self.XYZ.split('/')[1], self.XYZ.split('/')[2], None, smiles) + '/' + smiles + '.xyz'):
                cartlines = io.db_get_opt_prop(smiles, 'xyz', None, self.XYZ.split('/')[0], self.XYZ.split('/')[1], self.XYZ.split('/')[2]).split('\n\n')[1]
                cartlines = 'Geometry ' + str(len(cartlines.split('\n'))-1) + ' Angstrom\n' + cartlines
                io.write_file(cartlines,smiles + '.xyz')
            else:
                print ('\nERROR: No geometry found at ' + io.db_opt_path(self.XYZ.split('/')[0], self.XYZ.split('/')[1], self.XYZ.split('/')[2], None, smiles)
                               + smiles + '.xyz' + '\n...building using OpenBabel instead')
                self.build_cart(smiles)
        else:
            print('\nERROR: You have not specified a valid way to get the coordinates.  Use false, true, smiles.log, smiles.geo, smiles.xyz, or prog/method/basis')
        mol         = ob.get_mol(smiles)
        self.charge = ob.get_charge(mol)
        self.mult   = ob.get_multiplicity(mol)
        self.stoich = ob.get_formula(mol)
        #Peform Test_Chem#########################
        tempfile = 'temp'
        cart     = os.getcwd() + '/' +  smiles.replace('[','\[').replace(']','\]') + '.xyz'
        cart     = cart.replace('(','\(').replace(')','\)')  
        os.system(self.convert + ' ' + cart + ' > ' + tempfile)
        if os.stat(tempfile).st_size < 80:
            print('Failed')
            print('Please check that directory name and cartesian coordinate file name are equivalent')
            print('Please check that test_chem is in location: ' +  self.convert)
            print('Using OpenBabel zmat: no rotation dihedrals will be specified')
            atoms, measure = self.ob_zmat(smiles)                     #If test_chem fails use openbabel to get zmat 
            self.symnum = ' 1'
            self.ilin   = ' 0'
            if len(atoms) < 3:                                        #Linear if dihedral (NEEDS TO ACTUALLY BE COMPUTER)
                self.ilin =' 1' 
            return atoms, measure, angles

        #Get relevant data from Test_Chem output file########
        props,lines = io.read_file(tempfile).split('Z-Matrix:\n')
        props,order = props.split('Z-matrix atom order:')
        self.sort = []
        for index in order.split('\n')[1:-2]:
            self.sort.append(index.split('>')[1])
        io.rm(tempfile)
        lines = lines.split('\n')
        for i,line in enumerate(lines):
            if line == '':
                break
            atoms.extend([line.rstrip('\n').split(',')])                     #Gets connectivity information

        for j in range(i+1,len(lines)):
            if 'Rot' in lines[j]:
                break
            measure.extend(lines[j].replace('=',' ').rstrip('\n').split())   #Gets bond lengths, bond angles, and dihedral angles

        if lines[j].split(':')[1].rstrip('\n').strip() != '':
            angles  = lines[j].split(':')[1].strip().split(',')              #Gets rotational angles to scan

        self.ilin =' 0'
        self.symnum = re.search("symmetry number = (\w+)", props).groups()[0]#Symmetry factor
    
        measure = np.array(measure)                     
        measure = measure.reshape( len(measure)/2, 2)                        #Puts measurements into two columns
        for angle in measure:
            if 'R' in angle[0]:
                angle[1] = str(float(angle[1]) * 0.529177)                   #bohr to angstrom
        
        return atoms, measure, angles 
         

    def build(self, n, smiles, angles, atoms = [], measure = []):
        """ 
        Builds reacn.dat or prodn.dat for EStokTP withh user defined nosmps (Monte Carlo sampling points
        for geometry search) and nhindsteps (number of points on the PES) 
        """
        #Stochastic Geometry Search############
        zmatstring  = 'nosmp dthresh ethresh\n'     
        zmatstring += self.nsamps + '  1.0  0.00001\n'

        if self.typemol == 'reac' or self.typemol == 'prod':

            atoms, measure, angles  = update_interns(atoms,measure,angles)

            #Torsional Scan Parameters#############
            zmatstring += tau_hind_str(atoms, angles, self.interval, self.nsteps, self.MDTAU)
        
            #Size and linearity of molecule###########
            zmatstring += '\nnatom natomt ilin\n'
            ndummy = count_dummy(atoms)
            zmatstring += str(len(atoms)-ndummy) + ' ' + str(len(atoms)) + self.ilin + '\n'


        else:

            sys.path.insert(0, '/home/elliott/scripts')
            import get_sites
            #Torsional Scan Parameters#############
            zmatstring += tau_hind_str(atoms, angles, self.interval, self.nsteps, self.MDTAU)

            if n == 'ts':
                #i,j,k sites###########################
                zmatstring += '\nisite jsite ksite\n'
                if self.ijk[0] != 0:
                    zmatstring += ' '.join(self.ijk)
                else:
                    lines = io.read_file('../rmg.dat')
                    zmatstring += ' '.join(get_sites.sites(lines))
                if self.sort:
                    for i in range(3):
                        self.ijk[i] = str(int(self.sort[int(self.ijk[i])-1])+1)
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
            for hin in angles:                                       #Deletes rotational angles from
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
        if '0' in self.xyzstart:
            optim = 'opt geom           1'
            for meas in measure:
                optim += '\n\t' + meas[1]              
            for i in range(len(angles)):
                optim += '\n\t{}'.format(deletedangles[i])             
            if '.log' in self.XYZ:
                optim += '\n\t' + str(pa.gaussian_energy(io.read_file('../' + self.XYZ)))
            elif io.check_file('../' + smiles + '.ene'):
                E = io.read_file('../' + smiles + '.ene')
            elif len(self.XYZ.split('/')) > 2 :
                E = io.db_get_sp_prop(smiles, 'ene', None, self.XYZ.split('/')[0], self.XYZ.split('/')[1], self.XYZ.split('/')[2]
                                                          ,self.XYZ.split('/')[0], self.XYZ.split('/')[1], self.XYZ.split('/')[2])
                if E == None:
                    print('\nNo energy found at ' + io.db_sp_path(self.XYZ.split('/')[0], self.XYZ.split('/')[1], self.XYZ.split('/')[2]                                                                                                , None, smiles,self.XYZ.split('/')[0], self.XYZ.split('/')[1], self.XYZ.split('/')[2]) + smiles + '.ene')
            else:
                E = '0'
                optim += '\n\t' + E 
            io.write_file(optim, '../output/' + self.typemol + str(n) + '_opt.out')
        elif '1' in self.xyzstart:
            optim = 'opt level1 0'
            for meas in measure:
                optim += '\n\t' + meas[1]                
            for i in range(len(angles)):
                optim += '\n\t{}'.format(deletedangles[i])              
            io.write_file(optim, '../output/' + self.typemol + str(n) + '_opt.out')

        return zmatstring

def build_molpro(meth,freqcalc,opt):
    """
    Builds molpro theory file 
    """
    method, basis = meth[2].split('/')
    nondft = ['hf', 'ccsd', 'cisd']
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


def build_theory(meths,nTS):
    """
    Builds theory.dat 
    meth[0] is module, meth[1] is program, and meth[2] is theory/basis
    """
    theory    = ''
    tsopt  = ' opt=(ts,calcfc,noeig,intern,maxcyc=50)\n '
    rpopt  = ' opt=internal\n '
    vwopt  = ' opt(internal,calcall) scf=qc\n '
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
        else:
            print meth[0] + ' is not a recognized program.\n'
            
    theory += 'End'

    return theory 

def build_estoktp(params, jobs, nreacs, nprods, nTS):
    """
    Builds estoktp.dat
    """
    stoich    = params[0]
    reactype  = params[1]
    coresh    = params[2]
    coresl    = params[3]
    
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
    eststring += ' 200MW  300MW\n gmemll gmemhl\n'
 
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

def update_interns(atoms, measure ,angles):
    """
    Converts internal coordinate information from test_chem to the form for
    a Z-Matrix required to run EStokTP
    """
    if len(atoms) == 0:
        return atoms, measure, angles 
    for index,atom in enumerate(atoms):
        atoms[index][0] = atom[0].lower() + str(index+1)
        if len(atom) > 1:
            atoms[index][1] = atoms[int(atom[1])-1][0]
            if len(atom) > 3:
                atoms[index][3] = atoms[int(atom[3])-1][0]
                if len(atom) > 5:
                    atoms[index][5] = atoms[int(atom[5])-1][0]
    return atoms, measure, angles

def find_period(zmat,hin):
    """
    Rough way of determining internal symmetry number (Hydrogen counting)
    """
    sym1 = sym2 = k = 0
    for i, row in enumerate(zmat[3:]):
        if hin.upper() in row[6]:
            k = i

    atom1 = zmat[k][1]
    atom2 = zmat[k][3]

    for row in zmat[1:]:
        if atom1 in row[1]:
            if 'h' in row[0]:
                sym1 += 1
        elif atom2 in row[1]:
            if 'h' in row[0]:
                sym2 += 1
    if sym1 == 1 or sym2 ==1:
        return max(sym1, sym2)
    elif sym1 == sym2:
        return sym1 + sym2
    else:
        return sym1 * sym2

def tau_hind_str(atoms, angles, interval, nsteps, mdtau):
    #TAU
    string  = '\nntau number of sampled coordinates\n'
    string += str(len(angles)) + '\n'
    string += ' -->nametau, taumin, taumax\n'
    for angle in angles:
        periodicity = find_period(atoms, angle)
        string += angle + ' 0 ' + interval + '\n'

    #1 and 2D HIND
    string += '\nnhind\n'
    string += str(len(angles)) + '\n'
    string += ' -->namehind,hindmin,hindmax,nhindsteps,period\n'
    for hin in angles:
        periodicity = find_period(atoms, hin)
        string += hin + ' 0 ' + str(float(interval)/periodicity)  + ' ' + str(int(round(float(nsteps)/periodicity))) + ' ' + str(periodicity)  + '\n'   

    if mdtau:
        mdtau   = mdtau.strip('D').strip('d')
        string += '\nnhind' + mdtau + 'D\n'
        string += '1\n'
        string += ' -->namehind,hindmin,hindmax,nhindsteps,period\n'
        for i in range(int(mdtau)):
            periodicity = find_period(atoms, angles[i])
            string += angles[i] + ' 0 ' + str(float(interval)/periodicity) + ' ' + str(int(round(float(nsteps)/periodicity))) + ' ' + str(periodicity)  + '\n'   
    return string


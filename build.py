#!/usr/bin/python

import os
import numpy as np
import sys

sys.path.insert(0, '/home/elliott/Packages/QTC/')
import iotools as io

class MOL:
    def __init__(self,opts,typemol = 'reac'):
       
        self.convert    = io.get_path('/home/elliott/Packages/TorsScan/test_chem') 
        self.typemol    = typemol

        #OPTIONS##############################
        self.nsamps     = opts[0]    #number of MonteCarlo sampling points
        self.interval   = opts[1]    #Interval for tors scan (should be same geometry as 0 degree)
        self.nsteps     = opts[2]    #Number of points to take on PES
        self.QTC        = opts[3]    #if QTC, openbabel, and pybel are not importable, self.QTC = 'False' 
                                     #be sure to have cartesian <SMILE string>.xyz in data/ (ex. C-C-C.xyz)
        ######################################
 

    def build_cart(self,smiles):
        """
        Uses QTC interface by Murat Keceli to Openbabel to generate cartesian coorinate file based 
        on SMILE string
        """
        import obtools as ob
        
        mol = ob.get_mol(smiles)

        self.charge = ob.get_charge(mol)
        self.mult   = ob.get_multiplicity(mol)
        stoich = ob.get_formula(mol)

        lines    =  ob.get_xyz(mol).split('\n')
        lines[0] = 'Geometry ' + lines[0] + ' Angstrom'
        del lines[1]

        io.write_file('\n'.join(lines), smiles + '.xyz')

        return stoich

    def read_cart(self,smiles): 
        """
        Runs Test_Chem by Yuri Georgievski on a file of Cartesian coordinates and collects
        the internal coordinate and torsional angle information that it outputs
        """
        import re

        #initialize
        atoms   = []
        measure = []
        angles  = []

        if self.QTC.lower() == 'true':
            stoich = self.build_cart(smiles)
        else:     #If we can't use QTC to build_cart we assume an xyz is given and hope default charge/mult work
            self.charge = 0
            self.mult   = 1
            stoich = smiles

        #Peform Test_Chem
        tempfile = 'temp'
        cart       = os.getcwd() + '/' +  smiles + '.xyz'
        os.system(self.convert + ' ' + cart + ' > ' + tempfile)

        if os.stat(tempfile).st_size == 0:
            print('failed')
            print('Please check that directory name and cartesian coordinate file name are equivalent')
            print('Please check that test_chem is in location: ' +  self.convert)
            return atoms, measure, angles
        
        #Get relevant data from Test_Chem output file
        props,lines = io.read_file(tempfile).split('Z-Matrix:\n')
        io.rm(tempfile)

        lines = lines.split('\n')
        for i,line in enumerate(lines):
            if line == '':
                break
            atoms.extend([line.rstrip('\n').split(',')])                        #connectivity

        for j in range(i+1,len(lines)):
            if 'Rot' in lines[j]:
                break
            measure.extend(lines[j].replace('=',' ').rstrip('\n').split()) #bond lengths, bond angles, and dihedral angles

        if lines[j].split(':')[1].rstrip('\n').strip() != '':
            angles  = lines[j].split(':')[1].strip().split(',')                #angles to scan

        if re.search("molecule is (\w+)", props).groups()[0] == 'nonlinear':        #Linearity
            self.ilin = ' 0'
        else:
            self.ilin = ' 1'
           
        self.symnum = re.search("symmetry number = (\w+)", props).groups()[0] #Symmetry factor
            
        measure = np.array(measure)                     
        measure = measure.reshape( len(measure)/2, 2)      #Puts measurements into two columns

        return stoich, atoms, measure, angles
    
    def update_interns(self,smiles):
        """
        Converts internal coordinate information from Test_Chem to the form for
        a Z-Matrix required to run EStokTP
        """
        stoich, atoms, measure, angles = self.read_cart(smiles)             #Converts ZMAT to EStokTP format

        if atoms == 0:
            return  

        for index,atom in enumerate(atoms):
            atoms[index][0] = atoms[index][0].lower() + str(index+1)
            if len(atom) > 1:
                atoms[index][1] = atoms[int(atoms[index][1])-1][0]
                if len(atom) > 3:
                    atoms[index][3] = atoms[int(atoms[index][3])-1][0]
                    if len(atom) > 5:
                        atoms[index][5] = atoms[int(atoms[index][5])-1][0]

        return stoich, atoms, measure, angles

    def find_period(self,zmat,hin):
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
        return max(sym1,sym2)

    def build(self,smiles,n):
        """ 
        Builds reacn.dat or prodn.dat for EStokTP withh user defined nosmps (Monte Carlo sampling points
        for geometry search) and nhindsteps (number of points on the PES) 
        """

        stoich, atoms, measure, angles = self.update_interns(smiles)
        #Stochastic Geometry Search############
        zmatstring  = 'nosmp dthresh ethresh\n'     
        zmatstring += self.nsamps + '  1.0  0.00001\n'


        #Torsional Scan Parameters#############
        zmatstring += '\nntau number of sampled coordinates\n'
        zmatstring += str(len(angles)) + '\n'
        zmatstring += ' -->nametau, taumin, taumax\n'
        for angle in angles:
            periodicity = self.find_period(atoms, angle)
            zmatstring += angle + ' 0 ' + str(int(self.interval)/periodicity) + '\n'

        hind = angles

        zmatstring += '\nnhind\n'
        zmatstring += str(len(angles)) + '\n'
        zmatstring += ' -->namehind,hindmin,hindmax,nhindsteps,period\n'
        for hin in hind:
            periodicity = self.find_period(atoms, hin)
            zmatstring += hin + ' 0 ' + str(int(self.interval)/periodicity) + ' ' + self.nsteps + ' ' + str(periodicity)  + '\n'   
            for i,meas in enumerate(measure):
                if hin.lower().strip() == meas[0].lower().strip(): #
                    measure = np.array([np.delete(measure.T[0],i),np.delete(measure.T[1],i)]).T
        #Size and linearity of molecule###########
        zmatstring += '\nnatom natomt ilin\n'
        zmatstring += str(len(atoms)) + ' ' + str(len(atoms)) + self.ilin + '\n'

        #Typical Z-Matrix########################
        zmatstring += '\ncharge  spin  atomlabel\n'
        zmatstring += str(self.charge) + ' ' + str(self.mult) + '\n'

        for row in atoms:
            for j in range(len(row)):
                zmatstring += row[j] + ' '
            zmatstring +='\n'
    
        zmatstring += '\nintcoor'

        for meas in measure:
            zmatstring += '\n' + meas[0] + ' ' + meas[1]

        #Sym factor and no. of electronic states#
        zmatstring += '\n\nSymmetryFactor\n' + self.symnum + '\n'
        zmatstring += '\nnelec\n1\n 0.  1.\n\nend\n'

        zmat       = self.typemol +str(n) + '.dat'
        io.write_file(zmatstring, zmat)

        return stoich

class THEORY:
    def __init__(self,meths):
        self.meths = meths
 
    def build(self):
        """
        Builds theory.dat 
        """
        meths    = self.meths

        theory  = ''
        for meth in meths:
            theory += meth[0] + ' ' + meth[1] + '\n '
            theory += meth[2] + ' opt=internal\n'
            theory += ' int=ultrafine nosym ' 
            if meth[0] == 'level1':
                theory += ' freq'
            theory += '\n\n'
        theory += 'End'

        io.write_file(theory, 'theory.dat')

        return 

class ESTOKTP:
    def __init__(self,stoich,jobs,opts,nreacs,nprods):

        self.stoich  = stoich
        self.jobs    = jobs
        self.nreacs  = nreacs
        self.nprods  = nprods
        self.coresh  = opts[0]
        self.coresl  = opts[1]
                 
    def build(self,reactype,nTS):
       """
       Builds esktoktp.dat
       """

       eststring  = ' Stoichiometry\t' + self.stoich.upper()
     
       PossibleRxns = {'addition','abstraction','isomerization','betascission','well',''}
       if reactype.lower() not in PossibleRxns:
           print('ReactionType ' + reactype +' is unrecognized, please use: Addition, Abstraction, Isomerization, or Betascission')
       elif reactype != '' and reactype.lower() != 'well':
           eststring  += '\n ReactionType\t' + reactype
           if nTS != 0:
               eststring += ' ' + str(nTS) + 'TS'

       eststring +='\n Debug  2'

       for job in self.jobs:
           for n in range(self.nreacs):
               if job == 'Opt_1':
                   eststring += '\n ' + 'Opt_Reac' + str(n+1) + '_1'
               else:
                   eststring += '\n ' +job + '_Reac' + str(n+1)
           for n in range(self.nprods):
               if job == 'Opt_1':
                   eststring += '\n ' + 'Opt_Prod' + str(n+1) + '_1'
               else:
                   eststring += '\n ' +job + '_Prod' + str(n+1)

       eststring += '\nEnd'
       eststring += '\n ' + self.coresh + ',' + self.coresl + '\n numprocll,numprochl\n'
       eststring += ' 200MW  300MW\n gmemll gmemhl\n'
       io.write_file(eststring, 'estoktp.dat')

       return 


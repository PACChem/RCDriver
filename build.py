#!/usr/bin/python

import os
import numpy as np
import sys

sys.path.insert(0, '/home/elliott/Packages/QTC/')
import iotools as io

class REAC:
    def __init__(self,smiles,opts):
       
        self.QTC = opts[3]#if QTC, openbabel, and pybel are not importable, self.QTC = 'False' 
                          #be sure to have cartesian <SMILE string>.xyz in data/ (ex. C-C-C.xyz)
      
        #OPTIONS##############################
        self.nsamps     = opts[0]    #number of MonteCarlo sampling points
        self.interval   = opts[1]    #Interval for tors scan (should be same geometry as 0 degree)
        self.nsteps     = opts[2]    #Number of points to take on PES
        ######################################
 
        self.smiles     = smiles  #Name of molecule
        self.cart       = os.getcwd() + '/' +  smiles + '.xyz'
        self.convert    = io.get_path('/home/elliott/Packages/TorsScan/test_chem') #self.convert    = io.get_path('test_chem')
        self.zmat       = 'reac1.dat'

    def build_cart(self):
        """
        Uses QTC interface by Murat Keceli to Openbabel to generate cartesian coorinate file based 
        on SMILE string
        """
        import obtools as ob
        
        mol = ob.get_mol(self.smiles)
        self.charge = ob.get_charge(mol)
        self.mult   = ob.get_multiplicity(mol)
        self.stoich = ob.get_formula(mol)

        lines =  ob.get_xyz(mol).split('\n')
        lines[0] = 'Geometry ' + lines[0] + ' Angstrom'
        del lines[1]
        io.write_file('\n'.join(lines), self.smiles + '.xyz')

        return
    def get_stoich(self):
        return self.stoich

    def read_cart(self): 
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
            self.build_cart()
        else: #we assume an xyz is given and hope default charge/mult work
            self.charge = 0
            self.mult   = 1
            self.stoich = self.smiles

        #open cartesian input file
        tempfile = 'temp'
        os.system(self.convert + ' ' + self.cart + ' > ' + tempfile)

        if os.stat(tempfile).st_size == 0:
            print('failed')
            print('Please check that directory name and cartesian coordinate file name are equivalent')
            print('Please check that test_chem is in location: ' +  self.convert)
            return atoms, measure, angles
        
        props,lines = io.read_file(tempfile).split('Z-Matrix:\n')
        lines = lines.split('\n')
        for i,line in enumerate(lines):
            #collect the connectivity once Z-Matrix is begun
            if line == '':
                break
            atoms.extend([line.rstrip('\n').split(',')])
        for j in range(i+1,len(lines)):
            #collect the bond lengths, bond angles, and dihedral angles
            if 'Rot' in lines[j]:
                break
            measure.extend(lines[j].replace('=',' ').rstrip('\n').split())
        #record angles to scan
        if lines[j].split(':')[1].rstrip('\n').strip() != '':
            angles  = lines[j].split(':')[1].strip().split(',')

        #Linearity
        if re.search("molecule is (\w+)", props).groups()[0] == 'nonlinear':
            self.ilin = ' 0'
        else:
            self.ilin = ' 1'
           
        #Symmetry factor
        self.symnum = re.search("symmetry number = (\w+)", props).groups()[0]
            
        ##Beta-scission, may use eventually
        #self.beta = re.search("Beta-scission bonds: (\w+)", line).groups()[0]

        measure = np.array(measure)
        measure = measure.reshape( len(measure)/2, 2)

        return atoms, measure, angles
    
    def update_interns(self):
        """
        Converts internal coordinate information from Test_Chem to the form for
        a Z-Matrix required to run EStokTP
        """
        #Converts ZMAT into the format required for EStokTP
        atoms, measure, angles = self.read_cart()

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

        return atoms, measure, angles

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

    def build(self):
        """ 
        Builds reac1.dat for EStokTP withh user defined nosmps (Monte Carlo sampling points
        for geometry search) and nhindsteps (number of points on the PES) 
        """

        atoms, measure, angles = self.update_interns()
        #Stochastic Geometry Search############
        zmatstring  = 'nosmp dthresh ethresh\n'     
        zmatstring += self.nsamps + '  1.0  0.00001\n'


        #Torsional Scan Parameters#############
        zmatstring += '\nntau number of sampled coordinates\n'
        zmatstring += str(len(angles)) + '\n'
        zmatstring += ' -->nametau, taumin, taumax\n'
        for angle in angles:
            periodicity = self.find_period(atoms, angle)
            zmatstring += angle + ' 0 ' + str(self.interval/periodicity) + '\n'

        hind = angles

        zmatstring += '\nnhind\n'
        zmatstring += str(len(angles)) + '\n'
        zmatstring += ' -->namehind,hindmin,hindmax,nhindsteps,period\n'
        for hin in hind:
            periodicity = self.find_period(atoms, hin)
            zmatstring += hin + ' 0 ' + str(self.interval/periodicity) + ' ' + self.nsteps + ' ' + str(periodicity)  + '\n'   
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

        io.write_file(zmatstring, self.zmat)

        return 

class THEORY:
    def __init__(self,meth,jobs):
        start = 0
        self.meth = meth
        self.jobs = jobs
        self.oth  = ('','freq','','')
 
    def build(self):
        """
        Builds theory.dat 
        """

        meth    = self.meth
        jobs    = self.jobs
        theory  = ''

        for i,job in enumerate(jobs):
            if meth[i][1] != '':
                theory += job + ' ' + meth[i][0] + '\n '
                theory += self.meth[i][1] + ' opt=internal\n'
                theory += ' int=ultrafine nosym ' + self.oth[i]+'\n\n'

        theory += 'End'

        io.write_file(theory, 'theory.dat')

        return 

class ESTOKTP:
    def __init__(self,stoich,methods):
        self.stoich  = stoich
        self.methods = methods
                 
    def build(self):
       """
       Builds esktoktp.dat
       """

       jobs  = ('Opt_Reac1','Opt_Reac1_1','1dTau_Reac1','HL_Reac1','Symm_reac1','kTP')
       eststring  = ' Stoichiometry\t' + self.stoich.upper()
       eststring +='\n Debug  2'
       for i,meth in enumerate(self.methods):
           if meth[1] != '':
               eststring += '\n ' + jobs[i]
       eststring += '\nEnd'
       eststring += '\n 10,6\n numprocll,numprochl\n'
       eststring += ' 200MW  300MW\n gmemll gmemhl\n'
       io.write_file(eststring, 'estoktp.dat')

       return 


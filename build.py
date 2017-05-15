#!/usr/bin/python

import os
import re
import numpy as np

class REAC:
    def __init__(self,stoich):
        
        #OPTIONS##############################
        self.stoic      = stoich #Name of molecule
        self.charge     = 0      #Charge
        self.mult       = 1      #Multiplicity
        self.interval   = '360'  #Interval for torsional scan (should be same geometry as 0 degree)
        self.nsteps     = '10'   #Number of points to take on PES
        self.nsamps     = '0'    #number of MonteCarlo sampling points
        ######################################

        self.cart       = os.getcwd() + '/../' + stoich + '.inp'
        self.convert    = '~/Tscan/TorsScan/test_chem'
        self.tempfile   = 'tempfile'
        self.zmat       = 'reac1.dat'


    def read_cart(self): 
       """
       Runs Test_Chem by Yuri Georgievski on a file of Cartesian coordinates and collects
       the internal coordinate and torsional angle information that it outputs
       """
       #open cartesian input file
       os.system(self.convert + ' ' + self.cart + ' > ' + self.tempfile)
       file = open(self.tempfile,'r')
       
       #initialize
       atoms   = []
       measure = []
       angles  = 0
       collect = 0

       for line in file:

           #exits loop if the converter failed
           if 'terminate' in line:
               return 0,0,0

           #collect the bond lengths, bond angles, and dihedral angles
           if collect == 2:
               
               #record angles to scan
               if 'Rot' in line:
                   angles  = line.split(':')[1].rstrip('\n').split(',')
                   #Reach end of ZMAT, stop collection
                   break
               measure.extend(line.replace('=',' ').rstrip('\n').split())

           #collect the connectivity once Z-Matrix is begun
           if collect == 1:
               if line == '\n':
                   collect = 2
               else:
                   atoms.extend([line.rstrip('\n').split(',')])

           #find beginning of zmatrix
           if collect == 0:
               if 'Z-Mat' in line:
                   collect = 1

               #Linearity
               if "molecule is" in line:
                   self.ilin =  re.search("molecule is (\w+)", line).groups()[0]
               
               #Symmetry factor
               if "symmetry number" in line:
                   self.symnum = re.search("symmetry number = (\w+)", line).groups()[0]
               
               #Symmetry factor
               if "Beta-scission" in line:
                   self.beta = re.search("Beta-scission bonds: (\w+)", line).groups()[0]

       file.close
       os.remove(self.tempfile)

       if self.ilin == 'nonlinear':
           self.ilin = ' 0'
       else:
           self.ilin = ' 1'
       
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

    def build_zmat(self):
       """ 
       Builds reac1.dat for EStokTP withh user defined nosmps (Monte Carlo sampling points
       for geometry search), interval (degrees scanned for torsional scan), nhindsteps (number of
       points on the PES) 
       """
       zmat  = open(self.zmat,'w')
       atoms, measure, angles = self.update_interns()
         
       #Stochastic Geometry Search############
       zmat.write('nosmp dthresh ethresh\n')      
       zmat.write(self.nsamps + '  1.0  0.00001\n')


       #Torsional Scan Parameters#############
       zmat.write('\nntau number of sampled coordinates\n')
       zmat.write(str(len(angles)) + '\n')
       zmat.write(' -->nametau, taumin, taumax\n') 
       for angle in angles:
           zmat.write(angle + ' 0 ' + self.interval + '\n')       

       hind = angles
       periodicity = '3'

       zmat.write('\nnhind\n')
       zmat.write(str(len(angles)) + '\n')
       zmat.write(' -->namehind,hindmin,hindmax,nhindsteps,period\n') 
       for hin in hind:
           zmat.write(hin + ' 0 ' + self.interval + ' ' + self.nsteps + ' ' + periodicity  + '\n')     
           for i,meas in enumerate(measure):
               if hin.lower().strip() == meas[0].lower().strip(): #
                   measure = np.array([np.delete(measure.T[0],i),np.delete(measure.T[1],i)]).T
       #Size and linearity of molecule###########
       zmat.write('\nnatom natomt ilin\n')
       zmat.write(str(len(atoms)) + ' ' + str(len(atoms)) + self.ilin + '\n')


       #Typical Z-Matrix########################
       zmat.write('\ncharge  spin  atomlabel\n')
       zmat.write(str(self.charge) + ' ' + str(self.mult) + '\n')

       for row in atoms:
           for j in range(len(row)):
               zmat.write(row[j] + ' ')
           zmat.write('\n')
    
       zmat.write('\nintcoor')

       for meas in measure:
           zmat.write('\n' + meas[0] + ' ' + meas[1])

       #Sym factor and no. of electronic states#
       zmat.write('\n\nSymmetryFactor\n' + self.symnum + '\n')
       zmat.write('\nnelec\n1\n 0.  1.\n\nend\n')


       zmat.close

       return 

class THEORY:
    def __init__(self,prog,meth):
       start = 0
       self.prog = prog
       self.meth = meth
       self.oth  = ('','freq','','')
 
    def build_theory(self):
       """
       Builds theory.dat 
       """
       theory  = open('theory.dat','w')
   
       level = (' level0', ' level1', ' hind_rotor',' symmetry')
       for lev in range(4):
           theory.write(level[lev] + ' ' + self.prog[lev] + '\n ')
           theory.write(self.meth[lev] + ' opt=internal\n')
           theory.write(' int=ultrafine nosym ' + self.oth[lev]+'\n\n')

       theory.write(' hlevel ' + self.prog[4] + '\n\n')
       theory.write('End')

       theory.close
       return 

class ESTOKTP:
    def __init__(self,stoich,jobs):
        self.stoich = stoich
        self.jobs   = jobs
                 
    def build_estoktp(self):
       """
       Builds esktoktp.dat
       """
       est  = open('estoktp.dat','w')
       est.write(' Stoichiometry\t' + self.stoich.upper())
       est.write('\n Debug  2')
       for job in self.jobs:
           est.write('\n ' + job)
       est.write('\nEnd')
       est.write('\n 12,4\n numprocll,numprochl\n')
       est.write(' 200MW  300MW\n gmemll gmemhl\n')
       est.close
       return 


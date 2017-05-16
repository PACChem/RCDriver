#!/usr/bin/python

import os
import sys
import re

sys.path.insert(0, './bin')

import build

class Tscan:
    def __init__(self,node):
      #  ####USER INPUTS###############################
      #  #stoichiometry       
      #  self.mol  = 'C-C-O'                   #smiles string
      #  #self.mol  = os.getcwd().split('/')[-1]  #based on directory name
      #  #program for torsional scan at level 1
      #  prog  = 'g09'
      #  #method  for torsional scan at level 1
      #  meth  = 'b3lyp/6-311+g(d,p)'
      #  #list the things you want EStokTP to do:
      #  #self.jobs  = ('Opt_React1','Opt_React1_1','1dTau_Reac1','kTP')
      #  self.jobs  = ('Opt_React1','1dTau_Reac1')
      #  ############################################

        #self.prog = ('g09',prog,'g09','g09','molpro')
        self.node = node
        self.get_options()
        #self.meth = ('b3lyp/6-31+g(d,p)', meth,'b3lyp/6-31+g(d,p)','b3lyp/6-31+g(d,p)')

    def get_options(self):

        options = open('input.dat','r')
        self.meth = [[],[],[],[],[]]
        for opts in options:
            if   'Molecule name' in opts:
                self.mol = opts.split(':')[1].strip('\n').strip()
            elif 'sampling point'in opts:
                self.nsamps   = opts.split(':')[1].strip('\n').strip()
            elif 'interval'      in opts:
                self.interval = int(opts.split(':')[1].strip('\n'))
            elif 'steps'         in opts:
                self.nsteps   = opts.split(':')[1].strip('\n').strip()
            elif 'Level 0'       in opts:
                self.meth[0] = [opts.split(':')[1].strip(),opts.split(':')[2].strip('\n').strip()]
            elif 'Level 1'       in opts:
                self.meth[1] = [opts.split(':')[1].strip(),opts.split(':')[2].strip('\n').strip()]
            elif 'Hind_Rotor'    in opts:
                self.meth[2] = [opts.split(':')[1].strip(),opts.split(':')[2].strip('\n').strip()]
            elif 'Symmetry'      in opts:
                self.meth[3] = [opts.split(':')[1].strip(),opts.split(':')[2].strip('\n').strip()]
            elif 'High Level'    in opts:
                self.meth[4] = [opts.split(':')[1].strip(),opts.split(':')[2].strip('\n').strip()]
        options.close
        return
       
    def build_subdirs(self):
        
        """
        Builds data and output subdirectories
        """
        print('Task: Building directories...')
        if not os.path.exists('./data'):
            os.makedirs('./data') 
        if not os.path.exists('./output'):
            os.makedirs('./output') 
        print('completed')
        return

    def build_files(self):

        """
        Runs the build functions for reac1.dat, theory.dat, and estoktp.dat
        """
        os.chdir('./data')

        print('Task: Building reac1.dat...')
        opts = (self.nsamps, self.interval,self.nsteps)
        reac = build.REAC(self.mol,opts)
        reac.build()
        print('completed')

        print('Task: Building theory.dat...')
        theory = build.THEORY(self.meth)
        theory.build()
        print('completed')

        print('Task: Building estoktp.dat...')
        estoktp = build.ESTOKTP(reac.get_stoich(),self.meth)
        estoktp.build()
        print('completed')

        os.chdir('..')

        return
  
    def execute(self):
        
        """
        Runs EStokTP
        """
        os.system('~/bin/run_estoktp.com ' + self.node)
        
        return

run = Tscan('b431')
run.build_subdirs()
run.build_files()
run.execute()

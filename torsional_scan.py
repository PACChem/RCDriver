#!/soft/python/2.6.4/bin/python

import os
import sys

sys.path.insert(0, './bin')

import build

class Tscan:
    def __init__(self,node):

        ####USER INPUTS###############################
        #stoichiometry       
        #self.mol  = 'c3h6ooh'
        self.mol  = os.getcwd().split('/')[-1]
        #program for torsional scan at level 1
        prog  = 'g09'
        #method  for torsional scan at level 1
        meth  = 'b3lyp/6-311+g(d,p)'
        #list the things you want EStokTP to do:
        #self.jobs  = ('Opt_React1','Opt_React1_1','1dTau_Reac1','kTP')
        self.jobs  = ('Opt_React1','1dTau_Reac1')
        ############################################

        #self.prog = ('g09',prog,'g09','g09','molpro')
        self.prog = ('g09',prog)
        self.meth = ('b3lyp/6-31+g(d,p)', meth)
        self.node = node
        #self.meth = ('b3lyp/6-31+g(d,p)', meth,'b3lyp/6-31+g(d,p)','b3lyp/6-31+g(d,p)')
       
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
        run = build.REAC(self.mol)
        run.build_zmat()
        print('completed')

        print('Task: Building theory.dat...')
        run = build.THEORY(self.prog,self.meth)
        run.build_theory()
        print('completed')

        print('Task: Building estoktp.dat...')
        run = build.ESTOKTP(self.mol,self.jobs)
        run.build_estoktp()
        print('completed')

        os.chdir('..')

        return
  
    def execute(self):
        
        """
        Runs EStokTP
        """
        os.system('~/bin/run_estoktp.com ' + self.node)
        
        return
\
#run = Tscan('b431')
#run.build_subdirs()
#run.build_files()
#run.execute()

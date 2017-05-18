#!/usr/bin/python

import os
import sys
import re

sys.path.insert(0, './bin')

import build

class Tscan:
    def __init__(self):
        ####DEFAULT INPUTS###############################
        self.mol  = 'C-C-C' #SMILE string
        self.QTC  = 'True'  #QTC, Openbabel, and Pybel installed
        self.node = 'debug' #Default node to run on in is debug (won't run)
        self.nsamps   = '5' #Number of MC sampling points
        self.interval = 360 #Interval to scan
        self.nsteps   = '4' #Number of steps on PES
        self.meth = [['g09', 'b3lyp/6-31+g(d,p)'], ['', ''], ['g09', 'b3lyp/6-31+g(d,p)'], ['', ''], ['', '']]
        self.jobs = ['level0', 'level1', 'hind_rotor', 'symmetry', 'hlevel']
        ############################################
        self.get_options()
        #self.jobs  = ('Opt_React1','Opt_React1_1','1dTau_Reac1','kTP')
       
    def get_options(self):

        options = open('input.dat','r')
        for opts in options:

            if   'node'      in opts:
                self.node     = opts.split(':')[1].strip('\n').strip()
            elif 'Use QTC'   in opts:
                self.QTC      = opts.split(':')[1].strip('\n').strip()
            elif 'Molecule'  in opts:
                self.mol      = opts.split(':')[1].strip('\n').strip()
            elif 'sampling'  in opts:
                self.nsamps   = opts.split(':')[1].strip('\n').strip()
            elif 'interval'  in opts:
                self.interval = int(opts.split(':')[1].strip('\n'))
            elif 'steps'     in opts:
                self.nstepsm  = opts.split(':')[1].strip('\n').strip()

            elif 'level0'    in opts:
                self.meth[0] = [opts.split(':')[1].strip(),opts.split(':')[2].strip('\n').strip()]
                self.jobs[0] = opts.split(':')[0].strip()
            elif 'level1'    in opts:
                self.meth[1] = [opts.split(':')[1].strip(),opts.split(':')[2].strip('\n').strip()]
                self.jobs[1] = opts.split(':')[0].strip()
            elif 'hind_rotor'in opts:
                self.meth[2] = [opts.split(':')[1].strip(),opts.split(':')[2].strip('\n').strip()]
                self.jobs[2] = opts.split(':')[0].strip()
            elif 'symmetry'  in opts:
                self.meth[3] = [opts.split(':')[1].strip(),opts.split(':')[2].strip('\n').strip()]
                self.jobs[3] = opts.split(':')[0].strip()
            elif 'hlevel'    in opts:
                self.meth[4] = [opts.split(':')[1].strip(),opts.split(':')[2].strip('\n').strip()]
                self.jobs[4] = opts.split(':')[0].strip()

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
        Requirements: QTC, Openbabel, Pybel
        """
        if self.QTC.lower() == 'skip': 
            print('task skipped')
            return

        os.chdir('./data')

        print('Task: Building reac1.dat...')
        opts = (self.nsamps, self.interval,self.nsteps,self.QTC)
        reac = build.REAC(self.mol,opts)
        reac.build()
        print('completed')

        print('Task: Building theory.dat...')
        theory = build.THEORY(self.meth,self.jobs)
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
        Runs EStokTP on a given blues node (default b431)
        use 0 in input file to run on login
        use d or debug to just make input files and not run EStokTP
        """
        print('Task: Submitting EStokTP job...')
        if self.node == '0':
            os.system('soft add +g09; soft add +gcc-5.3; ~/tscan/EStokTP/exe/estoktp.x > estotkp.log &')
        elif self.node == 'd' or self.node == 'debug':
            print('task skipped')
            return
        else:
            os.system('~/bin/run_estoktp.com ' + self.node)
        print('completed')
        return

    def extract(self):
        """
        Extracts necessary EStokTP output for MESS
        """
        print('Task: Extracting data for MESS...')
        if not os.path.exists('me_files'):
            print('failed -- me_files not found, check estoktp.log')
            return
        os.chdir('me_files')
        return
    
if __name__ == "__main__":

    run = Tscan()
    run.build_subdirs()
    run.build_files()
    run.execute()
    run.extract()

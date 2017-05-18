#!/usr/bin/python

import os
import sys
import build

sys.path.insert(0, './bin')
sys.path.insert(0, '/home/elliott/Packages/QTC/')
import iotools as io

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
        #clean this up when we decide what input should look like
        options = io.read_file('input.dat').split('\n')
        for opts in options:

            if   'node'      in opts:
                self.node     = self.get_val(opts,1)
            elif 'Use QTC'   in opts:
                self.QTC      = self.get_val(opts,1)
            elif 'Molecule'  in opts:
                self.mol      = self.get_val(opts,1)
            elif 'sampling'  in opts:
                self.nsamps   = self.get_val(opts,1)
            elif 'interval'  in opts:
                self.interval = self.get_val(opts,1)
            elif 'steps'     in opts:
                self.nsteps   = self.get_val(opts,1)

            elif 'level0'    in opts:
                self.meth[0] = [self.get_val(opts,1,'false'),self.get_val(opts,2)]
                self.jobs[0] =  self.get_val(opts,0)
            elif 'level1'    in opts:
                self.meth[1] = [self.get_val(opts,1,'false'),self.get_val(opts,2)]
                self.jobs[1] =  self.get_val(opts,0)
            elif 'hind_rotor'in opts:
                self.meth[2] = [self.get_val(opts,1,'false'),self.get_val(opts,2)]
                self.jobs[2] =  self.get_val(opts,0)
            elif 'symmetry'  in opts:
                self.meth[3] = [self.get_val(opts,1,'false'),self.get_val(opts,2)]
                self.jobs[3] =  self.get_val(opts,0)
            elif 'hlevel'    in opts:
                self.meth[4] = [self.get_val(opts,1,'false'),self.get_val(opts,2)]
                self.jobs[4] =  self.get_val(opts,0)

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

    def build_mess(self):
        """
        Compiles together the mess data extracted from EStokTP computations
        to form an mess input file
        """
        print('Task: Extracting data for MESS...')
        if not os.path.exists('me_files'):
            print('failed -- me_files not found, check estoktp.log')
            return

        #Copy first section of QTC pf file
        pfqtc = self.extract_mess(self.mol + '.pf').split('\n')
        pf   = ''
        for line in pfqtc:
            #if line == 'Frequencies': #use QTC geometry
            #   break
            if "Species" in line: #Use EStokTP geometry
               break
            pf += line + '\n'
        
        os.chdir('me_files')
        #EStokTP geometry
        ge = self.extract_mess('reac1_1dge.me')
        ge = " Species" +  ge.split("Species")[1]
        #Use EStokTP hr data
        hr = self.extract_mess('reac1_hr.me')
        #Use EStokTP frequencies (projected out scanned torsional modes)
        fr = self.extract_mess('reac1_fr.me')

        pf += ge + hr + fr
        os.chdir('..')
        io.write_file(pf,'new_' + self.mol + '.pf')

        return

    def extract_mess(self,filename):
        """
        Extracts necessary EStokTP output for MESS
        """
        lines = io.read_file(filename)
        if lines == '':
            print('failed -- ' + filename + ' is empty, check estoktp.log')
        return lines

    def get_val(self,opt,result,endl = 'True'):
        """
        Gets values we care about from a line of input file
        """
        if endl == 'True':
            return opt.split(':')[result].strip('\n').strip()
        else:
            return opt.split(':')[result].strip()
       


if __name__ == "__main__":

    run = Tscan()
    run.build_subdirs()
    run.build_files()
    run.execute()
    run.extract()

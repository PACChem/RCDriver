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
        ############################################

        self.get_options()
       
    def get_options(self):
        """
        Gets options from the input file
        """
        options = io.read_file('input.dat')

        self.set_theory(options)

        options      = options.split('\n')

        self.node    = self.set_param('node'    , options)
        self.QTC     = self.set_param('Use QTC' , options)
        self.mol     = self.set_param('Molecule', options)
        self.nsamps  = self.set_param('sampling', options)
        self.interval= self.set_param('interval', options)
        self.nsteps  = self.set_param('steps'   , options)
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
        theory = build.THEORY(self.meths)
        theory.build()
        print('completed')

        print('Task: Building estoktp.dat...')
        estoktp = build.ESTOKTP(reac.get_stoich(),self.jobs)
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
    def key_check(self,line,keyword):
        """                                                   
        Checks if a line of the input file has a keyword on it
        """                                                   
        if keyword in line:                                   
            return True                                       
        else:                                                 
            return False                                     
       
    def set_param(self,keyword,inputlines):
         """
         Sets parameter based on a keyword in inputfile
         """
         for line in inputlines:
             if self.key_check(line,keyword):
                 return  self.get_val(line,1)
         return
      
    def set_theory(self,inputlines):
        """
        Sets theory parameters
        """
        comps = {'Opt_Reac1':'level0','Opt_Reac1_1':'level1','1dTau_Reac1':'hind_rotor',
                           'HL_Reac1':'hind_rotor','Symm_Reac1':'symmetry','kTP':'hlevel'}

        inputlines = inputlines.replace(' ','')
        lines      = inputlines.split('--------------------------')[2].strip('-').split('\n')
        del lines[0]
        self.jobs  = []
        self.meths = []

        for line in lines:
            line = line.strip().split(':') 
            if self.key_check(comps,line[0]):
                if self.key_check(self.jobs,comps[line[0]]) and line[1] != '':
                    self.jobs.append(line[0])
                    self.meths.append([comps[line[0]],line[1],line[2]])
        return

 
if __name__ == "__main__":

    run = Tscan()
    run.build_subdirs()
    run.build_files()
    run.execute()
    run.extract()

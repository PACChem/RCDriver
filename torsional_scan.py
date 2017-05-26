#!/usr/bin/python

import os
import sys
import build
import numpy as np
sys.path.insert(0, './bin')
sys.path.insert(0, '/home/elliott/Packages/QTC/')
import iotools as io

class ES:
    def __init__(self):
    
        ####DEFAULT INPUTS#########################
        self.reacs    = 'CCC'   #list of SMILE strings of reactants
        self.prods    = ''      #list of SMILE strings of products
        self.reactype = ''      #type of reaction (default well)
        self.nTS      = 0       #Number of transition states (default 0)
        self.QTC      = 'True'  #QTC, Openbabel, and Pybel installed
        self.node     = 'debug' #Default node to run on in is debug (won't run)
        self.coresh   = '10'    #Default high number of cores is 10
        self.coresl   = '6'     #Default low number of cores is 10
        self.nsamps   = '5'     #Number of MC sampling points
        self.interval = 360     #Interval to scan
        self.nsteps   = '4'     #Number of steps on PES
        ###########################################

        self.get_options()      #Options from input file
       
    def get_options(self):
        """
        Gets options from the input file
        """
        options = io.read_file('input.dat')

        self.get_theory_params(options)

        options      = options.split('\n')

        self.reactype= self.get_param(self.reactype, 'Reaction type', options)
        self.nTS     = self.get_param(self.reactype, 'transition'   , options)
        self.reacs   = self.get_param(self.reacs   , 'Reactant'     , options)
        self.prods   = self.get_param(self.prods   , 'Product'      , options)
        self.node    = self.get_param(self.node    , 'node'         , options)
        self.coresh  = self.get_param(self.coresh  , 'cores high'   , options)
        self.coresl  = self.get_param(self.coresl  , 'cores low'    , options)
        self.QTC     = self.get_param(self.QTC     , 'Use QTC'      , options)
        self.nsamps  = self.get_param(self.nsamps  , 'sampling'     , options)
        self.interval= self.get_param(self.interval, 'interval'     , options)
        self.nsteps  = self.get_param(self.nsteps  , 'steps'        , options)

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
        Runs the build functions for reacn.dat, prodn.dat, theory.dat, and estoktp.dat
        Requirements: QTC, Openbabel, Pybel (OR prepared cartesian coordinate files) and Test_Chem
        """
        if self.QTC.lower() == 'skip': 
            print('Skipping build of data/*.dat')
            self.stoich = self.mol
            return

        os.chdir('./data')
        
        opts = (self.nsamps, self.interval,self.nsteps,self.QTC)
        stoich = []
        Reac = build.MOL(opts,'reac')
        Prod = build.MOL(opts,'prod')
        reacs = self.reacs.split(',')
        prods = self.prods.split(',')
        for i, reac in enumerate(reacs):
            print('Task: Building reac' + str(i+1) + '.dat...')
            stoich.append(Reac.build(reac.strip(),i+1))
            print('completed')
        for j, prod in enumerate(prods):
            print('Task: Building prod' + str(j+1) + '.dat...')
            stoich.append(Prod.build(prod.strip(),j+1))
            print('completed')

        self.stoich = stoich[0]
        self.mol    =   reac[0]

        print('Task: Building theory.dat...')
        theory = build.THEORY(self.meths)
        theory.build()
        print('completed')

        print('Task: Building estoktp.dat...')
        opts    = (self.coresh,self.coresl)
        estoktp = build.ESTOKTP(self.stoich,self.jobs,opts,i+1,j+1)
        estoktp.build(self.reactype,self.nTS)
        print('completed')

        os.chdir('..')

        return
  
    def execute(self):
        """
        Runs EStokTP on a given blues node (default b431)
        use 0 in input file to run on login
        use d or debug to just make input files and not run EStokTP
        Requirements: PACC member on blues
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
       
    def get_param(self,param,keyword,inputlines):
         """
         Sets parameter based on a keyword in inputfile
         """
         for line in inputlines:
             if self.key_check(line,keyword):
                 return  self.get_val(line,1)
         return param
      
    def get_theory_params(self,inputlines):
        """
        Sets theory parameters
        """
        comps = {'Opt':'level0','Opt_WellP':'level0','Opt_WellR':'level0','Grid_Opt_TS':'level0',
                           'Opt_TS_0':'level0_ts','TauO_TS':'level0_ts','Opt_1':'level1',
                           'Opts_TS_1':'level1_ts','1dTau':'hind_rotor','MdTau':'hind_rotor',
                           '1dTau_TS':'hind_rotor_ts','MdTau_TS':'hind_rotor_ts','Symm':'symmetry',
                           'Symm_TS':'symmetry_ts','HL':'hlevel',
                           'HL_TS':'hlevel_ts','IRC':'irc'}
                   #'Opt_WellP_1':'level1','Opt_WellR_1':'level1','1dTau_WellR':'hind_rotor', 
                   #'1dTau_WellP':'hind_rotor','MdTau_WellR':'hind_rotor','MdTau_WellP':'hind_rotor',
                   #'Symm_WellR':'symmetry','Symm_WellP':'symmetry','HL_WellP':'hlevel','HL_WellR':'hlevel',
        inputlines = inputlines.replace(' ','')
        lines      = inputlines.split('------------------------------')[2].strip('-').split('\n')
        del lines[0]
        self.jobs  = []
        self.meths = []
        templist   = []
        for line in lines:
            line = line.strip().split(':')
            if self.key_check(comps,line[0]) and line[1] != '':
                if self.key_check(templist,comps[line[0]]) ==  False:
                    self.meths.append([comps[line[0]],line[1],line[2]])
                    templist.append(comps[line[0]])
                self.jobs.append(line[0])
            elif line[0] != '':
                if line[1] != '':
                    print (line[0] + ' is not a recognized module')
        return


class MESS:
    def __init__(self,mol,stoich):
        self.mol      = mol
        self.stoich   = stoich

    def build(self):
        """
        Compiles together the mess data extracted from EStokTP computations
        to form an mess input file
        """
        if not os.path.exists('me_files'):
            print('failed -- me_files not found, check estoktp.log')
            return

        print('Task: Extracting data for MESS...')
        os.chdir('me_files')
        ge = self.extract_mess('reac1_1dge.me')                      #Copy EStokTP geometry
        ge = " Species " + self.stoich +  ge.split("Species")[1]
        hr = self.extract_mess('reac1_hr.me')                         #Copy EStokTP hr data
        fr = self.extract_mess('reac1_fr.me')                 #Copy EStokTP projfrequencies
        fr = fr.split('End')[0] + 'End  '
        os.chdir('..')                            # (projected out scanned torsional modes)
        print('completed')

        print('Task: Building MESS input file...')
        pfqtc = self.extract_mess(self.mol + '.pf').split('\n') #Copy 1st section of QTC pf file
        pf   = ''
        for line in pfqtc:
            #if line == 'Frequencies': #use QTC geometry
            #   break
            if "Species" in line:                                     #Use EStokTP geometry
               break
            pf += line + '\n'
        pf += ge + hr + fr
        io.write_file(pf,'new_' + self.mol + '.pf')
        print('completed')

        return

    def run(self):
        """
        Runs mess
        """
        print('Task: Running mess')
        os.system('/home/ygeorgi/build/crossrate/partition_function ' + 'new_' + self.mol + '.pf')
        print('completed')

    def extract_mess(self,filename):
        """
        Extracts necessary EStokTP output for MESS
        """
        lines = io.read_file(filename)
        if lines == '':
            print('failed -- ' + filename + ' is empty, check estoktp.log')
        return lines

 
if __name__ == "__main__":

    es = ES()
    es.build_subdirs()
    es.build_files()
    es.execute()
    mess = MESS(es.mol,es.stoich)
    mess.build()
    mess.run()

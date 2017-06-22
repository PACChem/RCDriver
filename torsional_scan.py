#!/usr/bin/python

import os
import sys
import build
import numpy as np
sys.path.insert(0, './bin')
sys.path.insert(0, '/home/elliott/Packages/QTC/')
import obtools as ob
import iotools as io
import tctools as tc
import rmg_reader as rg

class ES:
    def __init__(self,args):

        #for reac/prod/ts.dat
        self.restart =   args.restart 
        self.XYZ     =   args.XYZ      
        self.xyzstart=   args.xyzstart
        self.reacs   =   args.reacs    
        self.prods   =   args.prods    
        self.reactype=   args.reactype   
        self.nTS     =   args.nTS      
        self.TS      =   args.TS      
        #for theory.dat
        self.meths   =   args.meths
        #for estoktp.at
        self.jobs    =   args.jobs
        self.nsamps  =   args.nsamps   
        self.interval=   args.interval 
        self.nsteps  =   args.nsteps   
        self.node    =   args.node     
        self.coresh  =   args.coresh   
        self.coresl  =   args.coresl   
        
   
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
        os.chdir('./data')

        if self.restart < 10: 

            stoich = []

            params = (self.nsamps, self.interval,self.nsteps,self.XYZ,self.xyzstart,'MdTau' in self.jobs)
            Reac = build.MOL(params,'reac')
            Prod = build.MOL(params,'prod')
            
            params = (str(int(self.nsamps)*2), self.interval,self.nsteps,self.XYZ,'start','MdTau' in self.jobs)
            TS   = build.MOL(params,'ts') 

            reacs = self.reacs
            prods = self.prods
            
            if 'abstraction' in self.reactype.lower():
                key = ['[CH3]','[O]','[O][O]','O[O]','[OH]','[H]']
            elif  'addition' in self.reactype.lower():
                key = ['[O][O]']
            else: 
                key = []

            i,j,k = 0,0,0
            TScharge, TSspin = 0, 0

            for i, reac in enumerate(reacs,start=1):
                print('Task: Building reac' + str(i) + '.dat...')
                atoms, measure, angles = Reac.cart2zmat(reac)
                zmatstring = Reac.build(i, angles, atoms,  measure)
                zmat = Reac.typemol +str(i) + '.dat'
                io.write_file(zmatstring, zmat)
                if self.nTS > 0:
                    TScharge += Reac.charge
                    TSspin   += 1./2 * (Reac.mult - 1)
                    if i == 1:
                        TSangles, TSatoms  = angles, atoms
                    elif reac in key:
                        import shutil
                        shutil.copyfile('/home/elliott/Packages/TorsScan/abstractors/' + reac + '.dat','reac2.dat')
                        
                stoich.append(Reac.stoich)
                print('completed')

            for j, prod in enumerate(prods,start=1):
                print('Task: Building prod' + str(j) + '.dat...')
                atoms, measure, angles = Prod.cart2zmat(prod)
                zmatstring = Prod.build(j, angles, atoms, measure)
                zmat = Prod.typemol +str(j) + '.dat'
                io.write_file(zmatstring, zmat)
                stoich.append(Prod.stoich)
                print('completed')

            tstype = ['ts','wellr','wellp']
            for k in range(self.nTS):
                print('Task: Building ' + tstype[k] +  '.dat...')
                TS.charge = TScharge
                TS.mult   = int(2.*TSspin + 1)
                TS.symnum = ' 1'
                if k == 0:
                    zmatstring =TS.build(tstype[k], TSangles, TSatoms)
                else:
                    params = ('1', self.interval,self.nsteps,'False','start','MdTau' in self.jobs)
                    TS   = build.MOL(params,'ts') 
                    TS.charge = TScharge
                    TS.mult   = int(2.*TSspin + 1)
                    TS.symnum = ' 1'
                    zmatstring =TS.build(tstype[k], [], [])
                zmat = tstype[k] + '.dat'
                io.write_file(zmatstring, zmat)

            self.stoich = stoich[0]
            self.mol    =   reac[0]

        if self.restart < 7:

            print('Task: Building theory.dat...')
            theostring = build.build_theory(self.meths,self.nTS)
            io.write_file(theostring, 'theory.dat')
            print('completed')

            print('Task: Building estoktp.dat...')
            for l,job in enumerate(self.jobs):
                if job == 'Opt'   and self.restart > 2:
                    self.jobs[l]  = 'n' + job
                if job == 'Opt_1' and self.restart > 3:
                    self.jobs[l]  = 'n' + job
                if job == '1dTau' and self.restart > 4:
                    self.jobs[l]  = 'n' + job
                if job == 'MdTau' and self.restart > 5:
                    self.jobs[l]  = 'n' + job
            params    = (self.stoich, self.reactype, self.coresh,self.coresl,self)
            eststring = build.build_estoktp(params,self.jobs,i,j,self.nTS)
            io.write_file(eststring, 'estoktp.dat')
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
            os.system('soft add +g09; soft add +gcc-5.3; ~/Packages/EStokTP/exe/estoktp.x > estotkp.log &')
        elif self.node == 'd' or self.node == 'debug':
            print('task skipped')
            return
        else:
            os.system('~/bin/run_estoktp.com ' + self.node)
        print('completed')
        return
    

class MESS:
    def __init__(self):
        ladi = 'dida'
 
    def build(self,reacs,prods,anharm,anovrwrt,node,meths):
        """
        Compiles together the mess data extracted from EStokTP computations
        to form an mess input file
        """
        species = []
        speclist = []
        if not os.path.exists('me_files'):
            print('failed -- me_files not found, check estoktp.log')
            return
        if os.path.exists(reacs[0].strip() + '.pf'):
            pfqtc = self.extract_mess(reacs[0].strip() + '.pf').split('\n') #Copy 1st section of QTC pf file
            tf   = ''
            for line in pfqtc:
                #if line == 'Frequencies': #use QTC geometry
                #   break
                if "Species" in line:                                     #Use EStokTP geometry
                   break
                tf += line + '\n'
        else:
            tf =" Temperature(step[K],size)\t100.\t30\n RelativeTemperatureIncrement\t\t 0.001\n"
 

        for n,reac in enumerate(reacs):
            if reac == '':
                break
            species.append(reac)
            speclist.append('reac' + str(n+1))
            natom = ob.get_natom(reac)
            print('Task: Extracting MESS data for reac' + str(n+1) + '...')
            os.chdir('me_files')
            ge = self.extract_mess('reac' + str(n+1) + '_1dge.me')                      #Copy EStokTP geometry
            if 'Species' in ge:
                ge = " Species " + reac.strip() +  ge.split("Species")[1]
            elif 'Fragment' in ge:
                ge = " Fragment " + reac.strip() +  ge.split("Fragment")[1]
            ge,ge1 = ge.split('Core')
            hr = self.extract_mess('reac' + str(n+1) + '_hr.me')                         #Copy EStokTP hr data
            if not 'Core' in hr:
                ge = ge + 'Core' +ge1
            else:
                hr = hr.split('Quantum')
                hr = hr[0].rstrip() + '\n    Quantum'  + hr[1].lstrip()
            if anharm.lower() == 'false':
                fr = self.extract_mess('reac' + str(n+1) + '_fr.me')                 #Copy EStokTP projfrequencies
                fr = fr.split('End')[0] + 'End  '
            else:
                for meth in meths:
                    if str(anharm) in meth[0]:
                        anlevel = meth[2]
                os.chdir('..')                          
                fr = get_anharm('reac', str(n+1), natom,node,anlevel,anovrwrt,reac)  #(xmat with projected out scanned torsional modes)
                
                os.chdir('me_files')
            os.chdir('..')                          
            print('completed')

            print('Task: Building MESS input file...')
            pf = tf + ge + hr + fr
            io.write_file(pf + '\n',reac.strip() + '.pf')
            print('completed')

        for n, prod in enumerate(prods):
            if prod == '':
                break
            species.append(prod)
            speclist.append('prod' + str(n+1))
            print('Task: Extracting MESS data for prod' + str(n+1) + '...')
            os.chdir('me_files')
            ge = self.extract_mess('prod' + str(n+1) + '_1dge.me')                      #Copy EStokTP geometry
            if 'Species' in ge:
                ge = " Species " + prod.strip() +  ge.split("Species")[1]
            elif 'Fragment' in ge:
                ge = " Species " + prod.strip() +  ge.split("Fragment")[1]
            ge,ge1 = ge.split('Core')
            hr = self.extract_mess('prod' + str(n+1) + '_hr.me')                         #Copy EStokTP hr data
            if not 'Core' in hr:
                ge = ge + 'Core' +ge1
            else:
                hr = hr.split('Quantum')
                hr = hr[0].rstrip() + '\n    Quantum'  + hr[1].lstrip()
            if anharm.lower() == 'false':
                fr = self.extract_mess('prod' + str(n+1) + '_fr.me')                 #Copy EStokTP projfrequencies
                fr = fr.split('End')[0] + 'End  '
            else:
                for meth in meths:
                    if str(anharm) in meth[0]:
                        anlevel = meth[2]
                os.chdir('..')                           
                fr = get_anharm('prod', str(n+1), natom,node,anlevel,anovrwrt,prod)  #(xmat with projected out scanned torsional modes)
                os.chdir('me_files')
            os.chdir('..')                        
            print('completed')

            print('Task: Building MESS input file...')
            pf = tf + ge + hr + fr
            io.write_file(pf+'\n', prod.strip() + '.pf')
            print('completed')

        return species, speclist

    def run(self,reacs,prods,anharm,anovrwrt,node,meths):
        """
        Runs mess
        """
        import shutil
        import heatform as hf
        import re
            
        speciess,speclist = self.build(reacs,prods,anharm,anovrwrt,node,meths)

        for i,species in enumerate(speciess):
            deltaH = hf.main(ob.get_formula(species),'geoms/'+speclist[i] + '_l1.log')

            os.system('soft add +intel-16.0.0; soft add +gcc-5.3')
            print('Task: Running mess')
            tc.run_pf('/home/ygeorgi/build/crossrate/partition_function', species + '.pf')
            print('Generating thermp input.\n')
            
            stoich = ob.get_formula(ob.get_mol(species))
            inp = tc.get_thermp_input(stoich,deltaH)
            print('Running thermp.\n')
            os.rename(species + '.pf.log','pf.dat')
            tc.run_thermp(inp,'thermp.dat','pf.dat','/home/elliott/Packages/therm/thermp.exe')
            lines = io.read_file('thermp.out')
            deltaH298 = ' h298 final\s*([\d,\-,\.]*)'
            deltaH298 = re.findall(deltaH298,lines)[-1]
            print ('Running pac99.\n')
            shutil.copyfile('/home/elliott/Packages/therm/new.groups','./new.groups')
            shutil.copyfile(stoich + '.i97',species + '.i97')
            tc.run_pac99(species,'/home/elliott/Packages/therm/pac99.x')
            print('Converting to chemkin format.\n')
            chemkinfile = species + '.ckin'
            print('Writing chemking file {0}.\n'.format(chemkinfile))
            method = meths[-1][2]
            tc.write_chemkin_file(deltaH, method, species, chemkinfile)

            
            print('completed')
        return deltaH, deltaH298

    def extract_mess(self,filename):
        """
        Extracts necessary EStokTP output for MESS
        """
        lines = io.read_file(filename)
        if lines == '':
            print('failed -- ' + filename + ' is empty, check estoktp.log')
        return lines

class ARGS:
    def __init__(self,optionfile):  
        ####DEFAULT INPUTS#########################
        self.reacs    = 'CCC'   #list of SMILE strings of reactants
        self.prods    = ''      #list of SMILE strings of products
        self.TS       = ''      #list of SMILE strings of transition states
        self.reactype = ''      #type of reaction (default well)
        self.nTS      = '0'     #Number of transition states (default 0)
        self.XYZ      = 'True'  #Optimized XYZ provided
        self.xyzstart = 'start'  #Optimized XYZ provided
        self.node     = 'debug' #Default node to run on in is debug (won't run)
        self.coresh   = '10'    #Default high number of cores is 10
        self.coresl   = '6'     #Default low number of cores is 10
        self.nsamps   = '5'     #Number of MC sampling points
        self.interval = 360     #Interval to scan
        self.nsteps   = '4'     #Number of steps on PES
        self.rmg      = 'false' #RMG file to give input
        self.restart  = 'false' #Point at which to restart a computation
        self.anharm   = 'false' #Use and/or run anharmonic xmat computation
        self.anovrwrt = 'true' #Use and/or run anharmonic xmat computation
        self.alltherm = 'true' #Run all the thermochemistry scrips?
        ###########################################

        self.get_options(optionfile)      #Options from input file
        self.reacs = filter(None, self.reacs)  
        self.prods = filter(None, self.prods) 
 
    def get_theory_params(self,inputlines):
        """
        Sets theory parameters
        """
        comps = {'Opt':'level0','Opt_WellP':'level0','Opt_WellR':'level0','Grid_Opt_TS':'level0',
                           'Opt_TS_0':'level0_ts','TauO_TS':'level0_ts','Opt_1':'level1',
                           'Opts_TS_1':'level1_ts','1dTau':'hind_rotor','MdTau':'hind_rotor',
                           '1dTau_TS':'hind_rotor_ts','MdTau_TS':'hind_rotor_ts','Symm':'symmetry',
                           'Symm_TS':'symmetry_ts','HL':'hlevel',
                           'HL_TS':'hlevel_ts','Irc':'irc'}
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
            if key_check(comps,line[0]) and line[1] != '':
                if key_check(templist,comps[line[0]]) ==  False:
                    self.meths.append([comps[line[0]],line[1],line[2]])
                    templist.append(comps[line[0]])
                self.jobs.append(line[0])
            elif line[0] != '':
                if line[1] != '':
                    print (line[0] + ' is not a recognized module')
        return
    
    def get_options(self,optionfile):
        """
        Gets options from the input file
        """
        options = io.read_file(optionfile)
    
        self.get_theory_params(options)
    
        options      = options.split('\n')
      
        self.reactype= get_param(self.reactype, 'Reaction type', options)
        self.nTS     = int(get_param(self.nTS , 'of transition', options))
        self.reacs   = get_param(self.reacs   , 'Reactant'     , options).replace(' ','').split(',')
        self.prods   = get_param(self.prods   , 'Product'      , options).replace(' ','').split(',')
        self.TS      = get_param(self.TS      , 'Transition'   , options).replace(' ','').split(',')
        self.node    = get_param(self.node    , 'node'         , options)
        self.coresh  = get_param(self.coresh  , 'cores high'   , options)
        self.coresl  = get_param(self.coresl  , 'cores low'    , options)
        self.XYZ     = get_param(self.XYZ     , 'Use QTC'      , options)
        self.xyzstart= get_param(self.xyzstart, 'Use xyz as'   , options)
        self.nsamps  = get_param(self.nsamps  , 'sampling'     , options)
        self.interval= get_param(self.interval, 'interval'     , options)
        self.nsteps  = get_param(self.nsteps  , 'steps'        , options)
        self.restart = get_param(self.restart , 'Restart'      , options)
        self.anharm  = get_param(self.anharm  , 'Anharmonic'   , options)
        self.anovrwrt= get_param(self.anovrwrt, 'Overwrite an' , options)
        self.alltherm= get_param(self.alltherm, 'thermochemist', options)
    
        self.rmg     = get_param(self.rmg     , 'RMG input'    , options)
        if self.rmg.lower() != 'false' and self.rmg != '':
            self.rmg_params(self.rmg)
        if self.restart.lower() == 'false':
            self.restart = 0
        else:
            self.restart = int(self.restart)
    
        return
       
    def rmg_params(self,rmgfile):
    
        full = io.read_file(rmgfile)
        inputs = full.split('\r\n\r\n')
        dic ={}
        tsdic ={}
        for inp in inputs:
            if 'species' in inp:
                Spec = rg.SPECIES(inp)
                dic[Spec.label] = [Spec.smiles, Spec.mult]
            if 'transitionState' in inp:
                Trans = rg.TRANS(inp)
                tsdic[Trans.label] = [Trans.smiles, Trans.mult]
            if 'reaction' in inp:
                Reac = rg.REACTION(inp)
                self.reactype = Reac.reactype
                reactants     = Reac.reactants
                products      = Reac.products
                tss           = Reac.TS
                self.reacs = []
                self.prods = []
                self.TS    = []
                self.nTS   = Reac.nTS
                for reac in reactants:
                    if not reac in dic:
                        print 'incomplete RMG data'
                        break
                    else:
                        self.reacs.append(dic[reac][0])
                for prod in products:
                    if not prod in dic:
                        print 'incomplete RMG data'
                        break
                    else:
                        self.prods.append(dic[prod][0])
                for ts in tss:
                    if not ts in tsdic:
                        print 'incomplete RMG data'
                        break
                    else:
                        self.ts.append(tsdic[ts][0])
        return

def get_anharm(rorp,i,natom,node,anlevel,anovrwrt,species):

    import anharm

    opts= {}
    opts['node'      ] =  node
    opts['theory'    ] =  anlevel
    opts['logfile'   ] = 'geoms/' + rorp +  i + '_l1.log'
    opts['natoms'    ] =  natom 
    opts['freqfile'  ] = 'me_files/' + rorp +  i + '_fr.me' 
    opts['unprojfreq'] = 'me_files/' + rorp +  i + '_unpfr.me'
    if io.check_file(species + 'anharm.log') and not anovrwrt.lower() == 'true':
        opts['anharmlog' ] = species + 'anharm'
        opts['writegauss'] = 'false'
        opts['rungauss'  ] = 'false'
    else:
        opts['writegauss'] = 'true'
        opts['rungauss'  ] = 'true'
        opts['anharmlog' ] = species + 'anharm'

    return anharm.main(opts)

def get_val(opt,result,endl = 'True'):
    """
    Gets values we care about from a line of input file
    """
    if endl == 'True':
        return opt.split(':')[result].strip('\n').strip()
    else:
        return opt.split(':')[result].strip()

def key_check(line,keyword):
    """                                                   
    Checks if a line of the input file has a keyword on it
    """                                                   
    if keyword in line:                                   
        return True                                       
    else:                                                 
        return False                                    

def get_param(param,keyword,inputlines):
     """
     Sets parameter based on a keyword in inputfile
     """
     for line in inputlines:
         if key_check(line,keyword):
             return get_val(line,1)
     return param
 
if __name__ == "__main__":
   
    print("""\t\t   TORSSCAN
                    _,--._
                  ,'      `.
          |\     /          \     /|
          )o),/ ( ,--,  ,--, ) \.(o(
         /o/// /|            |\ \\ \\o\\
        / / |\ \(   .----,   )/ /| \ \\
        | | \o`-/    `--'    \-'o/ | |
        \ \  `,'              `.'  / /
     \.  \ `-'  ,'|   /\   |`.  `-' /  ,/
      \`. `.__,' /   /  \   \ `.__,' ,'/
       \o\     ,'  ,'    `.  `.     /o/
        \o`---'  ,'        `.  `---'o/
         `.____,'           `.____,'  """)


    args = ARGS('input.dat')
    es   = ES(args)             
    if es.restart < 7:
        es.build_subdirs()
        es.build_files()
        es.execute()
    import parse 
    if args.alltherm.lower() == 'true':
        mess = MESS()
        deltaH, deltaH298 = mess.run(args.reacs,args.prods,args.anharm,args.anovrwrt,args.node,args.meths)
    for i,reac in enumerate(args.reacs, start=1):
        lines = io.read_file('geoms/reac' + str(i) + '_l1.log')
        print('=====================\n          '+reac+'\n=====================')
        print 'Method: ' +      parse.gaussian_method(  lines)
        print 'Basis:  ' +      parse.gaussian_basisset(lines)
        print 'Energy: ' +  str(parse.gaussian_energy(  lines))
        print 'Zmatrix:' +      parse.gaussian_zmat(lines)
        print 'Heat of formation(  0K): ' + str(deltaH) + ' kcal /' + str(deltaH/.00038088/ 627.503) + ' kJ'
        print 'Heat of formation(298K): ' + deltaH298   + ' kcal /' + str(float(deltaH298)/.00038088/ 627.503) + ' kJ'
    for i,prod in enumerate(args.prods, start=1):
        lines = io.read_file('geoms/prod' + str(i) + '_l1.log')
        print('=====================\n          '+prod+'\n=====================')
        print 'Method: ' +      parse.gaussian_method(  lines)
        print 'Basis:  ' +      parse.gaussian_basisset(lines)
        print 'Energy: ' +  str(parse.gaussian_energy(  lines))
        print 'Zmatrix:' +      parse.gaussian_zmat(lines)
        print 'Heat of formation(  0K): ' + str(deltaH) + ' kcal /' + str(deltaH/.00038088/ 627.503) + ' kJ'
        print 'Heat of formation(298K): ' + deltaH298   + ' kcal /' + str(float(deltaH298)/.00038088/ 627.503) + ' kJ'
    if args.nTS > 0:
        lines = io.read_file('geoms/tsgta_l1.log')
        print('=====================\n        TS\n=====================')
        print 'Method: ' +      parse.gaussian_method(  lines)
        print 'Basis:  ' +      parse.gaussian_basisset(lines)
        print 'Energy: ' +  str(parse.gaussian_energy(  lines))
        print 'Zmatrix:' +      parse.gaussian_zmat(lines)
        print 'Heat of formation(  0K): ' + str(deltaH) + ' kcal /' + str(deltaH/.00038088/ 627.503) + ' kJ'
        print 'Heat of formation(298K): ' + deltaH298   + ' kcal /' + str(float(deltaH298)/.00038088/ 627.503) + ' kJ'

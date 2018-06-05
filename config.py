#!/usr/bin/python

import os

class CONFIG:
    def __init__(self,configfile):
        self.configfile = configfile
        
    def get_paths(self):
        """
        Read lines from configfile and return them as a list
        """
        paths = read_file(self.configfile)
        paths = paths.splitlines()

        return paths

    def path_dic(self):
        """
        Tranform lines of a configfile into a dictionary of paths
        """
        paths = self.get_paths()
        dic = {}
        for i in range(len(paths)):
            paths[i] = paths[i].split('=')
        for path in paths:
            if path[0].strip().startswith('#'):
                pass
            elif path[0].strip() == '':
                pass
            else:
                key = path[0].strip().lower()
                if len(path) > 1:
                    val = path[1].strip()
                    dic[key] = val
                else:
                    dic[key] = None

        return dic

class ARGS:
    def __init__(self,optionfile):  
        """
        Object to store input arguments from an option file
        """
        ####DEFAULT INPUTS#########################
        self.reacs    = 'CCC'   #list of SMILE strings of reactants
        self.prods    = ''      #list of SMILE strings of products
        self.reactype = ''      #type of reaction (default well)
        self.nTS      = '0'     #Number of transition states (default 0)

        self.restart  = 'false' #Point at which to restart a computation
        self.XYZ      = 'False' #Optimized XYZ provided
        self.xyzstart = 'start' #Optimized XYZ provided

        self.node     = 'debug' #Default node to run on in is debug (won't run)
        self.coresh   = '16'    #Default high number of cores is 10
        self.coresl   = '10'     #Default low number of cores is 10
        self.mem   = '200'     #Default low number of cores is 10

        self.optoptions   = 'internal'     #Guassian options
        self.nsamps   = ''     #Number of MC sampling points
        self.nrotor   = '0'      #Number of rotors
        self.abcd     = '3,1,3,100'      #ABCD params to calculate number of mc points
        self.interval = 360     #Interval to scan
        self.nsteps   = '4'     #Number of steps on PES
        self.mdtype   = '2'     #2 or 3D?

        self.anharm   = 'false' #Use and/or run anharmonic xmat computation
        self.anovrwrt = 'true' #Use and/or run anharmonic xmat computation
        self.alltherm = 'true' #Run all the thermochemistry scrips?
        self.qtchf    = 'false'#Enter precomputed heat of formation in a comma-seperated list
        self.hfbasis  = 'auto' #Specify basis for heat of formation?
        self.parseall = 'true' #Specify basis for heat of formation?
        self.rmg      = 'false' #RMG file to give input
        self.store    =  False   #RMG file to give input
        self.database = '/home/elliott/thermodb' #RMG file to give input
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
        inputlines = inputlines.replace(' ','')
        inputlines = inputlines.replace('gaussian','g09').replace('Gaussian','g09')
        lines      = inputlines.split('------------------------------')[2].strip('-').split('\n')
        del lines[0]
        self.jobs  = []
        self.meths = []
        templist   = []
        for line in lines:
            line = line.strip().split(':')
            if 'kTP' in line:
                self.jobs.append('kTP')
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
        options = read_file(optionfile)
    
        self.get_theory_params(options)
    
        options      = options.split('\n')
      
        self.reactype= get_param(self.reactype, 'Reaction type', options)
        self.nTS     = int(get_param(self.nTS , 'of transition', options))
        self.reacs   = get_param(self.reacs   , 'Reactant'     , options).replace(' ','').split(',')
        self.prods   = get_param(self.prods   , 'Product'      , options).replace(' ','').split(',')

        self.node    = get_param(self.node    , 'node'         , options)
        self.coresh  = get_param(self.coresh  , 'cores high'   , options)
        self.coresl  = get_param(self.coresl  , 'cores low'    , options)
        self.mem     = get_param(self.mem     , 'Memory'    , options)

        self.XYZ     = get_param(self.XYZ     , 'Use QTC'      , options)
        self.xyzstart= get_param(self.xyzstart, 'Use xyz as'   , options)

        self.optoptions  = get_param(self.optoptions  , 'Gaussian optim'     , options)
        self.nsamps  = get_param(self.nsamps  , 'sampling'     , options)
        self.nrotor  = get_param(self.nrotor  , 'Number of rotors'     , options)
        self.abcd    = get_param(self.abcd    , 'Calculate no. MC points'     , options)
        self.interval= get_param(self.interval, 'interval'     , options)
        self.nsteps  = get_param(self.nsteps  , 'steps'        , options)
        self.mdtype  = get_param(self.mdtype  , 'Multidim'     , options)

        self.restart = get_param(self.restart , 'Restart'      , options)

        self.anharm  = get_param(self.anharm  , 'Anharmonic'   , options)
        self.anovrwrt= get_param(self.anovrwrt, 'Overwrite an' , options)
        self.alltherm= get_param(self.alltherm, 'thermochemist', options)
        self.qtchf   = get_param(self.qtchf   , 'heat of formation', options).replace(' ','').split(',')
        self.hfbasis = get_param(self.hfbasis , 'Basis for hea', options)
        self.parseall= get_param(self.parseall, 'Parse all'    , options)

        self.rmg     = get_param(self.rmg     , 'RMG input'    , options)
        if self.rmg.lower() != 'false' and self.rmg != '':
            self.rmg_params(self.rmg)
        if self.restart.lower() == 'false':
            self.restart = 0
        else:
            self.restart = int(self.restart)
        if '1' in self.xyzstart and self.restart < 2:
            self.restart = 2
        if '0' in self.xyzstart and self.restart < 1:
            self.restart = 1
        return

    def rmg_params(self,rmgfile):
        """
        Reads the network style rmg output to build specieslist...
        this is likely not the input we will be using though, 
        so this function will likely never be used
        """ 
        full = read_file(rmgfile)
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

def read_file(filename):
    """
    Read a configfile, break program if there is not one
    """
    from os.path import isfile
    if isfile(filename): 
        with open(filename, 'r') as f:
            lines = f.read()
    else:
       print '{} not found'.format(filename)
       import sys
       sys.exit()
    return lines

def get_val(opt,result,endl = 'True'):
    """
    For ARGS, gets values we care about from a line of input file
    """
    if endl == 'True':
        return opt.split(':')[result].strip('\n').strip()
    else:
        return opt.split(':')[result].strip()

def key_check(line,keyword):
    """                                                   
    For ARGS, checks if a line of the input file has a keyword on it
    """                                                   
    if keyword in line:                                   
        return True                                       
    else:                                                 
        return False                                    

def get_param(param,keyword,inputlines):
     """
     For ARGS, sets parameter based on a keyword in inputfile
     """
     for line in inputlines:
         if key_check(line,keyword):
             if not line.startswith("#"):
                 return get_val(line,1)
     return param
 

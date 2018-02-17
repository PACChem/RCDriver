#!/home/keceli/anaconda2/bin/python

import os
import sys
import build
import numpy as np
sys.path.insert(0, './bin')
sys.path.insert(0, '/home/elliott/Packages/QTC')
import obtools as ob
import iotools as io
import tctools as tc
import rmg_reader as rg

class ES:
    def __init__(self,args):
        """
        Builds and runs EStokTP jobs based on args from ARGS class
        """
        #for reac/prod/ts.dat
        self.restart =   args.restart 
        self.XYZ     =   args.XYZ      
        self.xyzstart=   args.xyzstart
        self.reacs   =   args.reacs    
        self.prods   =   args.prods    
        self.reactype=   args.reactype   
        self.nTS     =   args.nTS      
        #for theory.dat
        self.meths   =   args.meths
        self.mdtype  =   args.mdtype
        self.optoptions=  args.optoptions
        #for estoktp.dat
        self.jobs    =   args.jobs
        self.nsamps  =   args.nsamps   
        self.abcd    =   args.abcd
        self.interval=   args.interval 
        self.nsteps  =   args.nsteps   
        self.node    =   args.node     
        self.coresh  =   args.coresh   
        self.coresl  =   args.coresl   
        self.mem     =   args.mem
        self.symnums =   [] 
   
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
        Requirements: QTC, Openbabel, Pybel (OR prepared cartesian coordinate files) and test_chem
        """
        os.chdir('./data')

        stoich = []

        #Create Read, Prod, and TS objects from parameters
        if not 'MdTau' in self.jobs: #and self.mdtype.lower() != 'auto':
            self.mdtype = ''
        params = (self.nsamps, self.abcd, self.interval,self.nsteps,self.XYZ,self.xyzstart,self.mdtype)
        Reac = build.MOL(params,'reac')
        Prod = build.MOL(params,'prod')
        params = (self.nsamps,self.abcd,self.interval,self.nsteps,self.XYZ,'start',self.mdtype)
        TS   = build.MOL(params,'ts') 

        reacs = self.reacs
        prods = self.prods
        
        #These keys will later replace generated reac file if there is a TS search
        if 'abstraction' in self.reactype.lower():
            key = ['[CH3]','[O]','[O][O]','O[O]','[OH]','[H]','O=O']
        elif  'addition' in self.reactype.lower():
            key = ['[O][O]']
        else: 
            key = []
        i,j,k = 0,0,0
        TScharge, TSspin = 0, 0

        #Build reacn.dat
        for i, reac in enumerate(reacs,start=1):
            print('Task: Building reac' + str(i) + '.dat...')
            atoms, measure, angles = Reac.cart2zmat(reac)
            self.reacs[i-1] = reac.split('_')[0]
            if self.mdtype.lower() == 'auto':
                Reac.MDTAU, self.jobs = prepare_mdtau(len(angles), self.jobs)
            zmatstring = Reac.build(i, reac, angles, atoms,  measure)
            zmat = Reac.typemol +str(i) + '.dat'
            io.write_file(zmatstring, zmat)
            if self.nTS > 0:
                if i == 1:
                    sort = Reac.sort
                TScharge += Reac.charge
                TSspin   += 1./2 * (Reac.mult - 1)
                if i == 1:
                    TSangles, TSatoms  = angles, atoms
                elif reac in key:
                    import shutil
                    shutil.copyfile('/home/elliott/Packages/TorsScan/abstractors/' + reac + '.dat','reac2.dat')
            self.nsamps = Reac.nsamps
            stoich.append(Reac.stoich)
            self.symnums.append(Reac.symnum)
       
            print('completed')

        #Build prodn.dat
        for j, prod in enumerate(prods,start=1):
            print('Task: Building prod' + str(j) + '.dat...')
            atoms, measure, angles = Prod.cart2zmat(prod)
            self.prods[j-1] = prod.split('_')[0]
            if self.mdtype.lower() == 'auto':
                Prod.MDTAU, self.jobs = prepare_mdtau(len(angles), self.jobs)
            zmatstring = Prod.build(j, prod, angles, atoms, measure)
            zmat = Prod.typemol +str(j) + '.dat'
            io.write_file(zmatstring, zmat)
            stoich.append(Prod.stoich)
            self.symnums.append(Prod.symnum)
            print('completed')
        
        #Build TS, wellr, and wellp.dat
        tstype = ['ts','wellr','wellp']
        TS.ijk  = Reac.ijk
        TS.sort = Reac.sort
        for k in range(self.nTS):
            print('Task: Building ' + tstype[k] +  '.dat...')
            TS.charge = TScharge
            TS.mult   = int(2.*TSspin + 1)
            TS.symnum = ' 1'
            if k == 0:
                zmatstring =TS.build(tstype[k], prod, TSangles, TSatoms)
            else:
                params = ('1', self.abcd,self.interval,self.nsteps,'False','start','MdTau' in self.jobs)
                TS   = build.MOL(params,'ts') 
                TS.charge = TScharge
                TS.mult   = int(2.*TSspin + 1)
                TS.symnum = ' 1'
                zmatstring =TS.build(tstype[k], [], [])
            zmat = tstype[k] + '.dat'
            io.write_file(zmatstring, zmat)

        self.stoich = stoich[0]
        self.mol    =   reac[0]
        #Builds theory.dat
        print('Task: Building theory.dat...')
        theostring = build.build_theory(self.meths,self.nTS,self.optoptions)
        io.write_file(theostring, 'theory.dat')
        print('completed')

        #Builds estoktp.dat to restart at any step
        print('Task: Building estoktp.dat...')
        self.jobs = update_jobs(self.jobs, self.restart)
        params    = (self.stoich, self.reactype, self.coresh,self.coresl,self.mem,self)
        eststring = build.build_estoktp(params,self.jobs,i,j,self.nTS)
        io.write_file(eststring, 'estoktp.dat')
        print('completed')

        os.chdir('..')

        return
 

    def execute(self):
        """
        Runs EStokTP on a given blues node (default debug)
        use 0 in input file to run on login
        use d or debug to just make input files and not run EStokTP
        Requirements: PACC member on blues
        """
        print('Task: Submitting EStokTP job...')
        if self.node == '0':
           # os.system('soft add +g09; soft add +gcc-5.3; /home/elliott/Packages/EStokTP/exe/estoktp.x >& estoktp.log')
            os.system('soft add +g09-e.01; soft add +gcc-5.3; /lcrc/project/PACC/codes/EStokTP/exe/estoktp.x >& estoktp.log')
        elif self.node == 'd' or self.node == 'debug':
            print('task skipped')
            return
        else:
            os.system('/home/elliott/bin/run_estoktp.com ' + self.node)
        print('completed')
        return
    

    def check_geoms(self,nsamps):
        """
        Checks MC geoms to make sure they are the same inchii as the starting species
        """
        print('Task: Checking level0 geometries...')
        n = 2
        filename =  'geoms/reac1_' + '1'.zfill(n) + '.xyz'
        lowfilename   = filename
        coords = io.read_file(filename)
        lowcoords = coords
        mol = ob.get_mol(coords)
        name =  ob.get_inchi_key(mol)
        energy = float(coords.split('\n')[1])
        for i in range(2, int(nsamps) + 1):
            filename =  'geoms/reac1_' + '{}'.format(i).zfill(n) + '.xyz'
            if io.check_file(filename):
                coords = io.read_file(filename)
                mol = ob.get_mol(coords)
                if name ==  ob.get_inchi_key(mol):
                    if float(coords.split('\n')[1]) < energy:
                       energy = float(coords.split('\n')[1]) 
                       lowcoords = coords
                       lowfilename   = filename
                else: 
                    print('Connectivity change after torsional optimization. (InChI mismatch) {}.')
        io.cp(lowfilename,'torsopt.xyz')
        #io.write_file("\n".join(lowcoords.split("\n")),'geom.xyz')
        io.write_file("\n".join(lowcoords.split("\n")[2:]),'geom.xyz')
        print('Monte Carlo sampling successfully found geom.xyz!')
        return 

class THERMO:
    def __init__(self,args):
        """
        Parses the EStokTP output for mess input files, runs heatform for the heat of formation, 
        thermp, pac99, and generates chemkin file
        """
        self.anfreqs  = []
        self.anxmat   = []
        self.reacs    = args.reacs
        self.prods    = args.prods    
        self.nTS      = args.nTS
        self.anharm   = args.anharm
        self.anovrwrt = args.anovrwrt
        self.node     = args.node
        self.meths    = args.meths
        self.hfbasis  = args.hfbasis
        self.hlen     = args.hlen
        self.enlevel  = args.enlevel.replace('g09','gaussian')
        self.optlevel = args.optlevel.replace('g09','gaussian')
        self.taulevel = args.taulevel
        self.qtchf    = args.qtchf
        self.symnums  = args.symnums

    def build(self):
        """
        Compiles together the mess data extracted from EStokTP computations
        to form an mess input file
        """
        reacs    = self.reacs
        prods    = self.prods    
        anharm   = self.anharm
        anovrwrt = self.anovrwrt
        node     = self.node
        meths    = self.meths
        symnums  = self.symnums
 
        if not self.taulevel:
            self.taulevel = self.optlevel
        optprog, optmethod, optbasis = self.optlevel.split('/')
        prog, method, basis = self.taulevel.split('/')
        optlevel = self.optlevel
        species  = []  #list of all smiles
        speclist = []  #list of reac1, reac2 etc..
        if not os.path.exists('me_files'):
            print('failed -- me_files not found, check estoktp.log')
            return
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
                ge = " Species " + reac.strip() +  ge.split("Fragment")[1]
            if 'Core' in ge:
                ge,ge1 = ge.split('Core')
            else:
                ge1 = ''
            hr = self.extract_mess('reac' + str(n+1) + '_hr.me')                         #Copy EStokTP hr data
            if not 'Core' in hr:
                ge = ge + 'Core' +ge1
            elif 'ge' != '':
                hr = hr.split('Quantum')
                hr = hr[0].rstrip() + '\n    Quantum'  + hr[1].lstrip()
            if len(symnums) > n:
                symnum = symnums[n]
            else:
                symnum = 1
            ge = ge.replace('SymmetryFactor    1.0000000000000','SymmetryFactor     ' + str(symnum) + '.0')
            if anharm.lower() == 'false':
                fr = self.extract_mess('reac' + str(n+1) + '_fr.me')                 #Copy EStokTP projfrequencies
                fr = fr.split('End')[0] + 'End  '
            else:
                if len(anharm.split('/')) > 3:
                    anharm = anharm.replace('gaussian','g09')
                    split = anharm.split('/')
                    optlevel = '{}/{}/{}'.format(split[0],split[1],split[2])
                    anlevel =  '{}/{}/{}'.format(split[3],split[4],split[5])
                elif len(anharm.split('/')) == 3:
                    anharm = anharm.replace('gaussian','g09')
                    split = anharm.split('/')
                    anlevel =  '{}/{}/{}'.format(split[0],split[1],split[2])
                else:
                    for meth in meths:
                        if str(anharm) in meth[0]:
                            anlevel = meth[1] + '/' +  meth[2]
                            break
                        else:
                            anlevel = ''
                os.chdir('..')                          
                anfr,fr1, anx,fr2,fr3 = get_anharm('reac', str(n+1), natom,node,anlevel,anovrwrt,reac, optlevel.split('/'))  #(xmat with projected out scanned torsional modes)
                fr =  fr1 +  fr2 + fr3
                self.anfreqs.append(anfr)
                self.anxmat.append(anx)
                os.chdir('me_files')
            zpve = io.read_file('reac' + str(n+1) + '_zpe.me')
            io.db_store_sp_prop(zpve, reac, 'zpve', prog = prog, optprog = optprog, method= method, optmethod=optmethod, basis=basis, optbasis=optbasis)  

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
            if len(symnums) > n+len(reacs):
                symnum = symnums[n+len(reacs)]
            else:
                symnum = 1
            ge = ge.replace('SymmetryFactor    1.0000000000000','SymmetryFactor     ' + str(symnum) + '.0')
            if anharm.lower() == 'false':
                fr = self.extract_mess('prod' + str(n+1) + '_fr.me')                 #Copy EStokTP projfrequencies
                fr = fr.split('End')[0] + 'End  '
            else:
                if len(anharm.split('/')) > 3:
                    anharm = anharm.replace('gaussian','g09')
                    split = anharm.split('/')
                    optlevel = '{}/{}/{}'.format(split[0],split[1],split[2])
                    anlevel =  '{}/{}/{}'.format(split[3],split[4],split[5])
                elif len(anharm.split('/')) == 3:
                    anharm = anharm.replace('gaussian','g09')
                    split = anharm.split('/')
                    anlevel =  '{}/{}/{}'.format(split[0],split[1],split[2])
                else:
                    for meth in meths:
                        if str(anharm) in meth[0]:
                            anlevel = meth[1] + '/' + meth[2]
                            break
                        else:
                            anlevel = ''
                os.chdir('..')                           
                anfr,fr1, anx,fr2, fr3 = get_anharm('prod', str(n+1), natom,node,anlevel,anovrwrt,prod, optlevel.split('/'))  #(xmat with projected out scanned torsional modes)
                fr = fr1 +  fr2 + fr3
                self.anfreqs.append(anfr)
                self.anxmat.append(anx)
                os.chdir('me_files')
            zpve = io.read_file('reac' + str(n+1) + '_zpe.me')
            io.db_store_sp_prop(zpve, reac, 'zpve', prog = prog, optprog = optprog, method= method, optmethod=optmethod, basis= basis, optbasis=optbasis)  

            os.chdir('..')                        
            print('completed')

            print('Task: Building MESS input file...')
            pf = tf + ge + hr + fr
            io.write_file(pf+'\n', prod.strip() + '.pf')
            print('completed')

        if self.nTS > 0:
            ts = reacs[0] + '_' + reac[1]
            species.append(ts)
            speclist.append('ts')
            print('Task: Extracting MESS data for TS...')
            os.chdir('me_files')
            ge = self.extract_mess('ts_1dge.me')                      #Copy EStokTP geometry
            if 'Species' in ge:
                ge = " Species " + prod.strip() +  ge.split("Species")[1]
            elif 'Fragment' in ge:
                ge = " Species " + prod.strip() +  ge.split("Fragment")[1]
            ge,ge1 = ge.split('Core')
            hr = self.extract_mess('ts_hr.me')                         #Copy EStokTP hr data
            if not 'Core' in hr:
                ge = ge + 'Core' +ge1
            else:
                hr = hr.split('Quantum')
                hr = hr[0].rstrip() + '\n    Quantum'  + hr[1].lstrip()
            fr = self.extract_mess('ts_fr.me')                 #Copy EStokTP projfrequencies
            fr = fr.split('End')[0] + 'End  '
            #if anharm.lower() == 'false':
            #    fr = self.extract_mess('ts_fr.me')                 #Copy EStokTP projfrequencies
            #    fr = fr.split('End')[0] + 'End  '
            #else:
            #    if '/' in anharm:
            #        anlevel = anharm
            #    else:
            #        for meth in meths:
            #            if str(anharm) in meth[0]:
            #                anlevel = meth[1] + '/' + meth[2]
            #            else:
            #                anlevel = ''
            #    os.chdir('..')                           
            #    anfr,fr1, anx,fr2, fr3 = get_anharm('prod', str(n+1), natom,node,anlevel,anovrwrt,prod, optlevel.split('/'))  #(xmat with projected out scanned torsional modes)
            #    fr = fr1 +  fr2 + fr3
            #    self.anfreqs.append(anfr)
            #    self.anxmat.append(anx)
            #    os.chdir('me_files')
            zpve = io.read_file('ts_zpe.me')
            #io.db_store_sp_prop(zpve, ts, 'zpve', prog = prog, optprog = optprog, method= method, optmethod=optmethod, basis= basis, optbasis=optbasis)  
            os.chdir('..')                        
            print('completed')

            print('Task: Building MESS input file...')
            pf = tf + ge + hr + fr
            io.write_file(pf+'\n', ts.strip() + '.pf')
            print('completed')

        return species, speclist

    def run(self):
        """
        Runs heatform, partition_function, thermp, pac99, and write chemkin file
        """
        import shutil
        import heatform as hf
        import re
            
        reacs    = self.reacs
        prods    = self.prods    
        anharm   = self.anharm
        anovrwrt = self.anovrwrt
        node     = self.node
        meths    = self.meths
        hfbasis  = self.hfbasis
        qtchf    = self.qtchf
        enlevel  = self.enlevel
        hlen     = self.hlen
        self.hfbases = []
        speciess,speclist = self.build()
        self.dH0   = []
        self.dH298 = []
        anharmbool = False
        if anharm.lower() != 'false':
            anharmbool = True
        for i,species in enumerate(speciess):
            if len(qtchf) >= i and qtchf[i].lower() not in ['false', 'auto']:
                deltaH = float(qtchf[i])
                hfbasis = ['N/A']
                self.hfbases.append(hfbasis)
            else:
                logfile = 'geoms/'+speclist[i] + '_l1.log'
                if speclist[i] ==  'ts':
                    logfile = 'geoms/'+speclist[i] + 'gta_l1.log'
                if io.check_file(logfile):
                    lines = io.read_file(logfile)
                    energy=pa.energy(lines)[1]
                    zpve  = pa.zpve(lines)
                    printE = '{}-    E: {:5g} pulled from: {}'.format(species, energy, logfile)
                    printzpve = '{}- zpve: {:5g} pulled from: {}'.format(species, zpve, logfile)
                    if enlevel != 'optlevel':
                        energy = hlen[i]
                        printE = '{}-    E: {:5g} pulled from: {}'.format(species, energy, 'me_files/'+speclist[i] + '_en.me')
                    if io.check_file('me_files/'+speclist[i] + '_zpe.me'):
                        zpve = float(io.read_file('me_files/'+speclist[i] + '_zpe.me'))
                        printzpve = '{}- ZPVE: {:5g} pulled from: {}'.format(species, zpve, 'me_files/'+speclist[i] + '_zpe.me')
                    if zpve:
                        energy = energy + zpve
                    print printE + '\n' +  printzpve
                    deltaH, hfbasis = hf.main(species,logfile,E=energy,basis=self.hfbasis,anharm=anharmbool,enlevel=enlevel)
                    self.hfbases.append(hfbasis)
                else:
                    deltaH = 0.00
            self.dH0.append(deltaH)
            if not speclist[i] == 'ts':
                print('Task: Running mess')
                tc.run_pf('/home/ygeorgi/build/crossrate/partition_function', species + '.pf')
                print('Generating thermp input.\n')
                
                stoich = ob.get_formula(ob.get_mol(species))
                inp = tc.get_thermp_input(stoich,deltaH)
                print('Running thermp.\n')
                if io.check_file(species+'.pf.log'):
                    os.rename(species + '.pf.log','pf.dat')
                else:
                    print ('ERROR: no pf.dat produces, try soft adding gcc-5.3 and intel-16.0.0 and use Restart at: 5!')
                tc.run_thermp(inp,'thermp.dat','pf.dat','/home/elliott/Packages/therm/thermp.exe')
                lines = io.read_file('thermp.out')
                deltaH298 = ' h298 final\s*([\d,\-,\.]*)'
                deltaH298 = re.findall(deltaH298,lines)[-1]
                self.dH298.append(deltaH298)
                print ('Running pac99.\n')
                shutil.copyfile('/home/elliott/Packages/therm/new.groups','./new.groups')
                shutil.copyfile(stoich + '.i97',species + '.i97')
                tc.run_pac99(species,'/home/elliott/Packages/therm/pac99.x')
                print('Converting to chemkin format.\n')
                chemkinfile = species + '.ckin'
                print('Writing chemkin file {0}.\n'.format(chemkinfile))
                method = meths[-1][2]
                tc.write_chemkin_file(deltaH, method, species, chemkinfile)

            
            print('completed')
        return 

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
        options = io.read_file(optionfile)
    
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

def update_jobs(jobs, restart):
    """
    Updates the job list based on restart level
    """
    for l,job in enumerate(jobs):
        if job == 'Opt'   and restart > 0:
            jobs[l]  = 'n' + job
        if job == 'Opt_1' and restart > 1:
            jobs[l]  = 'n' + job
        if job == '1dTau' and restart > 2:
            jobs[l]  = 'n' + job
        if job == 'MdTau' and restart > 3:
            jobs[l]  = 'n' + job
    return jobs

def prepare_mdtau(nrot, jobs):
    """ 
    Returns what mdtau should be set to based on the number of hindered rotors and 
    inserts MdTau into the joblist if need be
    """
    mdtau = None
    if nrot > 0  and '1dTau' in jobs:
        mdtau = '1'
        if nrot > 1:
            mdtau = '2'
            if nrot > 2:
                mdtau = '3'
        if 'MdTau' not in jobs:
            index = jobs.index('1dTau')
            jobs.insert(index+1, 'MdTau')
    return mdtau, jobs

def get_anharm(rorp,i,natom,node,anlevel,anovrwrt,species, optlevel):
    """
    Runs the anharm module to project out torsional modes from xmatrix and
    find the updated vpt2 frequencies
    """
    import anharm
    opts= {}
    opts['smiles'    ] =  species
    opts['node'      ] =  node
    opts['theory'    ] =  anlevel
    opts['optlevel'  ] =  optlevel
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
 
def parse(n, species, lines, optlevel,enlevel,hlen,lines2):
    """
    Parses and prints just the quantum chemistry output
    """
    printstr= '\n=====================\n          '+species+'\n=====================\n'
    if enlevel == 'optlevel':
        prog   =  pa.get_prog(lines) 
        method =  pa.method(lines).lower().lstrip('r')
        basis  =  pa.basisset(lines).lower()
        energy = str(pa.energy(lines)[1]) 
    else:
        prog, method, basis = enlevel.split('/') 
        if len(hlen) > n-1:
            energy = str(hlen[n-1])
        else:
            energy = 'N/A'
    optprog, optmethod, optbasis = optlevel.split('/')
    zmat   =  pa.zmat(lines)   
    if prog == 'g09': 
        xyz    =  pa.gaussian_geo(lines) 
    else: 
        xyz    = None
    rotcon = ', '.join(pa.rotconsts(lines))
    freqs  =  pa.freqs(lines)
    xmat   = []
    if lines2 != '':
        pfreqs  = pa.EStokTP_freqs(lines2)
    else: 
        pfreqs = []
    
    freqs   = ', '.join(freq for freq in freqs[::-1])
    pfreqs  = ', '.join('{:>9}'.format(freq) for freq in pfreqs)

    printstr += 'Optimized at : {}\n'.format(optlevel)
    printstr += 'Prog  : ' + prog   + '\n'
    printstr += 'Method: ' + method + '\n' 
    printstr += 'Basis:  ' + basis  + '\n'
    printstr += 'Energy: ' + energy + ' A.U.\n'
    printstr += 'Rotational Constants:' + rotcon  + ' GHz\n'
    printstr += 'Zmatrix (Angstrom):' + zmat   + '\n'
    if xyz != None:
        printstr += 'Cartesian coordinates (Angstrom):\n' + pa.xyz(lines) 
    printstr += '\nUnprojected Frequencies (cm-1):\t'  + freqs
    if lines2 != '':
        printstr += '\nProjected Frequencies   (cm-1):\t'  + pfreqs
    return printstr

#def parse_qtc(species, jobs):
#    if 'Opt' in jobs:
#        if io.check_file('geoms/reac1_01.xyz'):
#            lines = io.read_lines('goems/reac1_01.xyz')
#            io.write_file(lines, species + '.xyz')
#    if 'Opt_1' in jobs:
#        if io.check_file('geoms/reac1_l1.xyz'):
#            lines = io.read_lines('goems/reac1_l1.xyz')
#            io.write_file(lines, species + '.xyz')
#    if 'HL' in jobs:
#        if io.check_file('geoms/reac1_l1.xyz'):
#            lines = io.read_lines('goems/reac1_l1.xyz')
#            io.write_file(lines, species + '.xyz')
#    return
      

def parse_all(n, species, lines, anharm,anfreqs,anxmat,deltaH0,deltaH298,hfbasis,optlevel,enlevel,taulevel,hlen,lines2):
    """
    Parses, prints, and stores the quantum chemistry and thermochemistry output
    """
    printstr= '=====================\n          '+species+'\n=====================\n'
    if enlevel == 'optlevel':
        prog   =  pa.get_prog(lines) 
        method =  pa.method(lines).lower().lstrip('r')
        basis  =  pa.basisset(lines).lower() 
        energy = str(pa.energy(lines)[1]) 
    else:
        prog, method, basis = enlevel.split('/') 
        if len(hlen) >= n+1:
            energy = str(hlen[n])
        else:
            energy = 'N/A'
    optprog, optmethod, optbasis = optlevel.split('/')
    zmat   =  pa.zmat(lines)   
    if prog == 'g09': 
        xyz    =  pa.gaussian_geo(lines) 
    else: 
        xyz    = None
    rotcon = ', '.join(pa.rotconsts(lines))
    freqs  =  pa.freqs(lines)
    xmat   = []
    if lines2 != '':
        pfreqs  = pa.EStokTP_freqs(lines2)
    else: 
        pfreqs = []
    
    #d = {optlevel:
    #     {prog:
    #      {'torsscan':
    #       {method:
    #        {basis:{
    #         'number of basis functions':0,
    #         'energy':energy,
    #          'geometry':{
    #           'xyz':xyz,
    #           'harmonic frequencies' : freqs,
    #           'xmat': xmat }}}}}}}
    #print d
    freqs   = ', '.join(freq for freq in freqs[::-1])
    pfreqs  = ', '.join('{:>9}'.format(freq) for freq in pfreqs)

    printstr += 'Optimized at : {}\n'.format(optlevel)
    printstr += 'Energy: ' + energy + ' A.U.\n'
    printstr += 'Prog  : ' + prog   + '\n'
    printstr += 'Method: ' + method + '\n' 
    printstr += 'Basis:  ' + basis  + '\n'
    printstr += 'Energy: ' + energy + ' A.U.\n'
    printstr += 'Rotational Constants:' + rotcon  + ' GHz\n'
    printstr += 'Zmatrix (Angstrom):' + zmat   + '\n'

    io.db_store_opt_prop(zmat, species, 'zmat', None, optprog, optmethod, optbasis)
    if xyz != None:
        printstr += 'Cartesian coordinates (Angstrom):\n' + pa.xyz(lines) 
        io.db_store_opt_prop(xyz, species, 'geo', None, optprog, optmethod, optbasis)
        xyz = str(len(xyz.split('\n'))) + '\n\n' + xyz
        io.db_store_opt_prop(xyz, species, 'xyz', None, prog, optmethod, optbasis)
    printstr += '\nUnprojected Frequencies (cm-1):\t'  + freqs
    if lines2 != '':
        printstr += '\nProjected Frequencies   (cm-1):\t'  + pfreqs
        io.db_store_sp_prop(pfreqs,species,'phrm', None, prog, method, basis, optprog, optmethod, optbasis)
    io.db_store_sp_prop(energy, species, 'ene', None, prog, method, basis, optprog, optmethod, optbasis)
    io.db_store_sp_prop(freqs,  species, 'hrm', None, prog, method, basis, optprog, optmethod, optbasis)
    io.db_store_sp_prop(rotcon, species,  'rc', None, prog, method, basis, optprog, optmethod, optbasis)
    if anharm != 'false':
        anpfr     = ', '.join('%4.4f' % freq for freq in anfreqs[i-1])
        pxmat     =('\n\t'.join(['\t'.join(['%3.3f'%item for item in row]) 
                    for row in anxmat[n-1]]))
        io.db_store_sp_prop(anpfr, species, 'panhrm', None, prog, method, basis, optprog, optmethod, optbasis)
        #io.db_store_sp_prop(pxmat, species, 'pxmat', None, prog, method, basis, prog, method, basis)   ##probably not at right level of theory
        printstr += '\nAnharmonic Frequencies  (cm-1):\t' + anpfr
        printstr += '\nX matrix:\n\t' + pxmat #+   anxmat[i-1]
     
    printstr += '\nHeat of formation(  0K): ' + str(deltaH[i-1]) + ' kcal /' + str(deltaH[n-1]/.00038088/ 627.503) + ' kJ\n'
    hf298k    = str(deltaH298[n-1])
    if not io.check_file(io.db_sp_path(prog, method, basis, None, species, optprog, optmethod, optbasis) + '/' + species + '.hf298k'):
        io.db_store_sp_prop('Energy (kcal)\tBasis\n----------------------------------\n',species,'hf298k',None,prog,method,basis, optprog, optmethod, optbasis)
    if len(hfbasis) >= n+1:
        io.db_append_sp_prop(hf298k + '\t' + ', '.join(hfbasis[n]) + '\n', species, 'hf298k',None, prog,method,basis, optprog, optmethod, optbasis)
    printstr += 'Heat of formation(298K): ' + hf298k   + ' kcal /' + str(float(deltaH298[n-1])/.00038088/ 627.503) + ' kJ\n'
    return printstr

def ts_parse(n, lines):
    """
    Similar to parse, but only for TS
    """
    tstype = ['TS', 'WELLR', 'WELLP']
    printstr= '=====================\n          '+tstype[n]+'\n=====================\n'
    prog   =  pa.get_prog(lines) 
    method =  pa.method(  lines).lower().lstrip('r')
    basis  =  pa.basisset(lines).lower() 
    energy = str(pa.energy(lines)[1]) 
    zmat   =  pa.zmat(    lines)    
    xyz    =  pa.xyz(     lines) 
    rotcon = ', '.join(pa.rotconsts(lines))
    freqs  = ', '.join(freq for freq in pa.freqs(lines)[::-1])

    printstr += 'Prog  : ' + prog   + '\n'
    printstr += 'Method: ' + method + '\n' 
    printstr += 'Basis:  ' + basis  + '\n'
    printstr += 'Energy: ' + energy + ' A.U.\n'
    printstr += 'Rotational Constants:' + rotcon  + ' GHz\n'
    printstr += 'Zmatrix (Angstrom):' + zmat   + '\n'

    if xyz != None:
        printstr += 'Cartesian coordinates (Angstrom):\n' + pa.xyz(lines) 
    printstr += '\nUnprojected Frequencies (cm-1):\t'  + freqs + '\n'
    return printstr 

def printheader():
    printstr  = '\n\n ______________________'
    printstr += '\n||                    ||'
    printstr += '\n||       OUTPUT       ||'
    printstr += '\n||____________________||\n\n'
    return printstr

def random_cute_animal():
    import random 
    print(random.choice(["""\t\t   TORSSCAN
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
         `.____,'           `.____,'  """,

     """    
            ,,,         ,,,
          ;"   ^;     ;'   ",
         ;    s$$$$$$$s      ;
          ,  ss$$$$$$$$$$s  ,'
          ;s$$$$$$$$$$$$$$$
          $$$$$$$$$$$$$$$$$$
         $$$$P""Y$$$Y""W$$$$$
         $$$$  0"$$$"0  $$$$$
         $$$$  .$$$$$.  $$$$
          $$$$$$$$$$$$$$$$$
     TORS  "Y$$$"'*'"$$$Y"  SCAN    
              "$$b.d$$"        """,
    """
                         _.---~-~-~~-..
     ..       __.    .-~               ~-.
     ((\     /   `}.~      Try            `.
      \\\\\   {     }     Sarah    /     \   \\
  (\   \\\\~~       }      opts   |       }   \\
   \`.-~-^~     }  ,-,.    and  |       )    \\
   (___,    ) _}  (    :  scans |    / /      `.
    `----._-~.     _\ \ |_       \   / /- _     -.
           ~~----~~  \ \| ~~--~~~(  + /     ~-.   '--~.
                     /  /         \  \         `~-..__ `~__
                  __/  /          _\  )               `~~---'
                .<___.'         .<___/
    """]))
    return

if __name__ == "__main__":
  
   
    random_cute_animal()
     
    #Get arguments##########
    import sys
    if len(sys.argv) > 1:
        optionfile = sys.argv[-1]
    else:
        optionfile = 'input.dat'
    args = ARGS(optionfile)

    #Build and Run EStokTP##
    es   = ES(args)             
    if args.restart < 5:
        es.build_subdirs()
        es.build_files()
        es.execute()
        if ("Opt" in args.jobs and not "Opt_1" in args.jobs):
            es.check_geoms(es.nsamps)

    #HEY ME! move this outside of main
    optlevel = ''
    taulevel = ''
    if '1' in args.xyzstart:
        if len(args.XYZ.split('/')) >  2:
            optlevel = args.XYZ
    else:
        for meth in args.meths:
            if meth[0] == 'level1':
                optlevel = '{}/{}'.format(meth[1],meth[2])
            if 'tau' in meth[0]:
                taulevel = '{}/{}'.format(meth[1], meth[2].split('/')[1])
    args.optlevel = optlevel
    args.taulevel = taulevel
    args.enlevel  = 'optlevel'
    args.hlen     = []
    for i in range(len(args.reacs)):
        if io.check_file('me_files/reac'+str(i+1) + '_en.me'):
            args.hlen.append(float(io.read_file('me_files/reac'+str(i+1) + '_en.me')))
            for meth in args.meths:
                if 'hlevel' in meth:
                    args.enlevel = '{}/{}'.format(meth[1],meth[2])
    for i in range(len(args.prods)):
        if io.check_file('me_files/prod'+str(i+1) + '_en.me'):
            args.hlen.append(float(io.read_file('me_files/prod'+str(i+1) + '_en.me')))
            for meth in args.meths:
                if 'hlevel' in meth:
                    args.enlevel = '{}/{}'.format(meth[1],meth[2])


    #Builds and runs the thermochemistry files
    if args.alltherm.lower() == 'true':
        import patools as pa
        args.symnums = es.symnums
        thermo = THERMO(args)
        thermo.run()
        anfreqs = thermo.anfreqs
        anxmat  = thermo.anxmat
        deltaH  = thermo.dH0
        delH298 = thermo.dH298
  
        #HEY ME! move parsing outside of main
        if args.parseall.lower() == 'true':
            printstr  = printheader()
            for i,reac in enumerate(args.reacs, start=1):
                lines   = io.read_file('geoms/reac' + str(i) + '_l1.log')
                lines2  = ''
                if io.check_file('me_files/reac' +  str(i) + '_fr.me'):
                    lines2  = io.read_file('me_files/reac' +  str(i) + '_fr.me')
                printstr += parse_all(i, reac, lines, args.anharm,anfreqs,anxmat,deltaH,delH298,thermo.hfbases,optlevel,args.enlevel,taulevel,args.hlen,lines2)
            for j,prod in enumerate(args.prods, start=1):
                lines = io.read_file('geoms/prod' + str(j) + '_l1.log')
                lines2 = ''
                if io.check_file('me_files/reac' +  str(j) + '_fr.me'):
                    lines2  = io.read_file('me_files/prod' +  str(j) + '_fr.me')
                printstr += parse_all(i+j-1, prod, lines, args.anharm,anfreqs,anxmat,deltaH,delH298,thermo.hfbases,optlevel,args.enlevel,taulevel,args.hlen,lines2)
            if args.nTS > 0:
                lines = io.read_file('geoms/tsgta_l1.log')
                printstr += ts_parse(0,lines)
                if args.nTS > 1:
                    lines = io.read_file('geoms/wellr_l1.log')
                    printstr += ts_parse(1,lines)
                    if args.nTS > 2:
                        lines = io.read_file('geoms/wellp_l1.log')
                        printstr += ts_parse(2,lines)
            print printstr

    elif args.parseall.lower() == 'true':
        import patools as pa
        printstr  = printheader()
        for i,reac in enumerate(args.reacs, start=1):
            lines   = io.read_file('geoms/reac' + str(i) + '_l1.log')
            lines2  = ''
            if io.check_file('me_files/reac' +  str(i) + '_fr.me'):
                lines2  = io.read_file('me_files/reac' +  str(i) + '_fr.me')
            printstr += parse(i, reac, lines, optlevel,args.enlevel,args.hlen,lines2)
        for j,prod in enumerate(args.prods, start=1):
            lines = io.read_file('geoms/prod' + str(j) + '_l1.log')
            lines2 = ''
            if io.check_file('me_files/reac' +  str(j) + '_fr.me'):
                lines2  = io.read_file('me_files/prod' +  str(j) + '_fr.me')
            printstr += parse(i+j-1, prod, lines, optlevel,args.enlevel,args.hlen,lines2)
        print printstr
#    if args.parseall.lower() == 'qtc':
#        parse_qtc(args.reacs[0],args.jobs)

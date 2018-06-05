#!/home/keceli/anaconda2/bin/python

import os
import sys
import build
import numpy as np
import rmg_reader as rg

class ES:
    def __init__(self,args,paths):
        """
        Builds and runs EStokTP jobs based on args from ARGS class
        """
        #Paths to executables and libraries
        self.paths   =   paths
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

        if not 'MdTau' in self.jobs: #and self.mdtype.lower() != 'auto':
            self.mdtype = ''
   
    def build_files(self, msg=''):

        """
        Runs the build functions for reacn.dat, prodn.dat, theory.dat, and estoktp.dat
        Requirements: QTC, Openbabel, Pybel (OR prepared cartesian coordinate files) and x2z
        """
        os.chdir('./data')
        self.stoichs= []


        #Create Read, Prod, and TS objects from parameters
        params = (self.nsamps, self.abcd, self.interval,self.nsteps,self.XYZ,self.xyzstart,self.mdtype)
        Reac   = build.MOL(self.paths, params, 'reac')
        Prod   = build.MOL(self.paths, params, 'prod')
        params = (self.nsamps,self.abcd,self.interval,self.nsteps,self.XYZ,'start',self.mdtype)
        TS     = build.MOL(self.paths, params,   'ts') 

        reacs = self.reacs
        prods = self.prods
        
        key = self.set_keys()
        i,j,k = 0,0,0
        TSprops = [0,0,[],[]] #charge, spin, angles, atoms

        #Build reacn.dat
        for i, reac in enumerate(reacs,start=1):

            msg += 'Task: Building reac{:g}.dat...'.format(i)
            msg  = log_msg(msg)
            self.reacs, angles, atoms = self.build_moldat(Reac, reacs, i)
            TSprops = self.prep_reacs4TS(Reac, reac, key, angles, atoms, TSprops)
            self.nsamps = Reac.nsamps
            msg += 'completed'
            msg  = log_msg(msg)

        #Build prodn.dat
        for j, prod in enumerate(prods,start=1):
            msg += 'Task: Building prod{:g}.dat...'.format(j)
            msg  = log_msg(msg)
            self.prods, angles, atoms = self.build_moldat(Prod, prods, j)
            msg += 'completed'
            msg  = log_msg(msg)
        
        #Build TS, wellr, and wellp.dat
        tstype = ['ts','wellr','wellp']
        TS.ijk  = Reac.ijk
        TS.sort = Reac.sort
        for k in range(self.nTS):
            msg += 'Task: Building ' + tstype[k] +  '.dat...'
            msg  = log_msg(msg)
            TS.charge = TSprops[0]
            TS.mult   = int(2.*TSprops[1] + 1)
            TS.symnum = ' 1'
            if k == 0:
                zmatstring = TS.build(tstype[k], prod, TSprops[2], TSprops[3])
            else:
                params = ('1', self.abcd,self.interval,self.nsteps,'False','start','MdTau' in self.jobs)
                TS   = build.MOL(params,'ts') 
                TS.charge = TScharge
                TS.mult   = int(2.*TSspin + 1)
                TS.symnum = ' 1'
                zmatstring =TS.build(tstype[k], [], [])
            zmat = tstype[k] + '.dat'
            io.write_file(zmatstring, zmat)
            msg += 'completed'
            msg  = log_msg(msg)

        self.stoich = self.stoichs[0]
        self.mol    =   reac[0]
        #Builds theory.dat
        msg += 'Task: Building theory.dat...'
        msg  = log_msg(msg)
        theostring = build.build_theory(self.meths,self.nTS,self.optoptions)
        io.write_file(theostring, 'theory.dat')
        msg += 'completed'
        msg  = log_msg(msg)

        #Builds estoktp.dat to restart at any step
        msg += 'Task: Building estoktp.dat...'
        msg  = log_msg(msg)
        self.jobs = update_jobs(self.jobs, self.restart)
        params    = (self.stoich, self.reactype, self.coresh,self.coresl,self.mem,self)
        eststring = build.build_estoktp(params,self.jobs,i,j,self.nTS)
        io.write_file(eststring, 'estoktp.dat')
        msg += 'completed'
        msg  = log_msg(msg)

        os.chdir('..')

        return
 

    def execute(self, msg = ''):
        """
        Runs EStokTP on a given blues node (default debug)
        use 0 in input file to run on login
        use d or debug to just make input files and not run EStokTP
        Requirements: PACC member on blues
        """
        paths = self.paths
        g09     = get_paths(paths,  'g09')
        gcc     = get_paths(paths,  'gcc')
        intel   = get_paths(paths,  'intel')
        estoktp = get_paths(paths,'estoktp')
        msg += 'Task: Submitting EStokTP job...'
        msg  = log_msg(msg)
        if self.node == '0':
           # os.system('soft add +g09; soft add +gcc-5.3; /home/elliott/Packages/EStokTP/exe/estoktp.x >& estoktp.log')
            os.system('{0}; {1}; {2}; {3}  >& estoktp.log'.format(gcc, intel, g09, estoktp))
            msg += 'completed'
        elif self.node == 'd' or self.node == 'debug':
            msg += 'task skipped'
        else:
            ssh = get_paths(paths, 'ssh')
            os.system('exec {3} -n {4} "cd `pwd`;{0}; {1}; {2}; {5} >& estoktp.log"'.format(gcc, intel, g09, ssh, self.node, estoktp))
            msg += 'completed'
        msg  = log_msg(msg)
        return
    
    def build_subdirs(self,msg=''):
        
        """
        Builds data and output subdirectories
        """
        msg += 'Task: Building directories...'
        if not os.path.exists('./data'):
            os.makedirs('./data') 
        if not os.path.exists('./output'):
            os.makedirs('./output') 
        msg += 'completed'
        msg = log_msg(msg)
        return 

    def build_moldat(self, MOL, mollist, n):
        """
        Builds reac1.dat, reac2.dat, prod1.dat etc
        """
        mol = mollist[n-1]
        atoms, measure, angles  = MOL.cart2zmat(mol)
        mollist[n-1] = mol.split('_')[0]

        if self.mdtype.lower() == 'auto':
            MOL.MDTAU, self.jobs = prepare_mdtau(len(angles), self.jobs)

        zmatstring = MOL.build(n, mol, angles, atoms,  measure)
        zmat = MOL.typemol + str(n) + '.dat'
        io.write_file(zmatstring, zmat)
        self.stoichs.append(MOL.stoich)
        self.symnums.append(MOL.symnum)
        return mollist, angles, atoms

    def prep_reacs4TS(self, MOL, reac, key, angles, atoms, tsprops):
        """
        Sets transition state charge, mult, angles, and atoms respectively based 
        on the reactants.  May use abstractor templates if key matches
        """
        if self.nTS > 0:
            if i == 1:
                sort = MOL.sort
            tsprops[0] += MOL.charge
            tsprops[1] += 1./2 * (Reac.mult - 1)
            if i == 1:
                tsprops[2], tsprops[3]  = angles, atoms
            elif reac in key:
                import shutil
                shutil.copyfile(self.paths['torsscan'] + '/abstractors/' + reac + '.dat','reac2.dat')
            return tsprops

    

    def check_geoms(self,nsamps,msg=''):
        """
        Checks MC geoms to make sure they are the same inchii as the starting species
        """
        msg += 'Task: Checking level0 geometries...'
        msg  = log_msg(msg)
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
        msg += 'Monte Carlo sampling successfully found geom.xyz!'
        msg  = log_msg(msg)
        return 

    def me_file_abs_path(self):
        """
        Replaces relative path in mdhr file with absolute path
        """
        if io.check_file('me_files/reac1_hr.me'):
            lines = io.read_file('me_files/reac1_hr.me')
            if "PotentialEnergySurface[kcal/mol]" in lines:
                before, after = lines.split("PotentialEnergySurface[kcal/mol]")
                after = after.split('\n')
                after[0] = after[0].replace('./',io.pwd() + '/')
                lines = before + "PotentialEnergySurface[kcal/mol]" + '\n'.join(after)
                io.write_file(lines,'me_files/reac1_hr.me')
        return

    def set_keys(self):
        """
        These keys will later replace generated reac file if there is a TS search
        """
        if 'abstraction' in self.reactype.lower():
            key = ['[CH3]','[O]','[O][O]','O[O]','[OH]','[H]','O=O']
        elif  'addition' in self.reactype.lower():
            key = ['[O][O]']
        else: 
            key = []
        return key

class THERMO:
    def __init__(self,args):
        """
        Parses the EStokTP output for mess input files, runs heatform for the heat of formation, 
        thermp, pac99, and generates chemkin file
        """
        self.anfreqs  = []
        self.anxmat   = []
        self.hfbasis  = args.hfbasis
        self.hlen     = args.hlen

    def build_pfinput(self, args, msg = ''):
        """
        Compiles together the mess data extracted from EStokTP computations
        to form an mess input file
        """
        reacs    = args.reacs
        prods    = args.prods    
        anharm   = args.anharm
        anovrwrt = args.anovrwrt
        node     = args.node
        meths    = args.meths
        symnums  = args.symnums

        optlevel = args.optlevel.replace('g09','gaussian') 
        if not args.taulevel:
            taulevel = optlevel
        optprog, optmethod, optbasis = optlevel.split('/')
        prog, method, basis = taulevel.split('/')
        species  = []  #list of all smiles
        speclist = []  #list of reac1, reac2 etc..

        if not os.path.exists('me_files'):
            print('failed -- me_files not found, check estoktp.log')
            return

        tf = " Temperature(step[K],size)\t100.\t30\n RelativeTemperatureIncrement\t\t 0.001\n"
        zp = "   ZeroPointEnergy[1/cm] 0\n   InterpolationEnergyMax[kcal/mol] 500\n"
        for n,reac in enumerate(reacs):
           
            if reac == '':
                break
            msg += 'Task: Extracting MESS data for reac{:g}...'.format(n+1)
            msg = log_msg(msg)

            try:
                os.chdir('me_files')
                natom = ob.get_natom(reac)
                ge, hr = self.read_gehr(reac, 'reac',  n)
                if len(symnums) > n:
                    symnum = symnums[n]
                else:
                    symnum = 1
                ge = ge.replace('SymmetryFactor    1.0000000000000','SymmetryFactor     {:.2f}'.format(float(symnum)))
                fr = self.get_fr(reac, natom, 'reac', anharm, anovrwrt, meths, node, n, args.store)
                os.chdir('..')                              
                msg += 'completed'
                msg  = log_msg(msg)
                pf = tf + ge + zp +  fr + hr
                io.write_file(pf + '\n',reac.strip() + '.pf')
                msg += 'Task: Building MESS input file...'
                msg += 'completed'
                msg  = log_msg(msg)
                species.append(reac)
                speclist.append('reac' + str(n+1))
            except IOError:
                msg += 'me_files are missing, check me_file/*, estoktp.log, and output/estoktp.out' 
                msg  = log_msg(msg)
                os.chdir('..')                         

        for n, prod in enumerate(prods):

            if prod == '':
                break
            msg += 'Task: Extracting MESS data for reac{:g}...'.format(n+1)
            msg = log_msg(msg)
 
            try:    
                os.chdir('me_files')
                natom = ob.get_natom(prod)
                ge, hr = self.read_gehr(prod, 'prod', n)
                if len(symnums) > n+len(reacs):
                    symnum = symnums[n+len(reacs)]
                else:
                    symnum = 1
                ge = ge.replace('SymmetryFactor    1.0000000000000','SymmetryFactor        {:.2f}'.format(float(symnum)))
                fr = self.get_fr(prod, natom, 'prod', anharm, anovrwrt, meths, node, n, args.store)
                os.chdir('..')                        
                msg += 'completed'
                msg  = log_msg(msg)
                pf = tf + ge + zp + fr + hr
                io.write_file(pf+'\n', prod.strip() + '.pf')
                msg += 'Task: Building MESS input file...'
                msg += 'completed'
                msg  = log_msg(msg)
                species.append(prod)
                speclist.append('prod' + str(n+1))
            except IOError:
                msg += 'me_files are missing, check me_file/*, estoktp.log, and output/estoktp.out' 
                msg  = log_msg(msg)
                os.chdir('..')                         

        if args.nTS > 0:
            ts = reacs[0] + '_' + reac[1]
            msg += 'Task: Extracting MESS data for TS...'
            msg  =  log_msg(msg)

            try:    
                os.chdir('me_files')
                ge, hr = self.read_gehr(ts, 'ts')
                fr = self.extract_mess('ts_fr.me')                 #Copy EStokTP projfrequencies
                fr = fr.split('End')[0] + 'End  '
                os.chdir('..')                        
                msg += 'completed'
                msg  = log_msg(msg)
                pf = tf + ge  + zp + fr + hr
                io.write_file(pf+'\n', ts.strip() + '.pf')
                msg += 'Task: Building MESS input file...'
                msg += 'completed'
                msg  = log_msg(msg)
                species.append(ts)
                speclist.append('ts')
            except IOError:
                msg += 'me_files are missing, check me_file/*, estoktp.log, and output/estoktp.out' 
                msg  = log_msg(msg)
                os.chdir('..')                         
        return species, speclist

    def run(self, args):
        """
        Runs heatform, partition_function, thermp, pac99, and write chemkin file
        """
        import shutil
        import heatform as hf
        import re
            
        reacs    = args.reacs
        prods    = args.prods    
        anharm   = args.anharm
        anovrwrt = args.anovrwrt
        node     = args.node
        meths    = args.meths
        hfbasis  = args.hfbasis
        qtchf    = args.qtchf
        enlevel  = args.enlevel
        hlen     = args.hlen
        self.hfbases = []
        speciess,speclist = self.build_pfinput(args)
        self.dH0   = []
        self.dH298 = []
        anharmbool = False
        deltaH = 0
        if anharm.lower() != 'false':
            anharmbool = True
        for i,species in enumerate(speciess):
            if qtchf[0].lower() not in ['false', 'auto']:
                if len(qtchf) >= i:
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
                if io.check_file(species+'.pf.dat'):
                    os.rename(species + '.pf.dat','pf.dat')
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
                c97file = species + '.c97'
                if io.check_file(c97file):
                    c97text  = io.read_file(c97file)
                    las, has, msg = tc.get_coefficients(c97text)
                chemkinfile = stoich + '.ckin'
                print('Writing chemkin file {0}.\n'.format(chemkinfile))
                method = meths[-1][2]
                chemininput = tc.write_chemkin_file(species, method, deltaH, float(deltaH298), stoich, 0, las, has, chemkinfile)

            
            print('completed')
        return 

    def extract_mess(self,filename):
        """
        Extracts necessary EStokTP output for MESS
        """
        lines = ''
        if io.check_file(filename):
            lines = io.read_file(filename)
        if lines == '':
            print('failed -- ' + filename + ' is empty, check estoktp.log')
        return lines

    def read_gehr(self, s, typ, n=-1):
        
        if n>-1: 
            name = '{0}{1:g}'.format(typ, n+1)
        else:
            name = typ
        ge = self.extract_mess('{}_1dge.me'.format(name))    #Copy EStokTP geometry
        if 'Species' in ge:
            ge = " Species " + s.strip() +  ge.split("Species")[1]
        elif 'Fragment' in ge:
            ge = " Species " + s.strip() +  ge.split("Fragment")[1]
        if 'Core' in ge:
            ge, ge1 = ge.split('Core')
        else:
            ge1 = ''
        hr = self.extract_mess('{}_hr.me'.format(name))          #Copy EStokTP hr data
        if not 'Core' in hr:
            ge = ge + 'Core' +ge1
        elif 'ge' != '':
            hr = hr.split('Quantum')
            hr = hr[0].rstrip() + '\n    Quantum'  + hr[1].lstrip()
        ge = ge.rstrip('End\n') 
        hr += '\nEnd'
        return ge, hr

    def get_fr(self,s, natom, typ, anharm, anovrwrt, meths, node, n=-1, store=False):
       
        if n>-1:
            name = '{0}{1:g}'.format(typ, n+1)
        else:
            name = typ
        if anharm.lower() == 'false':
            fr = self.extract_mess('{}_fr.me'.format(name)) #Copy EStokTP projfrequencies
            fr = fr.split('End')[0]
            fr = fr.replace('Zero','End\n   Zero') 
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
                optlevel = anlevel
            else:
                for meth in meths:
                    if str(anharm) in meth[0]:
                        anlevel = meth[1] + '/' +  meth[2]
                        optlevel = meth[1] + '/' +  meth[2]
                        break
                    else:
                        anlevel = ''
            os.chdir('..')                          
            anfr,fr1, anx,fr2,fr3 = get_anharm(typ, str(n+1), natom, node, anlevel, anovrwrt, s, optlevel.split('/'))  #(xmat with projected out scanned torsional modes)
            fr =  fr1 +  fr2 + fr3
            self.anfreqs.append(anfr)
            self.anxmat.append(anx)
            os.chdir('me_files')
        if store:
            zpve = io.read_file(name + '_zpe.me')
            io.db_store_sp_prop(zpve, mol, 'zpve', prog = prog, optprog = optprog, method= method, optmethod=optmethod, basis=basis, optbasis=optbasis)  
        fr = fr.rstrip('End') + '\n'
        return fr

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

def check_hrs(n, typ, msg=''):
    """
    Checks MC geoms to make sure they are the same inchii as the starting species
    """
    import sys
    import iotools as io

    msg += 'Task: Checking me_files/{1}{0}_hr.me...'.format(str(n),typ)
  
    filename = 'data/' + typ + str(n) + '.dat'
    nrotors = 0
    md = False
    if io.check_file(filename):
        data = io.read_file(filename)
        tmp = data.split('nhind')
        if len(tmp) > 2:
            nrotors = tmp[1].split('\n')[1]
        if len(tmp) > 3:
            md  = True

    data = ''
    filename =  'me_files/' + typ + str(n) +  '_hr.me'
    if io.check_file(filename):
        data = io.read_file(filename)
    else:
        msg += '\nERROR DETECTED: no hr me_file found'
        msg = log_msg(msg)
        return

    if md:
        if 'MultiRotor' in data:
            msg += '\n  MDTau successfully completed'
        else:
            msg += '\nERROR DETECTED: MD scan incomplete'
        filename =  'me_files/' + typ + str(n) +  '_1dhr.me'
        if io.check_file(filename):
            data = io.read_file(filename)
        else:
            msg += '\nERROR DETECTED: no 1dhr me_file found'
            msg = log_msg(msg)
            return

    data = data.split('Rotor') 
    ncomplete = len(data) - 1
    msg += '\n  {0} out of {1} rotors successfully scanned'.format(str(ncomplete), nrotors)
    if int(nrotors) == ncomplete:
        msg += '\n  1DTau has completed successfully'
    else:
        msg += '\nERROR DETECTED: scan incomplete'
    msg  = log_msg(msg)
    return

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
    opts['anlevel'  ] =  anlevel
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

def get_paths(dic, key):
    """
    Finds a value in a dic
    """
    val = None
    if key in dic:
        val = dic[key]
    else:
        print "path for {} not given in configfile".format(key)
    return val
        
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
      

def log_msg(msg):
    print(msg)
    return ''

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
     

    #####  Get arguments  ##########
    ################################
    import sys
    import os
    import config

    optionfile = 'input.dat'
    torspath   = os.path.dirname(os.path.realpath(sys.argv[0]))
    configfile = torspath + os.path.sep + 'configfile.txt'
    if len(sys.argv) > 1:
        optionfile = sys.argv[1]
        if len(sys.argv) > 2:
            configfile = sys.argv[2]
     
    args   = config.ARGS(  optionfile)
    Config = config.CONFIG(configfile)
    paths  = Config.path_dic()
    paths['torsscan'] = torspath

    sys.path.insert(0, get_paths(paths,     'bin'))
    sys.path.insert(0, get_paths(paths,     'qtc'))
    sys.path.insert(0, get_paths(paths,'torsscan'))

    #####  Build and Run EStokTP  ######
    ####################################
    import obtools as ob
    import iotools as io
    import tctools as tc

    es   = ES(args, paths)             
    if args.restart < 5:
        es.build_subdirs()
        es.build_files()
        es.execute()
     
    #check for failures
    if ("Opt" in args.jobs and not "Opt_1" in args.jobs):
        es.check_geoms(es.nsamps)
    if ("1dTau" in args.jobs or 'MdTau' in args.jobs):
        for i in range(len(args.reacs)):
            check_hrs(i+1,'reac')
        for i in range(len(args.prods)):
            check_hrs(i+1,'prod')
        es.me_file_abs_path()

    import results 
    rs = results.RESULTS(args)
    args.hlen = rs.get_hlen()
    args.optlevel = rs.optlevel
    args.enlevel = rs.enlevel
    args.taulevel = rs.taulevel
    #######  Build and run thermo  #########
    ########################################
    rs.thermo = False
    if args.alltherm.lower() == 'true':
        rs.thermo = True
        import patools as pa
        args.symnums = es.symnums
        thermo = THERMO(args)
        thermo.run(args)

    #######  Parse results  #########
    ########################################
    if args.parseall.lower() == 'true':
        if rs.thermo:
             rs.get_results(thermo)
        else:
             rs.get_results()

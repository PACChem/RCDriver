#!/usr/bin/python

import os
import build 

def build_files(args, paths, nodes = 1, msg=''):
    """
    Runs the build functions for reacn.dat, prodn.dat, theory.dat, and estoktp.dat
    Requirements: QTC, Openbabel, Pybel (OR prepared cartesian coordinate files) and x2z
    """

    import sys
    sys.path.insert(0, paths['qtc'])
    global io, ob
    import patools as pa
    import iotools as io
    import obtools as ob

    build_subdirs()
    os.chdir('./data')
    stoichs= []
    symnums= []

    if not 'MdTau' in args.jobs:
        args.mdtype = ''
    #Create Read, Prod, and TS objects from parameters
    params = (args.nsamps, args.abcd,nodes,args.interval,args.nsteps,args.XYZ,args.xyzstart,args.mdtype)
    Reac   = build.MOL(paths, params, 'reac')
    Prod   = build.MOL(paths, params, 'prod')
    params = (args.nsamps,args.abcd,nodes,args.interval,args.nsteps,args.XYZ,'start',args.mdtype)
    TS     = build.MOL(paths, params,   'ts') 

    reacs = args.reacs
    prods = args.prods
    
    key = set_keys(args.reactype)
    i,j,k = 0,0,0
    TSprops = [0, 0, [], []] #charge, spin, angles, atoms

    #Build reacn.dat
    for i, reac in enumerate(reacs,start=1):

        msg += 'Task: Building reac{:g}.dat...'.format(i)
        msg  = log_msg(msg)
        args.reacs, angles, atoms, args.jobs, stoichs, symnums = build_mol_dat(Reac, reacs, i, stoichs, symnums, args.jobs, args.mdtype)
        TSprops = prep_reacs4TS(Reac, reac, i, key, angles, atoms, TSprops, args.nTS, paths)
        nsamps = Reac.nsamps
        msg += 'completed'
        msg  = log_msg(msg)

    #Build prodn.dat
    for j, prod in enumerate(prods,start=1):
        msg += 'Task: Building prod{:g}.dat...'.format(j)
        msg  = log_msg(msg)
        args.prods, angles, atoms, args.jobs, stoichs, symnums = build_mol_dat(Prod, prods, j, stoichs, symnums, args.jobs, args.mdtype)
        msg += 'completed'
        msg  = log_msg(msg)
    
    #Build TS, wellr, and wellp.dat
    tstype = ['ts','wellr','wellp']
    TS.ijk  = Reac.ijk
    TS.sort = Reac.sort
    for k in range(args.nTS):
        msg += 'Task: Building ' + tstype[k] +  '.dat...'
        msg  = log_msg(msg)
        TS.charge = TSprops[0]
        TS.mult   = int(2.*TSprops[1] + 1)
        TS.symnum = ' 1'
        if k == 0:
            zmatstring = TS.build(tstype[k], prod, TSprops[2], TSprops[3])
        else:
            params = ('1', args.abcd,nodes,args.interval,args.nsteps,'False','start','MdTau' in args.jobs)
            TS   = build.MOL(paths, params,'ts') 
            TS.charge = TSprops[0]
            TS.mult   = int(2.*TSprops[1] + 1)
            TS.symnum = ' 1'
            zmatstring =TS.build(tstype[k], [], [])
        zmat = tstype[k] + '.dat'
        io.write_file(zmatstring, zmat)
        msg += 'completed'
        msg  = log_msg(msg)

    stoich = stoichs[0]
    mol    =   reac[0]
    #Builds theory.dat
    msg += 'Task: Building theory.dat...'
    msg  = log_msg(msg)
    theostring = build.build_theory(args.meths,args.nTS,args.optoptions)
    io.write_file(theostring, 'theory.dat')
    msg += 'completed'
    msg  = log_msg(msg)
    
    #Builds estoktp.dat to restart at any step
    msg += 'Task: Building estoktp.dat...'
    msg  = log_msg(msg)
    jobs = update_jobs(args.jobs, args.restart)
    params    = (stoich, args.reactype, args.coresh,args.coresl,args.mem)
    eststring = build.build_estoktp(params,jobs,i,j,args.nTS)
    io.write_file(eststring, 'estoktp.dat')
    msg += 'completed'
    msg  = log_msg(msg)

    os.chdir('..')

    return stoichs, symnums
 

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
        
def build_subdirs(msg=''):
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
 
def set_keys(reactype):
    """
    These keys will later replace generated reac file if there is a TS search
    """
    if 'abstraction' in reactype.lower():
        key = ['[CH3]','[O]','[O][O]','O[O]','[OH]','[H]','O=O']
    elif  'addition' in reactype.lower():
        key = ['[O][O]']
    else: 
        key = []
    return key

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

def prep_reacs4TS(MOL, reac, i, key, angles, atoms, tsprops, nTS, paths):
    """
    Sets transition state charge, mult, angles, and atoms respectively based 
    on the reactants.  May use abstractor templates if key matches
    """
    if nTS > 0:
        if i == 1:
            sort = MOL.sort
        tsprops[0] += MOL.charge
        tsprops[1] += 1./2 * (float(MOL.mult) - 1)
        if i == 1:
            tsprops[2], tsprops[3]  = angles, atoms
        elif reac in key:
            import shutil
            shutil.copyfile(paths['torsscan'] + '/abstractors/' + reac + '.dat','reac2.dat')
    return tsprops

    
def build_mol_dat(MOL, mollist, n, stoichs, symnums, jobs, mdtype='auto'):
    """
    Builds reac1.dat, reac2.dat, prod1.dat etc
    """

    mol = mollist[n-1]
    atoms, measure, angles  = MOL.cart2zmat(mol)
    mollist[n-1] = mol.split('_')[0]

    if mdtype.lower() == 'auto':
        MOL.MDTAU, jobs = prepare_mdtau(len(angles), jobs)

    zmatstring = MOL.build(n, mol, angles, atoms,  measure)
    zmat = MOL.typemol + str(n) + '.dat'
    io.write_file(zmatstring, zmat)
    stoichs.append(MOL.stoich)
    symnums.append(MOL.symnum)
    return mollist, angles, atoms, jobs, stoichs, symnums


def execute(paths, node, back = '', msg = ''):
    """
    Runs EStokTP on a given blues node (default debug)
    use 0 in input file to run on login
    use d or debug to just make input files and not run EStokTP
    Requirements: PACC member on blues
    """
    g09     = get_paths(paths,  'g09')
    gcc     = get_paths(paths,  'gcc')
    intel   = get_paths(paths,  'intel')
    estoktp = get_paths(paths,'estoktp')
    msg += 'Task: Submitting EStokTP job...'
    msg  = log_msg(msg)
    if node == '0':
       # os.system('soft add +g09; soft add +gcc-5.3; /home/elliott/Packages/EStokTP/exe/estoktp.x >& estoktp.log')
        os.system('{0}; {1}; {2}; {3}  >& estoktp.log'.format(gcc, intel, g09, estoktp))
        msg += 'completed'
    elif node == 'd' or node == 'debug':
        msg += 'task skipped'
    else:
        ssh = get_paths(paths, 'ssh')
        os.system('exec {3} -n {4} "cd `pwd`;{0}; {1}; {2}; {5} >& estoktp.log {6}"'.format(gcc, intel, g09, ssh, node, estoktp, back))
        msg += 'completed'
    msg  = log_msg(msg)
    return
    
def check_geoms(qtc, name, nsamps,msg=''):
    """
    Checks MC geoms to make sure they are the same inchii as the starting species
    """
    import sys
    sys.path.insert(0,qtc)
    import iotools as io
    import obtools as ob

    msg += 'Task: Checking level0 geometries...'
    msg  = log_msg(msg)
    n = 2
    filename =  'geoms/'  + name + '_' + '1'.zfill(n) + '.xyz'
    lowfilename   = filename
    coords = io.read_file(filename)
    lowcoords = coords
    mol = ob.get_mol(coords)
    name =  ob.get_inchi_key(mol)
    energy = float(coords.split('\n')[1])
    for i in range(2, int(nsamps) + 1):
        filename =  'geoms/' + name + '_' + '{}'.format(i).zfill(n) + '.xyz'
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
    return lowfilename

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

def me_file_abs_path():
    """
    Replaces relative path in mdhr file with absolute path
    """
    import iotools as io
    if io.check_file('me_files/reac1_hr.me'):
        lines = io.read_file('me_files/reac1_hr.me')
        if "PotentialEnergySurface[kcal/mol]" in lines:
            before, after = lines.split("PotentialEnergySurface[kcal/mol]")
            after = after.split('\n')
            after[0] = after[0].replace('./',io.pwd() + '/')
            lines = before + "PotentialEnergySurface[kcal/mol]" + '\n'.join(after)
            io.write_file(lines,'me_files/reac1_hr.me')
    return

def log_msg(msg):
    print(msg)
    return ''


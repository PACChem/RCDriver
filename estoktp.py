#!/usr/bin/python

import os
import build 
import logging
log = logging.getLogger(__name__)
import shutil 

def build_files(args, paths, nodes = 1):
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
    mdtype =  args.mdtype 
    if not 'MdTau' in args.jobs:
       mdtype = ''
    #Create Read, Prod, and TS objects from parameters
    params = (args.nsamps, args.abcd,nodes,args.interval,args.nsteps,args.XYZ,args.xyzstart,mdtype)
    Reac   = build.MOL(paths, params, 'reac')
    Prod   = build.MOL(paths, params, 'prod')
    params = (args.nsamps,args.abcd,nodes,args.interval,args.nsteps,args.XYZ,'start',mdtype)
    TS     = build.MOL(paths, params,   'ts') 

    reacs = args.reacs
    prods = args.prods
    
    key = set_keys(args.reactype)
    i,j,k = 0,0,0
    TSprops = [0, 0, [], [],[]] #charge, spin, angles, atoms, sort
    foundlist = []
    #Build reacn.dat
    for i, reac in enumerate(reacs,start=1):

        msg = 'Building reac{:g}.dat'.format(i)
        log.debug(msg)
        args.reacs, angles, atoms, args.jobs, foundlist, stoichs, symnums = build_mol_dat(Reac, reacs, i, stoichs, symnums, args.jobs, foundlist, args.select[i-1], mdtype)
        TSprops = prep_reacs4TS(Reac, reac, i, key, angles, atoms, Reac.sort,TSprops, args.nTS, paths)
        nsamps = Reac.nsamps
        msg = 'Completed'
        log.info(msg)

    #Build prodn.dat
    for j, prod in enumerate(prods,start=1):
        msg = 'Building prod{:g}.dat'.format(j)
        log.debug(msg)
        args.prods, angles, atoms, args.jobs, foundlist, stoichs, symnums = build_mol_dat(Prod, prods, j, stoichs, symnums, args.jobs, foundlist, args.select[i+j-1], mdtype)
        msg = 'Completed'
        log.info(msg)
    
    #Build TS, wellr, and wellp.dat
    tstype = ['ts','wellr','wellp']
    TS.ijk  = Reac.ijk
    TS.sort = TSprops[4]
    for k in range(args.nTS):
        msg = 'Building ' + tstype[k] +  '.dat'
        log.debug(msg)
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
        msg = 'Completed'
        log.info(msg)

    mol    =   reac[0]
    #Builds me_head.dat
    if args.nTS > 0:
        msg = 'Building me_head.dat'
        log.debug(msg)
        headstring = build.build_mehead()
        io.write_file(headstring, 'me_head.dat')
        msg = 'Completed'
        log.info(msg)
    #Builds theory.dat
    msg = 'Building theory.dat'
    log.debug(msg)
    theostring = build.build_theory(args.meths,args.nTS,args.zedoptions,args.oneoptions)
    io.write_file(theostring, 'theory.dat')
    msg = 'Completed'
    log.info(msg)
    
    #Builds estoktp.dat to restart at any step
    msg = 'Building estoktp.dat'
    log.debug(msg)
    jobs = update_jobs(args.jobs, args.restart)
    params    = (stoichs, args.reactype, args.coresh,args.coresl,args.mem,args.esoptions)
    eststring = build.build_estoktp(params,jobs,i,j,args.nTS,args.xyzstart,foundlist)
    io.write_file(eststring, 'estoktp.dat')
    msg = 'Completed'
    log.info(msg)

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
        
def build_subdirs():
    """
    Builds data and output subdirectories
    """
    msg = 'Building directories'
    log.debug(msg)
    if not os.path.exists('./data'):
        os.makedirs('./data') 
    if not os.path.exists('./output'):
        os.makedirs('./output') 
    msg = 'Completed'
    log.info(msg)
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

def prep_reacs4TS(MOL, reac, i, key, angles, atoms, sort, tsprops, nTS, paths):
    """
    Sets transition state charge, mult, angles, and atoms respectively based 
    on the reactants.  May use abstractor templates if key matches
    """
    if nTS > 0:
        if i == 1:
            sort = MOL.sort
        tsprops[0] += MOL.charge
        tsprops[1] = max(abs(tsprops[1] + 1./2 * (float(MOL.mult) - 1)), abs(tsprops[1] - 1./2 * (float(MOL.mult)-1)))
        tsprops[2].append(angles)
        tsprops[3].append(atoms)
        if i == 1:
            tsprops[4] = sort
        elif reac in key:
            shutil.copyfile(paths['torsscan'] + '/abstractors/' + reac + '.dat','reac2.dat')
    return tsprops

    
def build_mol_dat(MOL, mollist, n, stoichs, symnums, jobs, foundlist, select, mdtype='auto'):
    """
    Builds reac1.dat, reac2.dat, prod1.dat etc
    """

    mol = mollist[n-1]
    atoms, measure, angles, found, msg  = MOL.cart2zmat(mol, select)
    #mollist[n-1] = mol.split('_')[0]
    if mdtype.lower() == 'auto':
        MOL.MDTAU, jobs = prepare_mdtau(len(angles), jobs)
    if mdtype:
        log.info(msg)
    zmatstring = MOL.build(n, mol, angles, atoms,  measure)
    zmat = MOL.typemol + str(n) + '.dat'
    io.write_file(zmatstring, zmat)
    stoichs.append(MOL.stoich)
    symnums.append(MOL.symnum)
    foundlist.append(found)
    return mollist, angles, atoms, jobs, foundlist, stoichs, symnums


def execute(paths, node, back = ''):
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
    msg = 'Submitting EStokTP job to node {}'.format(node)
    log.debug(msg)
    if node == '0':
       # os.system('soft add +g09; soft add +gcc-5.3; /home/elliott/Packages/EStokTP/exe/estoktp.x >& estoktp.log')
        os.system('{0}; {1}; {2}; {3}  >& estoktp.log'.format(gcc, intel, g09, estoktp))
        msg = 'Completed'
    elif node == 'd' or node == 'debug':
        msg = 'Task skipped'
    else:
        ssh = get_paths(paths, 'ssh')
        os.system('exec {3} -n {4} "cd `pwd`;{0}; {1}; {2}; {5} >& estoktp.log {6}"'.format(gcc, intel, g09, ssh, node, estoktp, back))
        msg = 'Completed'
    log.info(msg)
    return
    
def check_geoms(qtc, name, nsamps):
    """
    Checks MC geoms to make sure they are the same inchii as the starting species
    """
    import sys
    sys.path.insert(0,qtc)
    import iotools as io
    import obtools as ob

    msg = 'Checking level0 geometries'
    log.debug(msg)
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
    msg = '\nMonte Carlo sampling successfully found geom.xyz!\n'
    log.info(msg)
    return lowfilename

def check_hrs(n, typ):
    """
    Checks MC geoms to make sure they are the same inchii as the starting species
    """
    import sys
    import iotools as io

    msg = 'Checking me_files/{1}{0}_hr.me'.format(str(n),typ)
  
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
        msg = '\nNo hr me_file found'
        log.error(msg)
        return

    if md:
        if 'MultiRotor' in data:
            msg = '\nMDTau successfully completed'
            log.info(msg)
        else:
            msg = '\nMD scan incomplete'
            log.error(msg)
        filename =  'me_files/' + typ + str(n) +  '_1dhr.me'
        if io.check_file(filename):
            data = io.read_file(filename)
        else:
            msg += '\nNo 1dhr me_file found'
            log.error(msg)
            return

    data = data.split('Rotor') 
    ncomplete = len(data) - 1
    msg = '\n{0} out of {1} rotors successfully scanned'.format(str(ncomplete), nrotors)
    if int(nrotors) == ncomplete:
        msg = '\n1DTau has completed successfully'
    else:
        msg = '\nScan incomplete'
        log.error(msg)
    log.info(msg)
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

def gather_mcgeoms(nodes):
    j,k = [0,0,0,0,0,0,0],[0,0,0,0,0,0,0]
    for i, node in enumerate(nodes):
        for geom in os.listdir(node + '/geoms'):
            if 'reac1' in geom:
                j[0] += 1
                index = j[0]
            elif 'reac2' in geom:
                j[1] += 1
                index = j[1]
            elif 'prod1' in geom:
                j[2] += 1
                index = j[2]
            elif 'prod2' in geom:
                j[3] += 1
                index = j[3]
            elif 'ts' in geom:
                j[4] += 1
                index = j[4]
            elif 'wellr' in geom:
                j[5] += 1
                index = j[5]
            elif 'wellp' in geom:
                j[6] += 1
                index = j[6]
            logging.info('Copying ' + node + '/geoms/' + geom +  ' to  ' + 'geoms/{}_{}.xyz'.format(geom.split('_')[0], str(index).zfill(2)))
            shutil.copy(node + '/geoms/' + geom, 'geoms/{}_{}.xyz'.format(geom.split('_')[0], str(index).zfill(2)))
        for geom in os.listdir(node + '/output'):
            if 'opt_' in geom:
                if 'reac1' in geom:
                    k[0] += 1
                    index = k[0]
                elif 'reac2' in geom:
                    k[1] += 1
                    index = k[1]
                elif 'prod1' in geom:
                    k[2] += 1
                    index = k[2]
                elif 'prod2' in geom:
                    k[3] += 1
                    index = k[3]
                elif 'ts' in geom:
                    k[4] += 1
                    index = k[4]
                elif 'wellr' in geom:
                    k[5] += 1
                    index = k[5]
                elif 'wellp' in geom:
                    k[6] += 1
                    index = k[6]
                logging.info('Copying ' + node + '/output/' + geom +  'to output/{}opt_{}.out'.format(geom.split('opt')[0], str(index).zfill(2)))
                shutil.copy(node + '/output/' + geom, 'output/{}opt_{}.out'.format(geom.split('opt')[0], str(index).zfill(2)))
    return j, k 

def run_level0(args, paths):
    stoichs, symnums = build_files(args, paths, len(args.nodes))
    if not 'd' in args.nodes[0]:
        if len(args.nodes) > 1:
            for i, node in enumerate(args.nodes):
                io.rmrf(node)
                io.mkdir(node)
                io.cd(node)
                shutil.copytree('../data', 'data')
                execute(paths, node, '&')
                io.cd('..')
            running = True
            import time
            outlength = {}
            for node in args.nodes:
                outlength[node] = 0
            while (running):
                running = False
                for node in args.nodes:
                    filename = node + '/geom.log'
                    if not io.check_file(filename):
                        running = True
                        #log.debug( 'waiting on node {}'.format(node))
                    else: 
                        lines = io.read_file(filename)
                        newlength = len(lines)
                        if 'termination' in lines and newlength == outlength[node]:
                            pass
                        else:
                            running = True
                            outlength[node] = newlength
                            #log.debug( 'waiting on node {}'.format(node))
                time.sleep(30)
            log.info('Jobs Completed')
            io.mkdir('geoms')
            j, k = gather_mcgeoms(args.nodes)
            for i in range(len(args.reacs)):
                filename = check_geoms(paths['qtc'], 'reac' + str(i+1), j[i])
                filename = filename.split('/')[1].split('_')[0] + '_opt_' +  filename.split('_')[1]
                filename = 'output/' + filename.replace('.xyz','.out')
                shutil.copy(filename, 'output/reac' + str(i+1) + '_opt.out')
            for i in range(len(args.prods)):
                filename = check_geoms(paths['qtc'], 'prod' + str(i+1), j[i+len(args.reacs)])
                filename = filename.split('/')[1].split('_')[0] + '_opt_' +  filename.split('_')[1]
                filename = 'output/' + filename.replace('.xyz','.out')
                shutil.copy(filename, 'output/prod' + str(i+1) + '_opt.out')
        else:
             execute(paths, args.nodes[0])
    return

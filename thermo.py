#!/home/keceli/anaconda2/bin/python
import os

def extract_mess(filename):
    """
    Extracts necessary EStokTP output for MESS
    """
    lines = ''
    if io.check_file(filename):
        lines = io.read_file(filename)
    if lines == '':
        print('failed -- ' + filename + ' is empty, check estoktp.log')
    return lines


def get_fr(s, natom, typ, anharm, anovrwrt, anfreqs, anxmat, meths, node, n=-1, store=False):
   
    if n>-1:
        name = '{0}{1:g}'.format(typ, n+1)
    else:
        name = typ
    if anharm.lower() == 'false':
        fr = extract_mess('{}_fr.me'.format(name)) #Copy EStokTP projfrequencies
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
        anfreqs.append(anfr)
        anxmat.append(anx)
        os.chdir('me_files')
    if store:
        zpve = io.read_file(name + '_zpe.me')
        io.db_store_sp_prop(zpve, mol, 'zpve', prog = prog, optprog = optprog, method= method, optmethod=optmethod, basis=basis, optbasis=optbasis)  
    fr = fr.rstrip('End') + '\n'
    return fr, anfreqs, anxmat

def read_gehr(s, typ, n=-1):
    
    if n>-1: 
        name = '{0}{1:g}'.format(typ, n+1)
    else:
        name = typ
    ge = extract_mess('{}_1dge.me'.format(name))    #Copy EStokTP geometry
    if 'Species' in ge:
        ge = " Species " + s.strip() +  ge.split("Species")[1]
    elif 'Fragment' in ge:
        ge = " Species " + s.strip() +  ge.split("Fragment")[1]
    if 'Core' in ge:
        ge, ge1 = ge.split('Core')
    else:
        ge1 = ''
    hr = extract_mess('{}_hr.me'.format(name))          #Copy EStokTP hr data
    if not 'Core' in hr:
        ge = ge + 'Core' +ge1
    elif 'ge' != '':
        hr = hr.split('Quantum')
        hr = hr[0].rstrip() + '\n    Quantum'  + hr[1].lstrip()
    ge = ge.rstrip('End\n') 
    hr += '\nEnd'
    return ge, hr

def build_pfinput(args, msg = ''):
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

    anfreqs  = []
    anxmat   = []
    hfbasis  = args.hfbasis

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
            ge, hr = read_gehr(reac, 'reac',  n)
            if len(symnums) > n:
                symnum = symnums[n]
            else:
                symnum = 1
            ge = ge.replace('SymmetryFactor    1.0000000000000','SymmetryFactor     {:.2f}'.format(float(symnum)))
            fr, anfreqs, anxmat = get_fr(reac, natom, 'reac', anharm, anovrwrt, anfreqs, anxmat, meths, node, n, args.store)
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
            ge, hr = read_gehr(prod, 'prod', n)
            if len(symnums) > n+len(reacs):
                symnum = symnums[n+len(reacs)]
            else:
                symnum = 1
            ge = ge.replace('SymmetryFactor    1.0000000000000','SymmetryFactor        {:.2f}'.format(float(symnum)))
            fr, anfreqs, anxmat = get_fr(prod, natom, 'prod', anharm, anovrwrt, anfreqs, anxmat, meths, node, n, args.store)
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
            ge, hr = read_gehr(ts, 'ts')
            fr = extract_mess('ts_fr.me')                 #Copy EStokTP projfrequencies
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
    return species, speclist, anfreqs, anxmat

def run(args, paths):
    """
    Runs heatform, partition_function, thermp, pac99, and write chemkin file
    """
    import sys
    sys.path.insert(0, paths['qtc'])
    global pa, io, ob
    import patools as pa
    import iotools as io
    import obtools as ob
    import heatform as hf
    import shutil
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
    hfbases = []
    speciess,speclist = build_pfinput(args)
    dH0   = []
    dH298 = []
    anharmbool = False
    deltaH = 0
    if anharm.lower() != 'false':
        anharmbool = True
    for i,species in enumerate(speciess):
        if qtchf[0].lower() not in ['false', 'auto']:
            if len(qtchf) >= i:
                deltaH = float(qtchf[i])
                hfbasis = ['N/A']
                hfbases.append(hfbasis)
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
                deltaH, hfbasis = hf.main(species,logfile,E=energy,basis=hfbasis,anharm=anharmbool,enlevel=enlevel)
                hfbases.append(hfbasis)
            else:
                deltaH = 0.00
        dH0.append(deltaH)
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
            dH298.append(deltaH298)
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
    return dH0, dH298, hfbases 


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


def log_msg(msg):
    print(msg)
    return ''

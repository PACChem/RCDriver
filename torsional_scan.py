#!/home/keceli/anaconda2/bin/python

import os
import sys
import numpy as np
import rmg_reader as rg

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
    import estoktp as es

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

    sys.path.insert(0, es.get_paths(paths,     'bin'))
    sys.path.insert(0, es.get_paths(paths,     'qtc'))
    sys.path.insert(0, es.get_paths(paths,'torsscan'))

    #####  Build and Run EStokTP  ######
    ####################################
    import obtools as ob
    import iotools as io
    import tctools as tc
    import shutil
 
    symnums = []
    samps = None
    if args.restart < 5:
        j,k = [0,0,0,0,0,0,0],[0,0,0,0,0,0,0]
        index = 0
        if len(args.nodes) > 0 and "Opt" in args.jobs and args.restart < 1:
            alljobs = args.jobs
            args.jobs = ["Opt"]
            stoichs, symnums = es.build_files(args, paths, len(args.nodes))
            if not 'd' in args.nodes[0]:
                for i, node in enumerate(args.nodes):
                    io.rmrf(node)
                    io.mkdir(node)
                    io.cd(node)
                    shutil.copytree('../data', 'data')
                    es.execute(paths, node, '&')
                    io.cd('..')
                running = True
                import time
                while (running):
                    running = False
                    for node in args.nodes:
                        filename = node + '/output/estoktp.out'
                        if not io.check_file(filename):
                            running = True
                            print 'waiting on node {}'.format(node)
                        else: 
                            if len(io.read_file(filename)) < 20:
                                running = True
                                print 'waiting on node {}'.format(node)
                    time.sleep(60)
                io.mkdir('geoms')
                for i, node in enumerate(args.nodes):
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
                            shutil.copy(node + '/output/' + geom, 'output/{}opt_{}.out'.format(geom.split('opt')[0], str(index).zfill(2)))
                args.restart = 1
                args.jobs = alljobs
                for i in range(len(args.reacs)):
                    filename = es.check_geoms(paths['qtc'], 'reac' + str(i+1), j[i])
                    filename = filename.split('/')[1].split('_')[0] + '_opt_' +  filename.split('_')[1]
                    filename = 'output/' + filename.replace('.xyz','.out')
                    shutil.copy(filename, 'output/reac' + str(i+1) + '_opt.out')
                for i in range(len(args.prods)):
                    filename = es.check_geoms(paths['qtc'], 'prod' + str(i+1), j[i+len(args.reacs)])
                    filename = filename.split('/')[1].split('_')[0] + '_opt_' +  filename.split('_')[1]
                    filename = 'output/' + filename.replace('.xyz','.out')
                    shutil.copy(filename, 'output/prod' + str(i+1) + '_opt.out')
        if "1dTau" in args.jobs:
            alljobs = args.jobs
            negvals = True
            while (negvals):
                args.jobs = ["Opt_1", "1dTau"]
                negvals = False
                stoichs, symnums = es.build_files(args, paths)
                es.execute(paths, args.nodes[0])
                args.restart = 3
                negvals = False
                for i in range(len(args.reacs)):
                    lowene = 0.
                    lowenefile = None
                    if io.check_file('me_files/reac' +  str(i+1) + '_hr.me'):
                        hr = io.read_file('me_files/reac' + str(i+1) + '_hr.me')
                        hr = hr.split('Rotor')
                        startkey = 'Potential' 
                        for j, rotor in enumerate(hr[1:]):
                            pot = rotor.split(startkey)[1]
                            pot = pot.splitlines()[1]
                            for k, ene in enumerate(pot.split()):
                                if float(ene) < lowene:
                                    lowene = float(ene)
                                    lowenefile = 'hr_geoms/geom_isp' + str(i+1).zfill(2) + '_hr' + str(j+1).zfill(2) + '_hpt' + str(k+1).zfill(2) + '.xyz'
                    if lowenefile:
                        negvals = True
                        xyz = io.read_file(lowenefile)
                        slabel = ob.get_slabel(ob.get_mol(xyz))
                        io.write_file(xyz,slabel + '.xyz')
                        args.restart = 1
                        args.XYZ = 'true'
                        args.xyzstart = '0'
                        print 'Lower configuration found in 1dTau. Restarting at Level1. Saved geometry to {}'.format(slabel + '.xyz')
                for l in range(len(args.prods)):
                    lowene = 0.
                    lowenefile = None
                    if io.check_file('me_files/prod' +  str(l+1) + '_hr.me'):
                        hr = io.read_file('me_files/prod' + str(l+1) + '_hr.me')
                        hr = hr.split('Rotor')
                        startkey = 'Potential' 
                        for j, rotor in enumerate(hr[1:]):
                            pot = rotor.split(startkey)[1]
                            pot = pot.splitlines()[1]
                            for k, ene in enumerate(pot.split()):
                                if float(ene) < lowene:
                                    lowene = ene
                                    lowenefile = 'hr_geoms/geom_isp' + str(i+l+1).zfill(2) + '_hr' + str(j+1).zfill(2) + '_hpt' + str(k+1).zfill(2) + '.xyz'
                    if lowenefile:
                        negvals = True
                        xyz = io.read_file(lowenefile)
                        slabel = ob.get_slabel(ob.get_mol(xyz))
                        io.write_file(xyz, slabel + '.xyz')
                        args.restart = 1
                        args.XYZ = 'true'
                        args.xyzstart = '0'
                        print 'Lower configuration found in 1dTau. Restarting at Level1. Saved geometry to {}'.format(slabel + '.xyz')
            args.jobs = alljobs
        stoichs, symnums = es.build_files(args, paths)
        es.execute(paths, args.nodes[0])
    if ("1dTau" in args.jobs or 'MdTau' in args.jobs):
        for i in range(len(args.reacs)):
            es.check_hrs(i+1,'reac')
        for i in range(len(args.prods)):
            es.check_hrs(i+1,'prod')
        es.me_file_abs_path()

    #######  Parse results  #########
    ########################################
    import results 
    rs = results.RESULTS(args, paths)
    args.hlen = rs.get_hlen()
    args.optlevel = rs.optlevel
    args.enlevel = rs.enlevel
    args.taulevel = rs.taulevel

    if args.parseall.lower() == 'true':
         rs.get_results()
    #######  Build and run thermo  #########
    ########################################
    import thermo
    rs.thermo = False
    if args.alltherm.lower() == 'true':
        rs.thermo = True
        args.symnums = symnums
        rs.dH0, rs.dH298, rs.hfbases, rs.anfreqs, rs.anxmat = thermo.run(args, paths)
        if args.parseall.lower() == 'true':
             rs.get_thermo_results()


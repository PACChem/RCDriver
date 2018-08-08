#!/home/keceli/anaconda2/bin/python

import os
import sys
import numpy as np
import rmg_reader as rg
import logging
import config  
import estoktp as es
log   = logging.getLogger(__name__)
import argparse

def random_cute_animal():
    import random 
    msg = random.choice(["""\n\t\t   TORSSCAN
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

     """\n    
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
    """\n
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
    """])
    return msg

if __name__ == "__main__":
  
   
    torspath   = os.path.dirname(os.path.realpath(sys.argv[0]))
    configfile = torspath + os.path.sep + 'configfile.txt'
    #####  Get arguments  ##########
    ################################
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                            description="TorsScan")
    parser.add_argument('-i','--inputfile',   type=str,   default='input.dat')
    parser.add_argument('-c','--configfile',   type=str,   default=configfile)
    parser.add_argument('-o','--outputfile',   type=str,   default='')
    inputs = parser.parse_args()

    inputfile  = inputs.inputfile
    configfile = inputs.configfile 
    outputfile = inputs.outputfile
   
    args   = config.ARGS(  inputfile)
    Config = config.CONFIG(configfile,outputfile)
    paths  = Config.path_dic()
    paths['torsscan'] = torspath

    sys.path.insert(0, es.get_paths(paths,     'bin'))
    sys.path.insert(0, es.get_paths(paths,     'qtc'))
    sys.path.insert(0, es.get_paths(paths,'torsscan'))

    log.info(random_cute_animal())
    #####  Build and Run EStokTP  ######
    ####################################
    import obtools as ob
    import iotools as io
    import tctools as tc
    import shutil
 
    symnums = []
    samps = None
    if args.restart < 5:
        index = 0
        if "Opt" in args.jobs and args.restart < 1:
            alljobs = args.jobs
            args.jobs = ["Opt"]
            es.run_level0(args, paths)
            args.restart = 1
            args.jobs = alljobs
        if "1dTau" in args.jobs and args.restart < 3:
            alljobs = args.jobs
            negvals = True
            attempt = 0
            while (negvals and attempt < 3):
                attempt += 1
                args.jobs = ["Opt_1", "1dTau"]
                negvals = False
                stoichs, symnums = es.build_files(args, paths)
                es.execute(paths, args.nodes[0])
                shutil.copy('output/estoktp.out','output/estoktp_l1.out')
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
                                if float(ene) - lowene < .2:
                                    lowene = float(ene)
                                    lowenefile = 'hr_geoms/geom_isp' + str(i+1).zfill(2) + '_hr' + str(j+1).zfill(2) + '_hpt' + str(k+1).zfill(2) + '.xyz'
                    if lowenefile:
                        xyz = io.read_file(lowenefile)
                        slabel = ob.get_slabel(ob.get_mol(xyz))
                        if slabel == ob.get_slabel(args.reacs[i]):
                            negvals = True
                            slabel = ob.get_smiles_filename(slabel)
                            io.write_file(xyz,slabel + '.xyz')
                            args.restart = 1
                            args.XYZ = 'true'
                            args.xyzstart = '0'
                            log.warning( 'Lower configuration found in 1dTau. Restarting at Level1. Saved geometry to {}'.format(slabel + '.xyz'))
                        else: 
                            log.warning( 'Lower configuration found in 1dTau. But has different smiles: {} vs. {}'.format(slabel, ob.get_slabel(args.reacs[i])))
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
                                if float(ene) -lowene < .2:
                                    lowene = ene
                                    lowenefile = 'hr_geoms/geom_isp' + str(i+l+2).zfill(2) + '_hr' + str(j+1).zfill(2) + '_hpt' + str(k+1).zfill(2) + '.xyz'
                    if lowenefile:
                        xyz = io.read_file(lowenefile)
                        slabel = ob.get_slabel(ob.get_mol(xyz))
                        if slabel == ob.get_slabel(args.prods[l]):
                            negvals = True
                            slabel = ob.get_smiles_filename(slabel)
                            io.write_file(xyz, slabel + '.xyz')
                            args.restart = 1
                            args.XYZ = 'true'
                            args.xyzstart = '0'
                            log.warning( 'Lower configuration found in 1dTau. Restarting at Level1. Saved geometry to {}'.format(slabel + '.xyz'))
                        else: 
                            log.warning( 'Lower configuration found in 1dTau. But has different smiles: {} vs. {}'.format(slabel, ob.get_slabel(args.prods[l])))
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

    if args.parseall.lower() == 'true' or args.alltherm.lower() == 'true':
         rs.get_results()
    #######  Build and run thermo  #########
    ########################################
    import thermo
    rs.thermo = False
    if args.alltherm.lower() == 'true':
        rs.thermo = True
        args.symnums = symnums
        rs.dH0, rs.dH298, rs.hfbases, rs.anfreqs, rs.anxmat = thermo.run(args, paths, rs.d)
        if args.parseall.lower() == 'true':
             rs.get_thermo_results()


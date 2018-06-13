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
        j,k = 0,0
        if len(args.nodes) > 1 and "Opt" in args.jobs and args.restart < 1:
            alljobs = args.jobs
            args.jobs = ["Opt"]
            stoichs, symnums = es.build_files(args, paths, len(args.nodes))
            for i, node in enumerate(args.nodes):
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
                time.sleep(30)
            io.mkdir('geoms')
            for i, node in enumerate(args.nodes):
                for geom in os.listdir(node + '/geoms'):
                    j += 1
                    shutil.copy(node + '/geoms/' + geom, 'geoms/{}_{}.xyz'.format(geom.split('_')[0], str(j).zfill(2)))
                for geom in os.listdir(node + '/output'):
                    if 'opt_' in geom:
                        k += 1
                        shutil.copy(node + '/output/' + geom, 'output/{}opt_{}.out'.format(geom.split('opt')[0], str(k).zfill(2)))
            samps = j
            args.restart = 1
            args.jobs = alljobs
            for i in range(len(args.reacs)):
                filename = es.check_geoms(paths['qtc'], 'reac' + str(i+1), samps)
                filename = filename.split('/')[1].split('_')[0] + '_opt_' +  filename.split('_')[1]
                filename = 'output/' + filename.replace('.xyz','.out')
                shutil.copy(filename, 'output/reac' + str(i+1) + '_opt.out')
            for i in range(len(args.prods)):
                filename = es.check_geoms(paths['qtc'], 'prod' + str(i+1), samps)
                filename = filename.split('/')[1].split('_')[0] + '_opt_' +  filename.split('_')[1]
                filename = 'output/' + filename.replace('.xyz','.out')
                shutil.copy(filename, 'output/prod' + str(i+1) + '_opt.out')
        stoichs, symnums = es.build_files(args, paths)
        es.execute(paths, args.nodes[0])
    if not samps:
        samps = 1000
    #check for failures
    if ("Opt" in args.jobs and not "Opt_1" in args.jobs):
        for i in range(len(args.reacs)):
            filename = es.check_geoms(paths['qtc'], 'reac' + str(i+1), samps)
            filename = filename.split('/')[1].split('_')[0] + '_opt_' +  filename.split('_')[1]
            filename = 'output/' + filename.replace('.xyz','.out')
            shutil.copy(filename, 'output/reac' + str(i+1) + '_opt.out')
        for i in range(len(args.prods)):
            filename = es.check_geoms(paths['qtc'], 'prod' + str(i+1), samps)
            filename = filename.split('/')[1].split('_')[0] + '_opt_' +  filename.split('_')[1]
            filename = 'output/' + filename.replace('.xyz','.out')
            shutil.copy(filename, 'output/prod' + str(i+1) + '_opt.out')
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


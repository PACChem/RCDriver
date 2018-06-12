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

    if args.restart < 5:
        stoichs, symnums = es.build_files(args, paths)
        es.execute(paths, args.node)
     
    #check for failures
    if ("Opt" in args.jobs and not "Opt_1" in args.jobs):
        es.check_geoms(es.nsamps)
    if ("1dTau" in args.jobs or 'MdTau' in args.jobs):
        for i in range(len(args.reacs)):
            es.check_hrs(i+1,'reac')
        for i in range(len(args.prods)):
            es.check_hrs(i+1,'prod')
        es.me_file_abs_path()

    import results 
    rs = results.RESULTS(args, paths)
    args.hlen = rs.get_hlen()
    args.optlevel = rs.optlevel
    args.enlevel = rs.enlevel
    args.taulevel = rs.taulevel

    #######  Build and run thermo  #########
    ########################################
    import thermo
    rs.thermo = False
    if args.alltherm.lower() == 'true':
        rs.thermo = True
        args.symnums = symnums
        rs.dH0, rs.dH298, rs.hfbasis, rs.hfbases = thermo.run(args, paths)

    #######  Parse results  #########
    ########################################
    if args.parseall.lower() == 'true':
         rs.get_results()

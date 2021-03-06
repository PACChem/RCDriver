# RCDriver
[//]: # (Badges)
[![CircleCI](https://circleci.com/gh/PACChem/RCDriver.svg?style=svg)](https://circleci.com/gh/PACChem/RCDriver)
[![Anaconda-Server Badge](https://anaconda.org/pacchem/rcdriver/badges/version.svg)](https://anaconda.org/pacchem/rcdriver)
[![Anaconda-Server Badge](https://anaconda.org/pacchem/rcdriver/badges/platforms.svg)](https://anaconda.org/pacchem/rcdriver)
[![Anaconda-Server Badge](https://anaconda.org/pacchem/rcdriver/badges/installer/conda.svg)](https://conda.anaconda.org/pacchem)

Sarah N. Elliott, Murat Keceli, and Stephen J. Klippenstein

RCDriver is a set of python modules that when given a short input file (that requests a SMILES molecule name and levels of theory) will perform hindered rotor scans and transition state searches.  Specifically it sets up and performs EStokTP computations and subsequently parses the resulting output files to compute 0 K heats of formation and updated anharmonic constants.  Given the option, RCDriver will use those computations to generate mess input and can run the mess partition 
function, thermp, pac99 executables to generate 298 K heats of formation, heat capacities, and NASA polynomials. 



## (1) GETTING TORSSCAN

Users can clone RCDriver from https://github.com/snelliott/RCDriver or simply use its location at  /home/elliott/Packages/RCDriver on Blues.  


## (2) DEPENDENCIES

RCDriver relies on EStokTP developed by Carlo Cavallotti, Matteo Pelucchi, and Stephen Klippenstein. RCDriver uses the iotools (input/output tools), obtools (openbabel tools), and patools (parsing tools), built by Murat Keceli and Sarah Elliott that can be cloned from  https://github.com/keceli/QTC or https://github.com/PACChem/QTC. Obtools, in turn, needs OpenBabel with pybel python bindings.  You will need to either install that by following these instructions https://pypi.python.org/pypi/openbabel.
The zmat builder uses x2z by Yuri Georgievski.  This can be cloned from https://github.com/PACChem/x2z. 
The thermochemistry computations in RCDriver uses heatform, anharm, and tctools also on the QTC github, pac99 by Bonnie Mcbride, thermp by Stephen Klippenstein, and mess by Yuri Georgievski. The locations of these files should be in the user's path or specified in a configfile (run rcdriver with -c <configfile>)


## (3) INPUT / OUTPUT

The main executable is rc_driver.py.  It will need an inputfile specified with the -i flag (an example is located at /home/elliott/Packages/RCDriver/input.dat).  If no inputfile is specified with the (e.g., the command is not  rc_driver.py myinputfile.txt)  the code will automatically look for a file named input.dat. The input file separates keywords from their input values with a colon.  The keywords are case sensitive but the values are not.  Lines can be commented out with a \#.  The output will print to terminal unless an output file is specified with a -o flag.


### SPECIES INPUT

**Reactant list:**

List up to two, comma-separated SMILES strings.  OpenBabel will generate an xyz geometry for x2z to transform into the zmat found in reac1.dat and reac2.dat unless a cartesian file is specified (see USE INPUT XYZ section). Multiplicity can be assigned by appending _m<mult> to the SMILES string.

*Reactant list (SMILES): C, [OH]_m2*



**Product list:**

List up to two, comma-separated SMILES strings 

*Product list:  [CH3], O_m1*



**Reaction type:**

Specify the type of reaction.  Currently RCDriver will perform Abstraction, Addition, and Isomerization reactions.
Default is blank, which means that the job is a well.

*Reaction type: Abstraction*



**No. of transition states:**

Specify the number of transition states. 0 and 1 are self explanatory, 2 and 3 will generate van der Waals wells. Default is 0.  When running a transition state, EStokTP needs to know the i,j,k sites.  Either have an rmg.dat file for reac1 in your working directory or add an xyz or geo file (and set Use QTC xyz in the next section to True)  for reac1 in your working directory and edit it so that there is a 1 and a space before the atom that will be losing the hydrogen, a 2 and a space before the hydrogen that will be abstracted, and a 3 and a space before a different atom that is connected to the atom that is losing the hydrogen.  If the second reactant is O, OH, H, HO2, or O2, you do not need an xyz or geo file for reac2.  Otherwise provide a xyz or geo for reac2 (labeled <SMILES>.xyz) and put a 4 and a space before the atom that will abstract the hydrogen OR put that line at the top of the geometry and a 4 will not be needed. An example of these files is provided at the end of this document.

*No. of transition states: 2*


### GEOMETRY OPTIONS

**Use input xyz:**

A user can choose to use an xyz or geo file instead of having OpenBabel generate it from the SMILES string by setting this keyword to True, cartesianfile.xyz, or logfilename.log. Using true will mean that the xyz or geo are in the working directory and are named <SMILES>.xyz or <SMILES>.geo (geo files do not have the extra two lines at the top with the first having the number of atoms and the second being blank/comment). Using <logfilename>.log will make RCDriver parse out the coordinates from a Molpro or Gaussian job (Note that using RCDriver for more than one species will break this method until I update it to use <SMILES>.log instead of <logfilename>.log). Default is false.

*Use input xyz: True*


**Use xyz as:**

EstokTP optimizes geometries in multiple of its modules.  This keyword tells RCDriver which of these modules the geometry from the previous keyword should replace.  Using start will use the geometry instead of generating an OpenBabel geometry and start EstokTP from the beginning (i.e., from level0 aka Opt_reac1).  Using 0 or level0 will use the geometry instead of optimizing at level0 so EstokTP will begin at level1 (aka Opt_reac1_1).  Finally 1 or level1 will use the geometry as the level1 optimization and skip straight to the hindered rotor scans BUT this is not quite working yet because it needs to also get force constant information in the output directory. Default is start.

*Use xyz as: 0*


### BLUES OPTIONS

**Run on node:**

This lets the user specify which blues node to run EstokTP on.  If the input is d or debug, RCDriver will build all the necessary EstokTP files but not run them.  Specifying b###  will submit the job to that blues node.  A comma seperated list of b###,b###,b### will submit the monte carlo level0 job to multiple nodes, and the remaining modules will be sent to teh first in the list. Using 0 will run EStokTP on your current node (on the login node if you haven’t sshed onto a blues node). 
Default is 0.

*Run on node: b431*


**No. of cores high:**

Specify the number of cores for your high level EStokTP computation to use. Default is 16.

*No. of  cores high: 20*


**No. of cores low:**

Specify the number of cores for your low level EstokTP computations. Default is 10.

*No. of cores low: 16*


**Memory high:**

Specify the memory for high level computations in MW. Default is 200.

*Memory: 500*

**Memory low:**

Specify the memory for low level computations in MW. Default is 200.

*Memory: 500*

### EStokTP OPTIONS

**No. of MC sampling points:**

Give the number of points for the Monte Carlo sampling (The level0 Opt).  Default is 5.

*No. MC sampling points: 5*



**Calculate no. of MC points:**

A, B, C, D input to compute A+B*numberOfTorsions^C a select the minimum of that and D.

*No. MC sampling points: 3, 1, 3, 100*


NOTE: if both MC point keywords are used, No. of MC sampling points takes priority.



**Scan interval:**

Specify the number of degrees for the hindered rotor scan to scan. Default is 360.

*Scan interval (degrees): 360*


**No. of steps on the PES:**

Specify the max number of points to take on the PES while doing a hindered rotor scan  (it will be divided by symmetry number). Default is 4.

*No. of steps on the PES: 8*


**Multidim scan:**

If a multidimensional scan (Mdtau) is specified in the module input section, this will tell it whether it should be a 2D or 3D scan by entering: 2, 2D, 3, 3D, or auto.  Default is 2.

*Multidim scan (2 or 3D): 2*



**Module table:**

The first column is a list of all of the EStokTP module names.  You can delete rows that you don’t need if desired, but don’t edit the module names that you will use.  If the program or level of theory isn’t specified in the next two columns (columns need to be separated by colons and spacing doesn’t actually matter) it will automatically skip those modules.  In order to run modules lower on the list, the modules above it have to successfully complete first.  The available modules are Opt (the level0 optimization and MC sampling, if you’ve used Geometry options to use xyz as level0 this will just be ignored), Opt_1 (the level1 optimization and frequency computation), 1dtau (The one dimensional hindered rotor scan), MdTau (the multidimensional hindered rotor scan – yes, 1dtau has to be run first), Symm (determines the rotational and optical symmetry numbers – BROKEN right now because it can’t find a MESS executable it needs), HL (high level optimization. TODO: set it up to use correlation additivity i.e., ccsd(t)/tz ~ ccsd(t)/dz + ccsd/tz – ccsd/tz, which should not be too hard for molpro), Irc (an irc), and kTP (which will set up the mess files and run it to get your rates if there is a reaction). The Program has to be either blank, g09, gaussian, or molpro.  Theory should be specified with method/basis.  I’ve found that EStokTP has trouble parsing the energies from completed computations for some methods and basis sets depending on which program you are using so be sure to check the output/estoktp.out file to see if your job died because of the method/basis selection.  


### THERMO OPTIONS

**Perform all thermochemistry?:**

True or false to run heat of formation code and if me_files are produced during EStokTP computation (the 1dTau has successfully completed) then it will execute mess’s partition function code, thermp, pac99, and produce a chemkin file.  The resulting files will have the name <SMILES>.pf (pf input) <SMILES>.dat (pf output), thermp.dat, thermp.out, <SMILES>.i97, <SMILES>.o97,  and <SMILES>.ckin. Default is true.

*Perform all thermochemistry? true*


**Anharmonic:**

Here we tell RCDriver if it should be using anharmonic corrected frequencies and anharmonic constants in the mess input file for partition function.  There are several options to tell it how to get the anharmonic constants.  Saying 0 will run an anharmonic frequency calculation at the same level of theory as the level0 optimization and saying 1 will do the same but for level1 optimization. Specifying prog/method/basis will have it run it at a specific program, method, and basis set.  All of these options are only set up for Gaussian computations at the moment. Inputing <SMILES>anharm.log will tell RCDriver that you’ve seperately run an anharmonic computation, placed it in the working directory, and named it <SMILES>anharm.log, and want RCDriver to parse the anharmonic constants from that file. The final way is to write optprog/optmethod/optbasis/prog/method/basis which will make it search in the PACC for anharmonic constants calculated at prog/method/basis for a geometry optimized at optprog/optmethod/optbasis.  Default is false.

*Anharmonic: COanharm.log*


**Overwrite anharmonic:**

If RCDriver already sees a <SMILES>anharm.log file it will just use that instead of computing on for itself.  If that’s not okay with you, input true to overwrite. Default is False.

*Overwrite anharmonic: true*


**Basis for heat of formation:**

This will be used in the heat of formation code.   The default is auto which will let the heat of formation determine the basis itself.  Otherwise the user can list molecules written as SMILES strings and separated by spaces (not commas, but if this gives enough problems, let me know and I can easily let let it read either option).  Be aware that if the user specified basis produces  a singular matrix during the heat of formation computation, the code will choose a new basis.  By basis, we mean the molecules whose heats of formations, electronic energies, and zero-point vibrational energies,  will be used in a linear combination to estimate the heat of formation for our species of interest. The coefficients for this combination are C: the coefficients [a,b,c,d..] for the basis molecules  [R1, R2, R3, R4...] that give CxHyOzNw =  aR1 + bR2 + cR3 + dR4 where [x,y,z,w...] is our stoichiometry vector for the molecule of interest.  See heat of formation documentation for more information.

*Basis for heat of formation: [O][O] [H][H] C*



### MISCELLANEOUS OPTIONS

**Restart at:**

Did something break in your computation and you don’t want to run it all all over again? No worries!  You can tell it to restart at almost any point. 0 – beginning (the default), 1 – if level0 optimization has successfully completed and you want to begin with level1, 2 – if level1 is completed and you want to start at the next module, 3 – The 1D hindered rotor scans are completed, 4 - the MD hindered rotor scans are completed (i.e., just symm, high level, irc, and kTP are left), 5 - just run thermo and no EStokTP computation.  If you’ve specified in geometry options that you want to use an xyz as the level0 geometry, then RCDriver will automatically set the restart to 1 if it is not at a higher restart level already.
Not enough options for you? If you have multiple reactants and/or products and want to specifically start at, say, a level1 computation for your second product you’ll have to go into data/EStokTP and delete or add a character to all of the modules prior to Opt_Reac2_1.  Then you have to run EStokTP directly instead of submitting it through RCDriver (which will overwrite those edits) by using the command /lcrc/projects/PACC/projects/codes/EstokTP/exe/estoktp.x &> estoktp.log.  

*Restart at: 1*


**Parse all:**

True will command RCDriver to parse all of the important quantum chemistry (and thermochemistry if Perform all thermochemistry? True)  properties, store them in the PACC database, and print them out at the end of the computation for all the species you’ve run the computation on.  There is not currently a way to parse partial.  

*Parse all: True*




### Example transition state input xyz files:

**In input.dat:**

Reactants list: CC, [CH3]

Use QTC xyz: True

Use xyz as: start


**In CC.xyz:**

8

1 C          1.09278       -0.02695       -0.03725

3 C          2.60483       -0.02695       -0.03725

2 H          0.70841        0.21668        0.95770

   H          0.70841       -1.01041       -0.32373
   
   H          0.70841        0.71289       -0.74571
   
   H          2.98920       -0.76678        0.67122
   
   H          2.98920        0.95652        0.24924
   
   H          2.98920       -0.27058       -1.03220

**In [CH3].geo:**

   C          1.15421        0.14357        0.13126
   
   H          0.64459       -0.52593        0.81229
   
   H          0.64457        1.01959       -0.24933
   
   H          2.17554       -0.05966       -0.16467
   
or

  H          0.64459       -0.52593        0.81229
  
4 C          1.15421        0.14357        0.13126

   H          0.64457        1.01959       -0.24933
   
   H          2.17554       -0.05966       -0.16467
   




**My job didn’t complete. What happened?**

So you followed all of the instructions and have a good input file but something still broke.  We have a few places to check to see what went wrong.  First I advise looking at the RCDriver output.  Did it break in the reac1-prod2.dat builds? Did you tell it you would provide an input file that you didn’t (<SMILES>.xyz, geo, anharm.log etc)? Make sure its there. Also, check to make sure your SMILES name is correct.  If it made it past that point: did it print your error and immediately stop? Maybe it was missing an executable (it should print that if so). Check the dependencies listed in section 3. Did it break when running mess?   Check to see if the me_files directory is populated. If that’s all good make sure you have the right executables for the thermochemistry computations. If it didn’t write the me_files, that means probably something went wrong in the EStokTP computation. Let’s check there… 
First go to output/estoktp.out and search for ‘starting’.  What was the last thing it started?  Nothing? Then go to estoktp.log and see if it tells you it can’t find an executable or specifies another error.  Otherwise was the last thing started in output/estoktp.out an optimization? Then check geom.log (gaussian job) or molpro.log (molpro job) and see why the optimization failed.  Edit your input file accordingly and restart at the right point (see Restart at).  Does it look the the optimization was complete but it didn’t move on to the next step? Maybe it had trouble parsing your output – check to see if output/estoktp.out says it couldn’t find a keyword – it has trouble with some methods and basis sets.  Still haven’t found your issue?  Keep perusing the data/estoktp.out and estoktp.log and if you still can’t figure it out you may have to go into the data directory and look at those files. You’ll need some knowledge about how EStokTP works to use these, however.

**Want to know more about what RCDriver is doing?**

The majority of the files RCDriver generates are in the data directory.  To understand the contents of these files please view the EStokTP manual (put link here when available).  All of the python files in RCDriver (rc_driver.py, build.py, anharm.py, and /home/elliott/Packages/QTC/heatform.py) are somewhat commented so feel free to dig around!  You can direct questions to Sarah Elliott.

## (4) ACKNOWLEDGEMENT
This research was supported by the Exascale Computing Project (ECP), Project Number: 17-SC-20-SC -- a collaborative effort of two DOE organizations, the Office of Science and the National Nuclear Security Administration, responsible for the planning and preparation of a capable exascale ecosystem including software, applications, hardware, advanced system engineering, and early test bed platforms to support the nation's exascale computing imperative -- and by the DOE Krell Computational Science Graduate Fellowship Grant Number DE-FG02-97ER25308. 


## (5) COPYRIGHT
Copyright 2018 Sarah N Elliott, Murat Keceli, and Stephen J Klippenstein

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

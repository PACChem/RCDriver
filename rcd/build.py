import os
import numpy as np
import sys
import logging
from qtc import obtools as ob
from qtc import iotools as io
from qtc import patools as pa
from qtc import qctools as qc
log = logging.getLogger(__name__)

class MOL:
    def __init__(self,paths, opts,typemol = 'reac', reactype = ''):
        """
        MOL object is reac, prod, or ts.
        """ 
        self.typemol    = typemol
        self.paths      = paths

        #OPTIONS##############################
        self.nsamps     = opts[0]    #number of MonteCarlo sampling points
        self.abcd       = opts[1]
        self.nodes      = opts[2]
        self.interval   = opts[3]    #Interval for tors scan (should be same geometry as 0 degree)
        self.nsteps     = opts[4]    #Number of points to take on PES
        self.XYZ        = opts[5]    #if QTC provides XYZ, 'true' (smiles.xyz), logfile 'anything.log', or just use smiles 'false'
        self.xyzstart   = opts[6]
        self.MDTAU      = opts[7]
        ######################################
        self.ijk        = [0, 0, 0, 0]
        self.sort       = None
        self.reactype   = reactype
 
    def build_xyzfile(self, smiles,findxyz=False):

        """
        Finds user-specified xyz or geo coordinates or generates openbabel coordinates and places
        them in data/<smiles>.xyz for x2z to use
 
        INPUT:
        smiles         -- SMILES name for molecule
        OUTPUT:
        smilesfilename -- name of cartesian coordinate file
        """
        cartlines = ''
        found = True
        XYZ = self.XYZ
        if findxyz:
            XYZ = 'true'
        if 'ts' in smiles and 'geomdir' in XYZ.lower():
            cartlines, ijk, found = read_cart('geomdir/' + smiles, None, self.reactype)
            smilesfilename = 'ts'
        else:
            smilesfilename = ob.get_smiles_filename(smiles)
            if '_m' in smiles:
                smiles, self.mult = smiles.split('_m')
            else:
                self.mult = ob.get_multiplicity(ob.get_mol(smiles))

            if XYZ.lower() == 'false':
                cartlines = build_obcart(smiles,self.mult)
                found = False

            elif '.log' in XYZ.lower():
                cartlines = io.read_file('../' + XYZ)
                cartlines = pa.gaussian_xyz_foresk(cartlines)

            elif '.xyz' in XYZ.lower():
                if io.check_file(XYZ):
                    cartlines = io.read_file(XYZ)
                else:
                    log.warning('no .geo or .xyz provided\n...Using openbabel instead')
                    cartlines = build_obcart(smiles, self.mult)
                    found = False
            elif len(XYZ.split('/') ) < 2:
                if 'geomdir' in XYZ.lower():
                    cartlines, ijk, found = read_cart('geomdir/' + smiles, self.mult, self.reactype)
                else:
                    cartlines, ijk, found = read_cart(smiles, self.mult, self.reactype)
                if ijk[0] != 0:
                    self.ijk = ijk
            elif len(XYZ.split('/')) > 2:  
                XYZ = self.XYZ.replace('g09','gaussian')

                if io.check_file(io.db_opt_path(XYZ.split('/')[0],XYZ.split('/')[1],XYZ.split('/')[2], None, smiles) + '/' + smilesfilename + '.geo'):
                    cartlines = io.db_get_opt_prop(smiles, 'geo', None, XYZ.split('/')[0], XYZ.split('/')[1], XYZ.split('/')[2])
                    cartlines = cartlines.replace('\n\n','')
                    cartlines = str(len(cartlines.split('\n'))) + ' \n\n' + cartlines
                elif io.check_file(io.db_opt_path(XYZ.split('/')[0], XYZ.split('/')[1], XYZ.split('/')[2], None, smiles) + '/' + smilesfilename + '.xyz'):
                    cartlines = io.db_get_opt_prop(smiles, 'xyz', None, XYZ.split('/')[0], XYZ.split('/')[1], XYZ.split('/')[2]).split('\n\n')[1]
                    cartlines = cartlines.split('\n')[0] + ' \n\n' + cartlines

                else:
                    log.warning('No geometry found at ' + io.db_opt_path(XYZ.split('/')[0], XYZ.split('/')[1], XYZ.split('/')[2], None, smiles)
                                   + smilesfilename + '.xyz' + '\n...building using OpenBabel instead')
                    cartlines = build_obcart(smiles,self.mult)
                    found = False
            else:
                log.error('You have not specified a valid way to get the coordinates.  Use false, true, smiles.log, smiles.geo, smiles.xyz, or prog/method/basis')
                found = False
        io.write_file(cartlines,smilesfilename + '.xyz')
        return smilesfilename, found

    def cart2zmat(self,smiles, select = []): 
        """
        `Runs x2z by Yuri Georgievski on a file of Cartesian coordinates and collects
        the internal coordinate and torsional angle information that it outputs
        """
        import re

        ######   initialize  ###############
        ####################################
        atoms   = []
        measure = []
        angles  = []
        consts  = []

        smilesfilename, found = self.build_xyzfile(smiles)
        if smiles != 'ts':
            mol         = ob.get_mol(smiles,make3D=True)
            self.charge = ob.get_charge(mol)
            self.stoich = ob.get_formula(mol)
        if not found and self.XYZ == 'geomdir':
            smilesfilename, found = self.build_xyzfile(smiles,True)
            found = False
        if found or smiles != 'ts':
            #######   RUN X2Z   #################
            #####################################
            tempfile = 'temp'
            paths = self.paths
            gcc     = paths['gcc']
            intel   = paths['intel']

            if io.check_file(smilesfilename + '.xyz'):
                os.system('{0}; {1}; {2} {3}.xyz > {4}'.format(gcc, intel, 'x2z', smilesfilename,tempfile))
            else:
                os.system('{0}; {1}; {2} {3}_m{4}.xyz >  {5}'.format(gcc, intel, 'x2z', smilesfilename, str(self.mult),tempfile))

            if os.stat(tempfile).st_size < 80 or 'not connected' in io.read_file(tempfile):
                log.warning('Failed')
                log.warning('Please check that directory name and cartesian coordinate file name are equivalent')
                log.warning('Please check that test_chem is in location: ' +  'x2z')
                log.warning('Using OpenBabel zmat: no hindered rotors will be specified')
                #Build openbabel zmat if x2z fails
                atoms, measure = build_obzmat(smiles) 
                self.symnum    = ' 1'
                self.ilin      = ' 0'

                return atoms, measure, angles, found, ''

            ####Get relevant data from x2z output file########
            ########################################################
            props, lines = io.read_file(tempfile).split('Z-Matrix:\n')
            io.rm(tempfile)
            
            #linearity
            self.ilin =' 1'
            result = re.search("molecule is (\w+)", props).groups()[0]
            if ('non' in result or 'plan' in result):
                self.ilin = ' 0'

            #enantiomers
            ent = False
            result = re.search("has enantiomer?(\w+)", props)
            if result: 
                result = result.groups()[0]
                if ('yes' in result):
                    ent = True

            #Symmetry number
            self.symnum = re.search("symmetry number = (\w+)", props).groups()[0]#Symmetry factor

            #Rearrangement
            props, order = props.split('Z-Matrix atom order:')
            self.sort = {}
            for index in order.split('\n')[1:-2]:
                index = index.replace('X','0')
                self.sort[int(index.split('>')[1])] = str(index.split('-->')[0])
            #rotational groups
            groups = ''
            if len(lines.split("Rotational groups:")) > 1:
                lines, groups = lines.split("Rotational groups:") #if x2z prints rot groups

            #zmatrix connectivity
            lines = lines.split('\n')
            twoorthree = 2
            if self.ilin == ' 1':
               twoorthree = 3
            twoorthree = 3 #seems to need three coordinates every time?
            for i,line in enumerate(lines):
                if len(consts)%twoorthree != 0:
                    consts.append(line.split(',')[-1].strip())
                if 'x' in line.lower():
                    consts.append(line.split(',')[2].strip())
                    if len(line.split(',')) > 4:
                        consts.append(line.split(',')[4].strip())
                if line == '':
                    break
                atoms.extend([line.rstrip('\n').replace(' ','').split(',')])     #Main part of ZMAT
            #zmatrix parameters
            for j in range(i+1,len(lines)):
                if lines[j]:
                    if 'const' in lines[j].lower() or 'Rot' in lines[j] or 'Beta' in lines[j]:
                        break
                measure.extend(lines[j].replace('=',' ').rstrip('\n').split())   #Gets parameters

            #Hindered rotor dihedral angles and constant angles
            for n in range(len(lines[j:])):
                if "Rot" in lines[j+n] and lines[j+n].split(':')[1].rstrip('\n').strip() != '':
                        angles  = lines[j+n].replace(" ","").upper().split(':')[1].strip().split(',')  
            
            #Reformat zmatrix parameters
            measure = np.array(measure)                     
            if  (len(measure)%2 != 0):
                measure = measure[:-1]
            measure = measure.reshape( len(measure)/2, 2)      #Puts measurements into two columns
            for angle in measure:
                if 'R' in angle[0]:
                    angle[1] = str(float(angle[1]) * 0.529177) #bohr to angstrom
            #Put dummy atoms paramters inside zmat
            #for i, row in enumerate(atoms):
            #    if 'X' in row[0]:
            #        j=0
            #        for meas in measure:
            #            if  meas[0] == row[2]:
            #                atoms[i][2] = '{:.6f}'.format(float(meas[1]))
            #                measure = np.delete(measure, j, axis=0)
            #                j-=1
            #            if len(row) > 3:
            #                if meas[0] == row[4]:
            #                    atoms[i][4] = '{:.2f}'.format(float(meas[1]))
            #                    measure = np.delete(measure, j, axis=0)
            #                    j-=1
            #            if len(row) > 5:
            #                 if meas[0] == row[6]:
            #                    atoms[i][6] = '{:.2f}'.format(float(meas[1]))
            #                    measure = np.delete(measure, j, axis=0)
            #                    j-=1
            #            j+=1
            #Put constant angles inside zmat 
            for angle in consts:
                for i, row in enumerate(atoms):
                    if len(row) > 6:
                        if angle == row[6]:
                            j=0
                            for meas in measure:
                                if  meas[0] == angle:
                                    atoms[i][6] = '{:.4f}'.format(float(meas[1]))
                                    measure = np.delete(measure, j, axis=0)
                                    j-=1
                                j+=1
                    if len(row) > 4:
                        if angle == row[4]:
                            j=0
                            for meas in measure:
                                if  meas[0] == angle:
                                    atoms[i][4] = '{:.4f}'.format(float(meas[1]))
                                    measure = np.delete(measure, j, axis=0)
                                    j-=1
                                j+=1
                    if len(row) > 2:
                        if angle == row[2]:
                            j=0
                            for meas in measure:
                                if  meas[0] == angle:
                                    atoms[i][2] = '{:.6f}'.format(float(meas[1]))
                                    measure = np.delete(measure, j, axis=0)
                                    j-=1
                                j+=1
            #Count rotors
            nmethylgroups = 0 
            rotors = {}
            groups = groups.split('Beta')[0]
            groups = groups.split('no')[0]
            for rotor in groups.split('\n'):
                if rotor:
                    grouplist = rotor.lower().split()[1:]
                    if 'c1h3' in grouplist:
                       nmethylgroups += 1
                       rotors[rotor.split()[0]] = [0, grouplist]
                    else:
                       rotors[rotor.split()[0]] = [1, grouplist]
            ##reorder angles to have nonmethyl rotors first
            methylrotors = []
            nonmethylrotors = []
            for rotor in rotors:
                if rotors[rotor][0] == 1:
                    nonmethylrotors.append(rotor)
                else:
                    methylrotors.append(rotor)
            angles = nonmethylrotors + methylrotors
            for i, angle in enumerate(angles):
                 if str(i+1) in select or angle.lower() in select:
                     del angles[i]
                     angles.insert(0, angle)
       
            msg = '' 
            if angles:
                msg += '{}\n  Num\tLabel\tgroup1\tgroup2\n'.format(smiles)
                for i, angle in enumerate(angles):
                    msg += '  {:g}\t{}\t{}\t{}\n'.format(i+1, angle, rotors[angle][1][0], rotors[angle][1][1])
            
            self.nrotors = len(angles)
            self.nrotors = self.nrotors - nmethylgroups
            
            #diatomics and triatomics must have 0 for ilin in EStokTP
            if len(atoms) < 5:                               
                self.ilin  =' 0' 
        return atoms, measure, angles, found, msg
         

    def build(self, n, smiles, jobs, found, angles, atoms = [], measure = [], restartts=False):
        """ 
        Builds reacn.dat or prodn.dat for EStokTP withh user defined nosmps (Monte Carlo sampling points
        for geometry search) and nhindsteps (number of points on the PES) 
        """
        if   n == 3:
            self.typemol = 'wellr'
        elif n == 4:
            self.typemol = 'wellp'
        if self.typemol == 'reac' or self.typemol == 'prod' or 'well' in self.typemol:
            if len(atoms) < 3 and self.ijk[0] and n == 1:
                atoms.append(list(atoms[1]))
                atoms[1][0] = 'X'
                atoms[1][2] = '2.0'
                atoms[2].extend(['2', '90.0'])
                #print self.ijk
                #for o in range(len(self.ijk)):
                #    if self.ijk[o] > 0:
                #       self.ijk[o] = str(int(self.ijk[o]) + 1 )
                self.ijk[2] = 'X2'
                #print self.ijk
            atoms, measure, angles  = update_interns(n,atoms,measure,angles)
            #Torsional Scan Parameters#############
            taustring, self.nrotors = tau_hind_str([atoms], jobs, [angles], self.interval, self.nsteps, self.MDTAU)
            if not self.nsamps:
                if len(self.abcd.split(',')) >3:
                    a, b, c, d = self.abcd.split(',')
                    a = int(a)
                    b = int(b)
                    c = int(c)
                    d = int(d)
                    self.nsamps = min(a + b * c**self.nrotors, d)
            self.nsamps = str(int(np.ceil(float(self.nsamps) / self.nodes)))
            #MC Parameters#############
            zmatstring  = 'nosmp dthresh ethresh\n'     
            zmatstring += self.nsamps + '  1.0  0.00001\n' + taustring
            #Size and linearity of molecule###########
            zmatstring += '\nnatom natomt ilin\n'
            ndummy = count_dummy(atoms)
            zmatstring += str(len(atoms)-ndummy) + ' ' + str(len(atoms)) + self.ilin + '\n'
            if n==1 and self.ijk[2]:
                if len(atoms) < 5 and ndummy > 0:
                    for o, row in enumerate(atoms):
                        if 'x' in row[0]:
                            self.ijk[2] = 'X' + str(o+1)
                elif int(self.ijk[2])-1 in self.sort:
                #if len(atoms) < 7 and ndummy > 1:
                    for col in atoms[int(self.sort[int(self.ijk[2])-1])-1]:
                        if 'x' in col.lower():
                            self.ijk[2] = 'X' + col.split('x')[1]
                            break
                if int(self.ijk[0])-1 in self.sort:
                    inchain = False
                    hasdumm = False
                    dumm    = ''
                    for col in atoms[int(self.sort[int(self.ijk[0])-1])-1]:
                        try:
                            if int(col[1:]) == int(self.sort[int(self.ijk[1])-1]):
                                inchain = True
                            if 'x' in col.lower(): 
                                hasdumm = True
                                dumm = 'X' + col.split('x')[1]
                        except:
                            pass 
                    if inchain and hasdumm:
                        self.ijk[2] = dumm
        else:
            #Torsional Scan Parameters#############
            if found: 
                atoms, measure, angles  = update_interns(1,atoms,measure,angles)
                if n == 'ts' and self.babs and 'isomerization' not in smiles.lower():
                    taustring, self.nrotors = tau_hind_str([atoms], jobs, [angles], self.interval, self.nsteps, self.MDTAU, self.ijk, self.babs)
                elif n!= 'wellr' and n!= 'wellp':
                    taustring, self.nrotors = tau_hind_str([atoms], jobs, [angles], self.interval, self.nsteps, self.MDTAU)
                else:
                    taustring, self.nrotors = tau_hind_str([], jobs, [], self.interval, self.nsteps, self.MDTAU)
            else:
                if n == 'ts' and self.babs and 'isomerization' not in smiles.lower():
                    taustring, self.nrotors = tau_hind_str(atoms, jobs, angles, self.interval, self.nsteps, self.MDTAU, self.ijk, self.babs)
                elif n != 'wellr' and n!= 'wellp':
                    taustring, self.nrotors = tau_hind_str(atoms, jobs, angles, self.interval, self.nsteps, self.MDTAU)
                else:
                    taustring, self.nrotors = tau_hind_str([], jobs, [], self.interval, self.nsteps, self.MDTAU)
            if not self.nsamps:
                if len(self.abcd.split(',')) >3:
                    a, b, c, d = self.abcd.split(',')
                    a = int(a)
                    b = int(b)
                    c = int(c)
                    d = int(d)
                    self.nsamps = min(a + b * c**1, d)
            self.nsamps = str(int(np.ceil(float(self.nsamps) / self.nodes)))
            #MC Parameters#############
            zmatstring  = 'nosmp dthresh ethresh\n'     
            zmatstring += self.nsamps + '  1.0  0.00001\n' + taustring
            if n == 'ts':
                import get_sites
                #i,j,k sites###########################
                if 'isomerization' in smiles.lower():
                    zmatstring += '\nisite ji ki\n'
                else:
                    zmatstring += '\nisite jsite ksite\n'

                if self.ijk[0] != 0:
                    zmatstring += ' '.join(self.ijk[:3])
                else:
                    lines = io.read_file('../rmg.dat')
                    zmatstring += ' '.join(get_sites.sites(lines))
                bond = ''.join(sorted(self.bond.lower()))
                bondlengths = {'cc':1.54, 'ch':1.09, 'hh':.74,'nn':1.45,'oo':1.48,'cn':1.47,'co':1.43,'ho':1.2,'hn':.99}
                if 'isomerization' in smiles.lower():
                    rmin1 = '2.0'
                    rmin2 = '1.0'
                    rstep1 = '.2'
                    rstep2 = '.1'
                    zmatstring += '\n\nrmin1 rstp1 rmin2 rstp2 ireact2\n {} {} {} {} {}\n'.format(rmin1,rstep1,rmin2,rstep2,self.ijk[3])
                else:
                    rmin = '1.0'
                    rmax = '2.5'
                    nr   = '8'
                    babs1 = '90.0'
                    if restartts:
                        babs1 = '180.0'
                    if smiles.lower() == 'abstraction':
                        rmin = '0.7'
                        rmax = '2.2'
                        if bond in bondlengths:
                            rmin = str(bondlengths[bond])    
                            rmax = str(bondlengths[bond] + 1.0)    
                        babs1 = '180.0' 
                        if restartts:
                            babs1 = '90.0'
                    elif 'addition' in smiles.lower():
                        rmin = '1.4'
                        rmax = '2.8'
                        if bond in bondlengths:
                            rmin = str(bondlengths[bond] + 0.2)    
                            rmax = str(bondlengths[bond] + 1.6)    
                    zmatstring += '\n\nrmin rmax nr\n {} {} {}\n  -->aabs1,babs1,aabs2,babs2,babs3\n 85., {}, 85., 90., 90.\n'.format(rmin,rmax,nr,babs1)
            if found:
                zmatstring += '\nnatom natomt ilin\n'
                ndummy = count_dummy(atoms)
                zmatstring += str(len(atoms)-ndummy) + ' ' + str(len(atoms)) + self.ilin + '\n'
        #Typical Z-Matrix########################
        zmatstring += '\ncharge  spin  atomlabel\n'
        zmatstring += str(self.charge) + ' ' + str(self.mult) + '\n'
        if self.typemol == 'reac' or self.typemol == 'prod' or 'well' in self.typemol or found:
            for row in atoms:
                for j in range(len(row)):
                    zmatstring += row[j] + ' '
                zmatstring +='\n'
            zmatstring += '\nintcoor'

            deletedangles=[]
            for hin in angles:                                       #Deletes hindered angles from
                for i,meas in enumerate(measure):                    #The internal coordinate list
                    if hin.lower().strip() == meas[0].lower().strip():
                        deletedangles.append(measure[i][1]) 
                        measure = np.array([np.delete(measure.T[0],i),np.delete(measure.T[1],i)]).T
            for meas in measure:
                zmatstring += '\n' + meas[0] + ' ' + meas[1]
            
            zmatstring += '\n'
        #Sym factor and no. of electronic states#
        excited = {'[O]_m3':[[3,158.265],[3,  226.977]]}
        nelec = 1
        if smiles in excited:
            nelec += len(excited[smiles])
        zmatstring += '\nSymmetryFactor\n' +  self.symnum + '\n'
        zmatstring += '\nnelec\n{:d}\n 0.  {}'.format(nelec,str(self.mult))
        if smiles in excited:
            for level in excited[smiles]:
                zmatstring += '\n {:.3f}  {:d}'.format(level[1], level[0])
        zmatstring += '\n\nend\n'
        #Build Reac/Prodnum_opt.out for starting after level0 or level1
        if (self.xyzstart == '0' or self.xyzstart == '1') and smiles and found:
            if self.typemol == 'reac' or self.typemol == 'prod':
                smilesfilename = ob.get_smiles_filename(smiles)
                optim = build_optout(self.xyzstart, self.XYZ, measure, angles, deletedangles, smiles, smilesfilename)
                io.write_file(optim, '../output/' + self.typemol + str(n) + '_opt.out')
            else:
                optim = build_optout(self.xyzstart, self.XYZ, measure, angles, deletedangles, smiles, self.typemol + '.xyz')
                io.write_file(optim, '../output/' + self.typemol  + '_opt.out')


        return zmatstring

def build_molpro(meth,freqcalc,opt):
    """
    Builds molpro theory file 
    """
    method, basis = meth[2].split('/')
    nondft = ['hf', 'ccsd', 'cisd','mp']
    dft = True
    for key in nondft:
        if key in method.lower():
            dft = False
    molstr  = 'nosym\nEnd1\n\n'
    molstr += '!closed shell input\n\nbasis=' + basis +'\n'
    if dft:
        molstr += 'dft=['+method+']\n'
        molstr += 'hf\ndft'
    else:
        molstr += 'hf\n' + method.lower()
    if opt:
        molstr += '\noptg'
    if freqcalc:
        molstr += '\nfrequencies,symm=auto,numerical'
    molstr += '\nENERGY=energy\n\n'
    molstr += 'CBSen=energy\n\n'
    molstr += '! these lines must be always included in molpro input\n'
    molstr += '! CBSen should be defined as desired\n'
    molstr += '! the molden line should be left as it is\n\nput,molden,molpro.molden\n\n---\n\n'
    molstr += 'End2\n\n!open shell input\n\n'
    molstr += 'basis=' + basis +'\n'
    if dft:
        molstr += 'dft=['+method+']\n'
        molstr += 'uhf\ndft'
    elif 'ccsd' in method:
        molstr += 'rhf\n' + 'u' + method.lower()
    else:
        molstr += 'uhf\n' + method.lower()
    if opt:
        molstr += '\noptg'
    if freqcalc:
        molstr += '\nfrequencies,symm=auto,numerical'
    molstr += '\nENERGY=energy\n\n'
    molstr += 'CBSen=energy\n\n'
    molstr += '! these lines must be always included in molpro input\n'
    molstr += '! CBSen should be defined as desired\n'
    molstr += '! the molden line should be left as it is\n\nput,molden,molpro.molden\n\n---\n\n'
    molstr += 'End3\n\n\n\n'

    return molstr

def build_mehead():
    head = '!***************************************************\n!               GLOBAL SECTION\n!***************************************************\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!\n!\nTemperatureList[K]                      300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500  \nPressureList[atm]                       1.0     10.   100.  1000. 10000.\n!\n!\nEnergyStepOverTemperature               .2         ! [Discretization energy step (global relax matrix)] / T                   \nExcessEnergyOverTemperature             30         ! [Highest barrier in the model (global relax matrix)] / T                      \nModelEnergyLimit[kcal/mol]              400        ! Highest reference energy used in the calculation ( or ReferenceEnergy[kcal/mol])\n!\nCalculationMethod                       direct     ! direct or low-eigenvalue                     \n!\nWellCutoff                              20         ! well truncation parameter : Max { dissociation limit (min barrier rel. to bottom of the well) / T }\nChemicalEigenvalueMax                   0.2        ! Max chemical eigenvalue / Lowest Collision relaxation eigenvalue \n!\nReductionMethod                         diagonalization ! [low eigenvalue method only] diagonalization or projection (default)\n!      \n!!!!!!!!!test!!!!!!!!!!!!!!!!!!!!\n!WellCutoff                            10\n!ChemicalEigenvalueMin                 1.e-6          #only for direct diagonalization method\n!!!!!!!!test!!!!!!!!!!!!!!!!!!!!!!!!!!\nAtomDistanceMin[bohr]                  1.3\n!!\nRateOutput                              rate.out                        ! output file name for rate coefficients                         \n!\n!\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!***************************************************\n!               MODEL SECTION\n!***************************************************\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!\n!                          \nModel        \n!                                                           \n  EnergyRelaxation                                                      ! Default collisional energy relaxation kernel                              \n    Exponential                                                         ! Currently the only possible energy relaxation model              \n       Factor[1/cm]                     260                             ! (Delta_E_down)^(0) @ standard T (300 K)                          \n       Power                            0.875                           ! Power n in the expression (Delta_E_down) = (Delta_E_down)^(0) (T/T0)^(n)\n       ExponentCutoff                   10                              ! if (Delta_E) / (Delta_E_down) > value  transition probability is zero     \n    End  \n!                                                               \n'
    head += '  CollisionFrequency                                                    ! Collision frequency model\n    LennardJones                                                        ! Currently the only possible collisional frequency model  based on LJ potential\n       Epsilons[K]                       90.58  617.0                   ! Epsilon_1 and Epsilon_2 (630.4 x kB x Na = 1.25)(cm-1 to K = x 1.4) Ar and c7h7o2 (from Murakami)        \n       Sigmas[angstrom]                  3.54    5.62                   ! Sigma_1 and Sigma_2 (from Murakami) \n       Masses[amu]                       39.948 29.0                    ! Masses of the buffer gas molecule and of the complex (check order)\n    End      \n!\n!*************************************************\n!\n!***************************************************\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!***************************************************\n!  REACTANTS\t\t\t\t\t\n!***************************************************\n\n\n'
    return head 

def build_theory(meths,tss, zedoptions, oneoptions, adiabatic):
    """
    Builds theory.dat 
    meth[0] is module, meth[1] is program, and meth[2] is theory/basis
    """
    theory    = ''
    tsopt  = ' opt=(ts,calcfc,noeig,intern,maxcyc=50) \n '
    rzpopt = ' opt=(' + zedoptions +  ') \n '
    rpopt  = ' opt=(' + oneoptions +  ') \n '
    vwopt  = ' opt=(internal,calcall) \n '
    allint = 'int=ultrafine nosym scf=xqc'
    ircfor = ' irc(forward,calcall,stepsize=3,maxpoints=10)\n int=ultrafine nosym iop(7/33=1)\n'
    ircrev = ' irc(reverse,calcall,stepsize=3,maxpoints=10)\n int=ultrafine nosym iop(7/33=1)\n'
    ircend = ' HRcc 1 1'
    ts    = ''
    wellr = ''
    wellp = ''
    if tss[0] != 'false':
        ts = tss[0]
        if tss[1] != 'false':
            wellr = tss[1]
        if tss[2] != 'false':
            wellp = tss[1]

    for meth in meths:
        if 'molpro' in meth[1]:
            prefix = meth[0]
            if 'hind' in prefix:
                prefix = 'onedtau'
            if 'hlevel' in prefix:
                prefix = 'hl'
            if 'symm' in prefix:
                prefix = 'symm'
              
            theory += meth[0] + ' ' + meth[1] + '\n\n'
            molpro  = ''
            if '.txt' in meth[2]:
                if not os.path.sep in meth[2]:
                    meth[2] = '..' + os.path.sep + meth[2]
                if io.check_file(meth[2]):
                    molpro = io.read_file(meth[2])
                else: 
                    log.info('molpro template at {} not found'.format(meth[2]))
            else:
                if meth[0] == 'level1':
                    freqcalc = True
                else:
                    freqcalc = False
                if meth[0] == 'hlevel':
                    opt = False
                else:
                    opt = True
                molpro  = build_molpro(meth,freqcalc, opt)
            molfile = prefix +  '_molpro.dat'
            io.write_file(molpro, molfile)

        elif  'g09' in meth[1]:   
            theory += meth[0] + ' ' + meth[1] + ' restart 1\n '
            if meth[0] == 'irc':
                theory +=  meth[2] + ircfor + meth[2] + ircrev + ircend
            elif 'hlevel' in meth[0]:
                theory += meth[2] + '\n' + allint
            elif meth[0] == 'level0':
                theory += meth[2] + rzpopt + allint
            elif meth[0] == 'hind_rotor' and adiabatic.lower() == 'true':
                theory += meth[2] + rpopt + 'int=ultrafine nosym \n mhr_freqs '
            else:
                theory += meth[2] + rpopt + allint
            if meth[0] == 'level1':
                theory += ' freq'
            theory += '\n\n'
            if wellr and meth[0] == 'level1':
                theory += meth[0] + '_61 ' + meth[1] + '\n '
                theory += meth[2] + vwopt + allint
                theory += ' freq\n\n'
            if wellp and meth[0] == 'level1':
                theory += meth[0] + '_51 ' + meth[1] + '\n '
                theory += meth[2] + vwopt + allint
                theory += ' freq\n\n'
            if ts and meth[0] != 'hlevel' and meth[0] != 'irc':
                theory += meth[0] + '_ts ' + meth[1] + ' restart 1\n '
                theory += meth[2] + tsopt + allint
                if meth[0] == 'level1':
                    theory += ' freq'
                if meth[0] == 'hind_rotor':
                    if adiabatic.lower() == 'true':
                        theory += '\n ' + meth[2] + rpopt + 'int=ultrafine nosym \n mhr_freqs ' 
                    else:
                        theory += '\n ' + meth[2] + rpopt + allint
                theory += '\n\n'
        elif meth[0].lower() != 'ktp':
            log.warning(meth[0] + ' is not a recognized program.\n')
            
    theory += 'End'

    return theory 

def build_estoktp(params, jobs, nreacs, nprods, tss, xyzstart, foundlist,isTS=False):
    """
    Builds estoktp.dat
    """
    stoichs   = params[0]
    reactype  = params[1]
    coresh    = params[2]
    coresl    = params[3]
    meml      = params[4]
    memh      = params[5]
    esoptions = params[6]
    ts    = ''
    wellr = ''
    wellp = ''
    if tss[0] != 'false':
        ts = tss[0]
        if tss[1] != 'false':
            wellr = tss[1]
        if tss[2] != 'false':
            wellp = tss[1]
    fullstoich = {}
    if len(stoichs) > 2:
        stoichs = stoichs[:2]
    for stoich in stoichs:
        stoich = qc.get_atom_stoich(stoich)
        if 'C' in stoich:
           if 'C' in fullstoich:
               fullstoich['C'] += stoich['C']
           else:
               fullstoich['C']  = stoich['C']
        if 'H' in stoich:
           if 'H' in fullstoich:
               fullstoich['H'] += stoich['H']
           else:
               fullstoich['H']  = stoich['H']
        if 'O' in stoich:
           if 'O' in fullstoich:
               fullstoich['O'] += stoich['O']
           else:
               fullstoich['O']  = stoich['O']
        if 'N' in stoich:
           if 'N' in fullstoich:
               fullstoich['N'] += stoich['N']
           else:
               fullstoich['N']  = stoich['N']
        if 'S' in stoich:
           if 'S' in fullstoich:
               fullstoich['S'] += stoich['S']
           else:
               fullstoich['S']  = stoich['S']
    stoich = ''
    if 'C' in fullstoich:
       if fullstoich['C'] > 1:
           stoich += 'C{:g}'.format(fullstoich['C'])
       else:
           stoich += 'C'
    if 'H' in fullstoich:
       if fullstoich['H'] > 1:
           stoich += 'H{:g}'.format(fullstoich['H'])
       else:
           stoich += 'H'
    if 'O' in fullstoich:
       if fullstoich['O'] > 1:
           stoich += 'O{:g}'.format(fullstoich['O'])
       else:
           stoich += 'O'
    if 'N' in fullstoich:
       if fullstoich['N'] > 1:
           stoich += 'N{:g}'.format(fullstoich['N'])
       else:
           stoich += 'N'
    if 'S' in fullstoich:
       if fullstoich['S'] > 1:
           stoich += 'S{:g}'.format(fullstoich['S'])
       else:
           stoich += 'S'
    eststring = ' Stoichiometry\t' + stoich.upper()
    
    PossibleRxns = ['addition','abstraction','isomerization','betascission','well','', 'addition_well', 'isomerization_well']
    Tstype       = ['TS','wellr','wellp']
  
    if reactype.lower() not in PossibleRxns:
        log.warning('ReactionType ' + reactype +' is unrecognized, please use: Addition, Abstraction, Isomerization, or Betascission')
    elif reactype != '' and reactype.lower() != 'well':
        eststring  += '\n ReactionType   ' + reactype.split('_')[0]
        nts = 1
        if reactype.lower() == 'addition_well' or reactype.lower() == 'isomerization_well':
            nts -= 1
        eststr = ''
        if wellr:
            if wellr.lower() == 'false':
                wellr = ''
            elif wellr.lower().startswith('find'):
                nts += 1
                eststr += '\n WellR findgeom level1'
            else:
                nts += 1
                eststr += '\n WellR'
        if wellp:
            if wellp.lower() == 'false':
                wellp = ''
            elif wellp.lower().startswith('find'):
                nts += 1
                eststr += '\n WellP findgeom level1'
            else:
                nts += 1
                eststr += '\n WellP'
        eststring += '  ' + str(nts) + 'TS' + eststr
    if 'Irc' in jobs:
        eststring += '\n Variational'
    for line in esoptions.split(','):
        eststring += '\n {}'.format(line.strip())
    if nprods > 0:
        eststring += '\n Prod'
        if nprods > 1:
            eststring += 's'
    eststring +='\n Debug  2\n'
    for job in jobs:
        if 'kTP' in job or 'irc' in job.lower():
            eststring += ' ' + job
        else:
            for n in range(nreacs):
                if 'Opt_1' in job:
                    eststring += '\n ' + job.rstrip('_1') + '_Reac' + str(n+1) + '_1'
                elif 'Opt' in job and xyzstart == '0' and foundlist[n] and not isTS:
                        eststring += '\nn' + job + '_Reac' + str(n+1)
                elif not isTS:
                    eststring += '\n ' + job + '_Reac' + str(n+1)
 
            for n in range(nprods):
                if 'Opt_1' in job:
                    eststring += '\n ' + job.rstrip('_1') + '_Prod' + str(n+1) + '_1'
                elif 'Opt' in job and xyzstart == '0' and foundlist[n+nreacs] and not isTS:
                        eststring += '\nn' + job + '_Prod' + str(n+1)
                elif not isTS:
                    eststring += '\n ' +job + '_Prod' + str(n+1)
 
            for n in range(3):
                if tss[n] and tss[n] != 'false':
                    if 'Opt_1' in job:
                        eststring += '\n ' + job.rstrip('_1') + '_' + Tstype[n]  + '_1'
                    elif 'Tau' in job and isTS:
                      #  if n < 1:
                            eststring += '\n ' +job + '_' + Tstype[n] 
                    elif 'Opt' in job and n < 1 and isTS:
                        eststring += '\n Grid_' + job + '_' + Tstype[n]
                        eststring += '\n ' + job + '_' + Tstype[n] + '_0'
                        if 'nOpt' in job:
                            eststring += '\n nTauo_' + Tstype[n]
                        else:
                            eststring += '\n Tauo_' + Tstype[n]
                    elif 'Opt' in job and (wellp or wellr) and isTS:
                        if n == 1 and wellr.lower() != 'find' and wellr.lower() != 'findgeom':
                            eststring += '\n ' + job + '_WellR'
                        elif n == 2 and  wellp.lower() != 'find' and wellp.lower() != 'findgeom':
                            eststring += '\n ' + job + '_WellP'
                    elif 'Opt' not in job:
                        eststring += '\n ' + job + '_' + Tstype[n]
                    
        eststring += '\n'
 
    eststring += '\nEnd'
    eststring += '\n ' + coresl + ',' + coresh + '\n numprocll,numprochl\n'
    eststring += ' ' + meml + 'MW  ' + memh + 'MW\n gmemll gmemhl\n'
 
    return eststring

def is_dummy(atom):
   if 'xe' in atom.lower():
      return False
   elif 'x' in atom.lower():
      return True
   return False

def count_dummy(atoms):
  dummy = 0
  for atom in atoms:
      if is_dummy(atom[0]):
          dummy += 1
  return dummy

def update_interns(n, atoms, measure ,angles):
    """
    Converts internal coordinate information from x2z to the form for
    a Z-Matrix required to run EStokTP
    """
    if len(atoms) == 0:
        return atoms, measure, angles
    for index,atom in enumerate(atoms):
        if n == 2:
            atoms[index][0] = atom[0].lower() + str(index+99)
        else:
            atoms[index][0] = atom[0].lower() + str( index+1)
        if len(atom) > 1:
            atoms[index][1] = atoms[int(atom[1])-1][0]
            if n > 1:
                atoms[index][2] = atoms[index][2].upper().replace('R','R10')
            if len(atom) > 3:
                atoms[index][3] = atoms[int(atom[3])-1][0]
                if n > 1:
                    atoms[index][4] = atoms[index][4].upper().replace('A','A10')
                if len(atom) > 5:
                    atoms[index][5] = atoms[int(atom[5])-1][0]
                    if n > 1:
                        atoms[index][6] = atoms[index][6].upper().replace('D','D10')
    if n == 2:
        for i in range(len(measure)):
            measure[i][0] = measure[i][0].upper().replace('R','R10')
            measure[i][0] = measure[i][0].upper().replace('A','A10')
            measure[i][0] = measure[i][0].upper().replace('D','D10')
        for i in range(len(angles)):
            angles[i] = angles[i].upper().replace('D','D10')
    return atoms, measure, angles

def find_period(i, zmat,hin, ijk):
    """
    Rough way of determining internal symmetry number (Hydrogen counting)
    """
    period = 0
    if ijk and i == 0:
        hinatom = zmat[int(ijk[0])-1][0] 
        for row in zmat:
            if len(row) > 6:
                if row[6] == hin:
                     if row[0] == hinatom:
                         period = 1
                     elif row[5] == hinatom:
                         period = 1
    if ijk and i == 1:
         hinatom = zmat[0][0]
         for row in zmat:
            if len(row) > 6:
                if row[6] == hin:
                     if row[1] == hinatom:
                         period = 1
                     elif row[3] == hinatom:
                         period = 1
    if not period:
        sym1 = sym2 = k = 0
        for i, row in enumerate(zmat[3:]):
            if hin.upper() in row[6]:
                k = i+3
        atom1 = zmat[k][1]
        atom2 = zmat[k][3]
        for row in zmat[1:]:
            if atom1 == row[1]:
                if 'h' in row[0]:
                    sym1 += 1
                elif atom2 != row[0]:
                    sym1 -= 100
            elif atom1 == row[0]:
                if 'h' in row[1]:
                    sym1 += 1
                elif atom2 != row[1]:
                    sym1 -= 100
            elif atom2 == row[1]:
                if 'h' in row[0]:
                    sym2 += 1
                elif atom1 != row[0]:
                    sym2 -= 100
            elif atom2 == row[0]:
                if 'h' in row[1]:
                    sym2 += 1
                elif atom1 != row[1]:
                    sym2 -= 100
        if sym1 <= 1 or sym2 <=1:
            period =  max(sym1, sym2)
        elif sym1 == sym2:
            period =  sym1 #+ sym2
        else:
            period = sym1 * sym2
        if period < 1:
            period = 1
    return period

def tau_hind_str(atomslist, jobs, angleslist, interval, nsteps, mdtau, ijk =[], babs=0):
    anglen = 0
    if angleslist:
        anglen = len(angleslist[0])
        if len(angleslist) > 1:
            anglen += len(angleslist[1])
    allangles = []
    periods   = []
    #TAU
    string  = '\nntau number of sampled coordinates\n'
    #if babs > 1:
    #    anglen += 2 
    #elif babs > 0:
    if babs > 0:
        anglen += 1
          
    string += str(anglen) + '\n'
    string += ' -->nametau, taumin, taumax\n'
    for i, angles in enumerate(angleslist):
        for angle in angles:
            periodicity = find_period(i, atomslist[i], angle, ijk)
            string += angle + ' 0 ' + str(interval) + '\n'
            allangles.append(angle)
            periods.append(periodicity)
    if babs > 0: 
        string += 'BABS2 0 ' + str(interval) +  '\n'  
        #if babs > 1:
        #    string += 'BABS3 0 ' + str(interval) +  '\n'  

    #1 and 2D HIND
    string += '\nnhind\n'
    if '1dTau' in jobs:
        string += str(anglen) + '\n'
        string += ' -->namehind,hindmin,hindmax,nhindsteps,period\n'
        for i, hin in enumerate(allangles):
            periodicity = periods[i]
            string += hin + ' 0 ' + str(float(interval)/periodicity)  + ' ' + str(int(round(float(nsteps)/periodicity))) + ' ' + str(periodicity)  + '\n'  
        if babs > 0: 
            string += 'BABS2 0 ' + str(float((interval)))   + ' ' + str(int(round(float(nsteps)))) + ' 1\n' 
            #if babs > 1: 
            #    string += 'BABS3 0 ' + str(float((interval)))   + ' ' + str(int(round(float(nsteps)))) + ' 1\n' 
    else:
        string += '0\n'
    if mdtau and anglen >= int(mdtau):
        mdtau   = mdtau.strip('D').strip('d')
        string += '\nnhind' + mdtau + 'D\n'
        string += '1\n'
        string += ' -->namehind,hindmin,hindmax,nhindsteps,period\n'
        #if babs > 1 and anglen == int(mdtau):
        #    string += 'BABS2 0 ' + str(float((interval)))   + ' ' + str(int(round(float(nsteps)))) + ' 1\n'
        #    string += 'BABS3 0 ' + str(float((interval)))   + ' ' + str(int(round(float(nsteps)))) + ' 1\n'
        #    mdtau = str( int(mdtau) - 2 )
        #elif babs > 0 and anglen == int(mdtau):
        if babs > 0 and anglen == int(mdtau):
            string += 'BABS2 0 ' + str(float((interval)))   + ' ' + str(int(round(float(nsteps)))) + ' 1\n'
            mdtau = str( int(mdtau) - 1 )
        if mdtau >= anglen:
            for i in range(int(mdtau)):
                periodicity = periods[i]
                string += allangles[i] + ' 0 ' + str(float(interval)/periodicity) + ' ' + str(int(round(float(nsteps)/periodicity))) + ' ' + str(periodicity)  + '\n'   
    #    if isTS: 
    #        string += 'BABS2 0 ' + str(interval) +  '\n'  
    return string, anglen

def build_obzmat(smiles):
    """
    Uses QTC interface by Murat Keceli to Openbabel to generate zmat file based on smile string
    
    INPUT:
    smiles -- SMILES string for molecule
    OUTPUT:
    atoms   -- The atoms and their connection to previous atoms (in EStokTP format)
    measure -- The distance and angle parameters for the zmat   (in EStokTP format)
    """
    
    mol   = ob.get_mol(smiles,make3D=True)
    zmat  = ob.get_zmat(  mol)
    atoms, measure = zmat.split('\nVariables:\n')
    atoms = atoms.split('\n')
    for i in range(len(atoms)):
        atoms[i] = atoms[i].split()
    measure = measure.split('\n')
    del measure[-1]
    for i in range(len(measure)):
        measure[i] = measure[i].upper().split('= ')
    return atoms, measure

def build_obcart(smiles, mult):
    """
    Uses QTC interface by Murat Keceli to Openbabel to generate cartesian coordinates based 
    on SMILES string

    INPUT:
    smiles -- SMILES string for molecule
    mult   -- mult for molecule
    OUTPUT:
    lines  -- cartesian coordinates, formatted for a xyz file
    """
    filename    =  ob.get_smiles_filename(smiles) + '_m' + str(mult) + '.xyz'
    mol         =  ob.get_mol(smiles,make3D=True)
    lines       =  ob.get_xyz(mol)#.split('\n')
    return lines
 
def build_optout(xyzstart, XYZ, measure, angles, deletedangles, smiles, smilesfilename):
    """
    Creates the output/reac1_opt.out input for when an initial xyz is specified to start 
    EStokTP with at either level0 or level1
    """
    optim = ''
    if '0' in xyzstart:
        XYZ = XYZ.replace('g09','gaussian')
        optim = 'opt geom           1'
        for meas in measure:
            optim += '\n\t' + meas[1]              
        for i in range(len(angles)):
            optim += '\n\t{}'.format(deletedangles[i])             

        if '.log' in XYZ:
            E= str(pa.gaussian_energy(io.read_file('../' + XYZ)))
        elif '.xyz' in XYZ:
            E = io.read_file(XYZ).split('\n')[1]
        elif io.check_file('../' + smilesfilename + '.ene'):
            E = io.read_file('../' + smilesfilename + '.ene')
        elif len(XYZ.split('/')) > 2 :
            split = XYZ>split('/')
            E = io.db_get_sp_prop(smilesfilename, 'ene', None, split[0], split[1], split[2],
                                                         split[0], split[1], split[2])
            if E == None:
                log.warning('No energy found at ' + io.db_sp_path(split[0], split[1], split[2], 
                        None, smiles, split[0], split[1], split[2]) + smilesfilename + '.ene')
        else:
            E = '0'

        if E == None:
           if io.check_file(XYZ.replace('xyz','ene')):
               E = io.read_file(XYZ.replace('xyz','ene'))
        optim += '\n\t' + E 

    elif '1' in xyzstart:
        optim = 'opt level1 0'
        for meas in measure:
            optim += '\n\t' + meas[1]                
        for i in range(len(angles)):
            optim += '\n\t{}'.format(deletedangles[i])             
 
    return optim

def read_cart(smiles, mult, reactype = ''):
    """
    Finds and Reads an initial XYZ file (and set i,j,k sites for a transition state)
    """

    if 'geomdir' in smiles:
        smilesfilename = 'geomdir/' +  ob.get_smiles_filename(smiles.split('/')[1])
    else:
        smilesfilename = ob.get_smiles_filename(smiles)
    cartlines = ''
    ijk = [0, 0, 0, 0]
    if reactype in ['wellp', 'wellr', 'ts']:
        if io.check_file('../' + smilesfilename + '.xyz'):               
            cartlines = io.read_file('../' + reactype + '.xyz')
    if io.check_file('../' + smilesfilename + '_m' + str(mult) + '.xyz'):               
      #  cartlines = io.read_file('../' + smilesfilename + '_m' + str(mult) + '.xyz').split('\n\n')[1]
      #  cartlines =  str(len(cartlines.split('\n'))-1) + ' \n\n' + cartlines
        cartlines = io.read_file('../' + smilesfilename + '_m' + str(mult) + '.xyz')
    elif io.check_file('../' + smilesfilename + '.xyz'):                
        #cartlines = io.read_file('../' + smilesfilename + '.xyz').split('\n\n')[1]
        #cartlines =  str(len(cartlines.split('\n'))-1) + ' \n\n' + cartlines
        cartlines = io.read_file('../' + smilesfilename +  '.xyz')
    elif io.check_file('../' + smilesfilename + '_m' + str(mult) + '.geo'):
        cartlines = io.read_file('../' + smilesfilename + '_m' + str(mult) + '.geo')
        cartlines = str(len(cartlines.split('\n'))-1) + ' \n\n' + cartlines
    elif io.check_file('../' + smilesfilename + '.geo'):
        cartlines = io.read_file('../' + smilesfilename + '.geo')
        cartlines = str(len(cartlines.split('\n'))-1) + ' \n\n' + cartlines
    else:
        try:
            cartlines = build_obcart(smiles, mult)
            log.warning('No .geo or .xyz provided\n...Using openbabel instead')
        except:
            cartlines = ''
        return cartlines, ijk, False

    cartlines = cartlines.splitlines()
    #Reorder cartesian coordinates for zmat if isomerization
    if 'isomerization' in reactype.lower():
        toplines = [0,0,0,0]
        botlines = []
        for i,line in enumerate(cartlines[2:], start=1):
            if len(line.split()) > 4:
                toplines[int(line.split()[0])-2] = line
                #toplines[int(line.split()[0])-1] = line
            else:
                botlines.append(line)
        if toplines[0]:
            cartlines = cartlines[:2]
            cartlines.extend(toplines)
            #cartlines.extend(toplines[::-1])
            cartlines.extend(botlines)
        #order  = [1, 2, 3, 4]
        #k = 0
        #while (order != sorted(order, reverse = True) and k < 100):
        #    j = 5
        #    m = 0
        #    k += 1
        #    newcart = cartlines[:2]
        #    for i,line in enumerate(cartlines[2:], start=1):
        #        if len(line.split()) > 4:
        #            order[int(line.split()[0])-1] = i
        #            if j < i and m < line.split()[0]:
        #                index = len(newcart)-1
        #                newcart.insert(index, line)
        #                j = i - 1 
        #            else:
        #                newcart.append(line) 
        #                j = i
        #            m = line.split()[0]
        #        else:
        #             newcart.append(line) 
        #    cartlines = newcart
    #Find if i,j,k site is specified:
    for i,line in enumerate(cartlines[1:], start=1):
        if len(line.split()) > 4:
            cartlines[i] = '   '.join(line.split()[1:])
            if 'isomerization' in reactype.lower():
                line = line.split()
                line[0] = line[0].replace('1','5').replace('3','1').replace('4','3')
                line = '  '.join(line)
            if line.split()[0] == '1':
                ijk[1] = str(i)
            elif line.split()[0] == '2':
                ijk[0] = str(i)
            elif line.split()[0] == '3':
                ijk[2] = str(i)
            elif line.split()[0] == '5':
                ijk[3] = str(i)
            elif line.split()[0] == '4':
                temp = cartlines[i]
                del cartlines[i]
                cartlines.insert(2,temp)
    cartlines = '\n'.join(cartlines)
    return cartlines, ijk, True

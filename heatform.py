#!usr/bin/python

import re
import numpy as np
import os

def get_atomlist(mol):
    """
    Makes a list of all atoms in a molecule
    INPUT:
    mol      - stoichiometry of molecule
    OUTPUT:
    atomlist - list of distinct atoms in that molecule
    """
    atomlist = []
    elements = {'He','Li','Be','Ne','Na','Mg','Al','Si','Cl','Ar',
      'Ca','Sc','Ti','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge', 
      'As','Se','Br','Kr','C','B','H','O','F','S','N','P','K','V'}
    for el in elements:
        if el in mol:
           atomlist.append(el)
           mol = mol.replace(el,'')
    return atomlist

def select_basis(atomlist):
    """
    Given a list of atoms, generates a list of molecules
    that is best suited to serve as a basis for those atoms
    INPUT:
    atomlist - list of atoms
    OUPUT:
    basis    - recommended basis as a list of stoichiometries
    """
    count = len(atomlist)-1
    basis = []
    i= 0
    if 'N' in atomlist and i <= count:
        basis.append('NH3')
        i += 1
    if 'S' in atomlist and i<= count:
        basis.append('SO2') 
        i += 1
    if 'C' in atomlist and i<= count:
        basis.append('CH4')
        i += 1
    if 'H' in atomlist and i<= count:
        basis.append('H2')
        i += 1
    if 'H' in atomlist and  i<= count:
        basis.append('H2O')
        i += 1
    if 'C' in atomlist and i<= count:
        basis.append('CO2')
        i += 1
    if 'O' in atomlist and i<= count:
        basis.append('O2')
        i += 1
    if 'C' in atomlist and i<= count:
        basis.append('H2CO')
        i += 1
    if 'C' in atomlist and i<= count:
        basis.append('CH3CH3')
        i += 1
    if 'O' in atomlist and i<= count:
        basis.append('CH3OH')
        i += 1
    return basis


def get_stoich(mol,atomlist):
    """
    Given a molecule's stoichiometry and a list of atoms, finds the
    number of each atom that the molecule contains
    INPUT:
    mol       - molecule stoichiometry
    atomlist  - list of atoms 
    OUTPUT:
    stoich    - list of numbers corresponding to the number of each atom 
                in the atomlist that the molecule contains
    """
    stoichlist = np.zeros(len(atomlist))

    for i, atom in enumerate(atomlist):
        val = 0
        if atom in mol:
            a = re.compile(atom + '(\d*)')
            a = a.findall(mol)
            for b in a:
                if b == '':
                    b = '1'
                val += float(b)
        stoichlist[i] = val

    return stoichlist


def form_mat(basis,atomlist):
    """
    Form a matrix for a given basis and atomlist
    INPUT:  
    basis     - basis of molecules
    atomlist  - list of atoms (all atoms that appear 
                in basis should be in atomlist)
    OUTPUT:
    mat       - matrix (length of basis by length of atomlist)
                (square if done right)
    """
    mat = np.zeros((len(atomlist),len(atomlist)))
    for i,mol in enumerate(basis):
        mat[i] = get_stoich(mol,atomlist)
    mat = mat.T
     
    return mat

def comp_coeff(mat,stoich):
    """
    Finds the coefficients that solve C = M^-1 S.  For our purposes C are the coefficients [a,b,c,d..] for 
    the basis [R1, R2, R3, R4...] that give CxHyOzNw =  aR1 + bR2 + cR3 + dR4 where [x,y,z,w...] is S, our 
    stoichiometry vector for a molecule.  M is a nonsingular matrix that puts the basis in terms of 
    a list of atoms
    INPUT:
    mat    - nonsingular matrix 
    stoich - list of numbers that give a molecules stoichiometry in terms of a list of atoms
    OUTPUT:
    coeff  - coefficients [a,b,c,d...] as described above
    """
    mati = np.linalg.inv(mat)
    coeff = np.dot(mati,stoich)

    return coeff
    
def update_dic(mol,dic):
    """
    Finds if we already have the energy for a molecule in our dictionary, and if not checks if we are in 
    its EStokTP directory that we can get that information from and adds it to the dictionary
    INPUT:
    mol   - molecule
    dic   - molecule/energy dictionary
    OUTPUT:
    dic   - updated molecule/energy dictionary
    """
    if mol in dic:
        return dic
    elif os.path.exists('geoms/reac1_l1.xyz'):
        lines = open('geoms/reac1_l1.xyz','r').read()
        E = float(lines.split('\n')[1])
        dic[mol] = E
        return dic
    
def E_dic(mol,dic):
    """
    More-or-less pointless functions that just returns a dictionary value but may be extended for some idiotproofing later
    """
    if mol in dic:
        return dic[mol]
    return 0

def comp_energy(mol,basis,coefflist,E):
    """
    Uses the coefficients [a,b,c...] obtained from C = M^-1 S to find 
    delH(CxHyOz) = adelH(R1) + bdelH(R2) + cdelH(R3) + Eo(CxHyOz) - aEo(R1) - bEo(R2) -cEo(R3)
    where Rn are our basis molecules, delH(Rn) are their heats of formation, and Eo(Rn) are their
    electronic energies computed at the same level of theory as Eo(CxHyOz)
    INPUTS:
    mol       - molecule named stoichiometrically
    basis     - selected basis molecule list
    coefflist - coefficients [a,b,c,d...] described above
    E         - electronic energy of molecule
    OUTPUTS:
    lE        - 0K heat of formation of molecule
    hE        - 298K heat of formation of molecule
    sigma     - uncertainty estimate (i.e., delH = lE +/- sigma)
   
    """
    zeroTheat = {'CH4':-66.550,'CH3CH3':-68.29,'CH3OH':-189.83,'NH3':-38.562,'H2O':-286.300,'CO2':-393.109,'H2':0,'H2CO':-105.349,'O2':0,'SO2':-296.810}
    highTheat = {'CH4':-74.520,'CH3CH3':-83.91,'CH3OH':-200.71,'NH3':-45.554,'H2O':-285.828,'CO2':-393.475,'H2':0,'H2CO':-109.188,'O2':0}
    uncert    = {'CH4':.057,'CH3CH3':.14,'CH3OH':.16,'NH3':0.03,'H2O':0.027,'CO2':.015,'H2':0,'H2CO':.099,'O2':0,'SO2':.2}
    elecenerg = {'CH4':-40.52614,'CH3CH3':-79.84164,'CH3OH':-115.73487,'H2CO': -114.51152,'CO2':-188.59039,'H2O':-76.434049,'H2':-1.178539,'O2':-150.26605}
    elecuncert= (-40.52614 - -40.496760)/(2*40.526) + (-114.511 - -114.487) / (2*114.511)

    if E == -1000.0:
        update_dic(mol,elecenerg)
    else:
        elecenerg[mol] = E

    lE        = E_dic(mol,elecenerg)
    hE        = E_dic(mol,elecenerg)
    var       = 0

    for i,bas in enumerate(basis):

        lE  +=  coefflist[i] * E_dic(bas,zeroTheat) * 0.00038088
        hE  +=  coefflist[i] * E_dic(bas,highTheat) * 0.00038088

        var += (coefflist[i] * E_dic(bas,uncert) * 0.00038088   )**2
 
        elecenerg  =  update_dic(bas,elecenerg)
        lE  -=  coefflist[i] * E_dic(bas,elecenerg) 
        hE  -=  coefflist[i] * E_dic(bas,elecenerg) 
                                                                
        var += (coefflist[i] * elecuncert * E_dic(bas,elecenerg) )**2

    sigma = np.sqrt(var)
    return lE, hE, sigma
    
def check(clist, basis,stoich):
    """
    Makes sure nothing funky happened while computing coefficients
    """
    check = np.zeros(len(clist))
    statement = 'Coefficients produce correct stoichiometry'
    for i, c in enumerate(clist):
       check += c * get_stoich(basis[i],atomlist)
    for i, sto in enumerate(stoich):
        if not check[i] == sto:
            statement = 'Coefficients do NOT produce correct stoichiometry'
            break
    return statement

if __name__ == '__main__': 
    """
    Run heat of formation code
    """
    #SET PARAMETERS############
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                            description="SARAH!!! you haven't done this yet!!!")

    parser.add_argument('-s','--stoichiometry',type=str,default='CH3CH2CH3')
    parser.add_argument('-e','--electronic_energy',type=float,default=-1000.0)
    parser.add_argument('-b','--select_basis',type=str,default='XXX')

    args  = parser.parse_args()
    mol   = args.stoichiometry
    E     = args.electronic_energy
    basis = args.select_basis.split()
    ###########################

    #AUTO SET NONUSERDEFINED PARAMETRS## 
    basprint = 'manually select basis'
    
    atomlist = get_atomlist(mol)
    if basis[0] == 'XXX':
        basis = select_basis(atomlist)
        basprint = 'automatically generate basis'
    elif basis[0] == 'basis.dat':
        basis = io.read_file('basis.dat').split()
        basprint = 'read basis from basis.dat'
    lines =  ('\n-------------------------------------------------\n\n' +
              'HEAT OF FORMATION FOR: ' + mol +
              '\n\n-------------------------------------------------\n\nYou have chosen to ' + 
              basprint + '\n\nBasis is: ' + ', '.join(basis))
    print lines  
    for bas in basis:
         atomlist.extend(get_atomlist(bas))
    ####################################

    #COMPUTE Atomlist, stoichlist, matrix, and coefficients
    atomlist = list(set(atomlist))
    stoich = get_stoich(mol,atomlist)
    mat = form_mat(basis,atomlist)

    if np.linalg.det(mat) == 0:
        print 'Matrix is singular -- select new basis'

    clist =  comp_coeff(mat,stoich)
    ######################################################
     
    ###PRINT STUFF OUT
    lines = '\n  ' + mol + '\t' +  '\t'.join(basis) 
    for i in range(len(mat)):
       lines += '\n' + atomlist[i] + '  '
       lines += str(stoich[i]) + '    \t'
       for el in mat[i]:
           lines += str(el) + '\t'
    lines +=  '\n\nCoefficients are: '
    for co in clist:  lines += str(co) + ' '
    print lines + '\n'
    print check(clist, basis,stoich)
    ##################

    #COMPUTE AND PRINT delH###
    E =  comp_energy(mol,basis,clist,E)
    lines =  '\n        delHf(0K) \t delHf(298K) \t Uncert'
    lines += '\nA.U. \t'
    for e in E:  lines += str(e) + '\t'
    lines += '\nkJ   \t'
    for e in E:  lines += str(e/ .00038088) + '\t'
    lines += '\n\n-------------------------------------------------\n\n'
    print lines
    return E[0]
    ##########################

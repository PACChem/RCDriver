import numpy as np

def sites(lines):
    
    atoms  = lines.split('\n')
    del atoms[-1]
    natoms = len(atoms)
    mat    = np.zeros((natoms,natoms))
    
    reaclist = np.zeros(natoms)
    atomlist = []
    for i, atom in enumerate(atoms):
        if atom[3] != ' ':
            reaclist[i] = int(atom[3])
        atomlist.append(atom[5])
        elems = atom.split('{')
        del elems[0]
        for elem in elems:
             ind,val = elem.strip().rstrip('}').split(',')
             if val.upper() == 'S':
                 val = 1
             elif val.upper() == 'D':
                 val = 2
             else: val = 3
             mat[i,int(ind)-1] = val
    m=[0]
    a = True
    bond = 1
    light = ['H']
    heavy = ['C','O','N','X']
    while a == True:
        order = 0
        hcount = 0
        M = m[-1]
        if atomlist[M] != 'H':
            for k, j in enumerate(mat[M]):
                if j != 0:
                    order += j
                    if atomlist[int(k)] not in heavy:
                        hcount += 1
                        if reaclist[int(k)] == 2:
                            print '----------------\nISITE: Hydrogen abstracted from site ' + str(M+1)
                            print 'JSITE: Hydrogen is site ' + str(k+1)
                            for l,L in enumerate(mat[M]):
                                if L == 1 and l != k:
                                    print 'KSITE: up the chain is site ' + str(l+1) + '\n---------------'
                                    return str(k+1), str(M+1), str(l+1) 
                    elif reaclist[int(k)] == 3:
                        print '----------------\nISITE: Addition at site ' + str(k+1)
                        print 'JSITE: Up the chain is ' + str(M+1)
                        for l,L in enumerate(mat[M]):
                            if L == 1 and l != k:
                                print 'KSITE: further up the chain is site ' + str(l+1) + '\n---------------'
                                return str(k+1), str(M+1), str(l+1)
                    elif int(k) in m:
                        a = False
                    else:
                        bond = j
                        m.append(int(k))
    return

import os
import sys
import numpy as np

sys.path.insert(0, '/home/elliott/Packages/QTC/')
import iotools as io


def gauss_xmat(filename,natoms):
    """
    Retrieves the anharmonic constant matrix from Gaussian logfile 
    INPUTS:
    filename - name of gaussian logfile
    natoms   - number of atoms in molecule
    OUTPUT:
    xmat     - anharmonic constant matrix (nmode by nmode)
    """ 
    full = io.read_file(filename)
    nmodes = 3*natoms-6 
    lines = full.split('X matrix')[1].split('Resonance')[0]
    lines = lines.split('\n')
    del lines[0]
    del lines[-1]
    
    xmat = np.zeros((nmodes, nmodes))
    
    rangemod = 1
    if nmodes%5 == 0:
       rangemod = 0
    marker = 0

    for m in range(0,nmodes/5+rangemod):
        length = nmodes - m * 5 
        a= np.array( lines[marker+1:marker+length+1])
        for i in range(length):
            for j in range(0,len(a[i].split())-1):
                xmat[m*5 + i,m*5 + j] = a[i].split()[j+1]
                xmat[m*5 + j,m*5 + i] = a[i].split()[j+1]
        marker += length+1

    return xmat


def get_freqs(filename):
    """
    Pulls the frequencies out from EStokTP me output file 
    INPUT:
    filename - name of EStokTP output file (reac1_fr.me or reac1_unpfr.me)
    OUTPUT:
    freqs    - frequencies obtained from output file
    order    - in case the frequencies were reordered when sorting, keeps 
               track of which index of freqs corresponds to which normal mode
    """ 
    full = io.read_file(filename)
    full = full.strip('\n')
    full = full.split('[1/cm]')[1].split('Zero')[0] 
    full = full.split()
    nfreqs = full[0]
    freqs = full[1:]
    #[freq=float(freq) for freq in freqs]
    freqs = np.array(map(float, freqs))
    a= freqs.argsort()[::-1]
    freqs = np.sort(freqs)[::-1]
    return freqs.tolist(), a.tolist()


def find_hinfreqs(proj,unproj,order):
    """
    Compares the frequencies from EStokTP projected and unprojected frequency
    output to determine which normal modes are hindered rotors
    INPUTS:
    proj   -  frequencies after projection
    unproj -  unprojected frequencies
    order  -  in case the frequencies were reordered when sorting, keeps track of 
              which index of unproj corresponds to which normal mode
    """
    for i in range(len(proj)):
        length = len(unproj)-1
        k = 0
        closeenough = 5
        if abs(proj[0]-unproj[k]) < closeenough:
            del proj[0]
            del unproj[k]
            del order[k]
        elif k < length:
            k+=1
            if abs(proj[0]-unproj[k]) < closeenough:
                del proj[0]
                del unproj[k]
                del order[k]
            elif k < length:
                k+=1
                if abs(proj[0]-unproj[k]) < closeenough:
                    del proj[0]
                    del unproj[k]
                    del order[k]
        modes = [mode+1 for mode in order]
    return modes

def remove_modes(xmat,modes):
    """
    Removes specified modes from anharmonic constant matrix
    INPUTS:
    xmat  - anharmonic constant matrix
    modes - the modes to delete from the matrix (with 1 being the first mode)
    OUTPUTS:
    xmat  - anharmonic constant matrix with columns and rows deleted for specified modes
    """ 
    modes.sort(reverse=True)
    modeindex = [mode-1 for mode in modes]
    
    for index in modeindex:
        xmat = np.delete(xmat,index,0)
        xmat = np.delete(xmat,index,1)
    return xmat

def gauss_anharm_inp(filename):
    """
    Forms the Gaussian input file for anharmonic frequency computation following an EStokTP 
    level 1 computation on a molecule
    INPUT:
    filename - EStokTP output file to read (reac1_l1.log)
    OUTPUT:
    zmat     - lines for entire guassian input file (not just the zmat part, its poorly named)
    """
    full = io.read_file(filename)
    full = full.split('Z-matrix:')
    zmat = full[0].split('***************************')[2].replace('*','')
    zmat = zmat.split('Will')[0]
    zmat = ' ' + zmat.lstrip() 
    zmat += full[0].split('-------------------------------------------')[3].replace('--','')
    zmat += '# scf=verytight nosym Freq=Anharmonic Freq=Vibrot\n'
    zmat += '\nAnharmonic computation\n'
    zmat += full[1].split('       Variables:')[0]
    zmat += 'Variables:\n'
    zmat = zmat.replace('Charge = ','')
    zmat = zmat.replace('Multiplicity =','')
    varis = full[1].split('Optimized Parameters')[1].split('--------------------------------------')[1]
    varis = varis.split('\n')
    del varis[0]
    del varis[-1]
    for var in varis:
        var = var.split()
        zmat += ' '+  var[1] + '\t' + var[2] + '\n'
    return zmat

def write_anharm_inp(readfile='reac1_l1.log',writefile='anharm.inp'):
    
    """
    Writes Guassian input to a file given an EStokTP G09 output file name
    INPUT:
    readfile  - EStokTP output file to read (reac1_l1.log)
    writefile - name of Gaussian input file to write
    """
    zmat = gauss_anharm_inp(readfile)
    io.write_file(zmat,writefile)
    return

def run_gauss(filename, node='b456'):
    """
    Executes Guassian 
    INPUT:
    filename - name of Guassian input file
    node     - node to run it on
    """
    if io.check_file(filename):
        execute = 'cd `pwd`; export PATH=$PATH:~/bin; soft add +gcc-5.3; soft add +g09; g09 ' + filename + ' &'
        ssh ='/usr/bin/ssh'
        host =node
        os.system('exec ' + ssh + ' -n ' + host +' \"' + execute + '\"')
    
    return

def anharm_freq(freqs,xmat):
    """
    Uses anharmonic frequency matrix and harmonic frequencies to compute VPT2 anharmonic frequencies
    INPUT:
    freqs   - harmonic frequencies
    xmat    - anharmonic constant matrix
    OUTPUT:
    anharms - VPT2 anharmonic frequencies
    """
    anharms = np.zeros(len(freqs))

    for i, freq in enumerate(freqs):
        anharms[i]  = freq
        anharms[i] += 2. * xmat[i,i]

        tmp = 0
        for j in range(len(freqs)):
            if j != i:
                tmp += xmat[i,j]

        anharms[i] += 1./2 * tmp

    return anharms

if __name__ == '__main__':

    #SET PARAMETERS############
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
             description="""
       This module computes anharmonic corrections to the projected
       frequencies produced during an EStokTP 1D or MD torsional scan.  

       Requirements: 
       In addition to both the projected frequency and unprojected frequency files that 
       EStokTP puts in me_files, it requires EITHER a g09 anharmonic logfile OR for g09 
       to be available so that this module can execute a g09 anharmonic computation 
       written using an optimization logfile (usually taken from geoms/reac1_l1.log)
       """)

    parser.add_argument('-n',         '--natoms',type=int,help = 'number of atoms in the molecule. Required.',                              required=True)
    parser.add_argument('-a',      '--anharmlog',type=str,help = 'location of g09 anharmonic logfile IF unavailable, use next 3 options',   default='anharm.log')
    parser.add_argument('-l',        '--logfile',type=str,help = 'path to  optimization logfile (required if no g09 anharmfile available)', default='geoms/reac1_l1.log')
    parser.add_argument('-w',     '--writegauss',type=str,help = 'if true will write gaussian anharmonic input file',                       default='false')
    parser.add_argument('-r',       '--rungauss',type=str,help = 'if true will execute guassian anharmonic computation',                    default='false')
    parser.add_argument('-freq',    '--freqfile',type=str,help = 'path to estoktp UNprojected frequency file found in me_files',            default='me_files/reac1_fr.me')
    parser.add_argument('-unfreq','--unprojfreq',type=str,help = 'path to estoktp   projected frequency file foudn in me_files',            default='me_files/reac1_unprfr.me')
    parser.add_argument('-x',  '--computeanharm',type=str,help = 'specify false to avoid computing anharmonic correction',                  default='true')

    args      = parser.parse_args()
    eskfile   = args.logfile
    natoms    = args.natoms
    eskproj   = args.freqfile
    eskunproj = args.unprojfreq
    anharmlog = args.anharmlog
    ##########################
    if args.writegauss.lower() == 'true':
        write_anharm_inp(eskfile,'anharm.inp')
    if args.rungauss.lower() == 'true':
        run_gauss('anharm.inp')
    if args.computeanharm.lower() == 'true':
        xmat = gauss_xmat(anharmlog,natoms)
        proj, b   = get_freqs(eskproj)
        unproj, a = get_freqs(eskunproj)
        print anharm_freq(unproj,xmat)
        modes     = find_hinfreqs(proj,unproj,a)
        xmat      = remove_modes(xmat,modes)
        proj, b   = get_freqs(eskproj)
        print anharm_freq(proj,xmat)

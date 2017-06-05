
import os
import sys
sys.path.insert(0, '/home/elliott/Packages/QTC/')
import iotools as io
import re

class REACTION:
    def __init__(self,lines):
        self.label     = get_label(lines)
        self.reactype  = get_type(lines)
        self.reactants = get_reactants(lines)
        self.products  = get_products(lines)
        self.TS        = get_ts(lines)
        self.nTS       = len(self.TS)
        #self.arr       = get_Arr(lines)   

class SPECIES:
    def __init__(self,lines):
        self.label   = get_label(lines)
        self.smiles  = get_smiles(lines)
        self.mult    = get_multiplicity(lines)

class TRANS:
    def __init__(self,lines):
        self.label   = get_label(lines)
        #self.smiles  = get_smiles(lines)   coming soon to RMG
        self.mult    = get_multiplicity(lines)

def get_label(data):
    search = re.compile('label = \\\'(.*)\\\'')
    return search.findall(data)[0]

def get_type(data):
    search = re.compile('reactype = \\\'(.*)\\\'')
    return search.findall(data)[0]
     
def get_reactants(data):
    search = re.compile('reactants = \[(.*)\]')
    return search.findall(data)[0].replace("'","").replace(" ","").split(',')

def get_products(data):
    search = re.compile('products = \[(.*)\]')
    return search.findall(data)[0].replace("'","").replace(" ","").split(',')

def get_ts(data):
    search = re.compile('transitionState = \\\'(.*)\\\'')
    return search.findall(data)

def get_smiles(data):
    search = re.compile('structure = SMILES(.*),')
    search = search.findall(data)
    return search[0].replace('(','').replace(')','').replace("'",'')

def get_multiplicity(data):
    search = re.compile('spinMultiplicity = (.*),')
    return int(search.findall(data)[0])

def get_Arr(data):
    search = re.compile('Arrhenius\(A=\((.*)\\\'')
    return search.findall(data)

#if __name__ == '__main__':
#    rmgfile = 'rmg_input_example.py'
#    full    = io.read_file(rmgfile)
#    inputs  = full.split('\r\n\r\n')
#    dic   ={}
#    tsdic ={}
#    for inp in inputs:
#        if 'species' in inp:
#            Spec = rg.SPECIES(inp)
#            dic[Spec.label] = [Spec.smiles, Spec.mult]
#        if 'transitionState' in inp:
#            Trans = rg.TRANS(inp)
#            tsdic[Trans.label] = [Trans.smiles, Trans.mult]
#        if 'reaction' in inp:
#            Reac = rg.REACTION(inp)
#            self.reactype = Reac.reactype
#            reactants     = Reac.reactants
#            products      = Reac.products
#            tss           = Reac.TS
#            self.reacs = []
#            self.prods = []
#            self.TS    = []
#            self.nTS   = Reac.nTS
#            for reac in reactants:
#                if not reac in dic:
#                    print 'incomplete RMG data'
#                    break
#                else:
#                    self.reacs.append(dic[reac][0])
#            for prod in products:
#                if not prod in dic:
#                    print 'incomplete RMG data'
#                    break
#                else:
#                    self.prods.append(dic[prod][0])
#            for ts in tss:
#                if not ts in tsdic:
#                    print 'incomplete RMG data'
#                    break
#                else:
#                    self.ts.append(tsdic[ts][0])
#    return

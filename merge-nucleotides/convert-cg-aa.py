#!/usr/bin/env python

import numpy as np
import sys, argparse, os
sys.path.insert(0, os.environ["ATTRACTTOOLS"])
sys.path.insert(0, os.environ["SCRIPTS"])
from rmsdlib import multifit

def npy2to3(npy):
    if len(npy.shape) == 2:
        npy = npy.reshape(npy.shape[0], npy.shape[1]/3, 3)
    else:
        assert len(npy.shape) == 3
    return npy

def npy3to2(npy):
    if len(npy.shape) == 3:
        npy = npy.reshape(npy.shape[0], 3*npy.shape[1])
    else:
        assert len(npy.shape) == 2 and npy.shape[1]%3 == 0
    return npy

def fit_another_multi_npy(a, cg, ref):
    a = npy2to3(a)
    cg = npy2to3(cg)
    COM = cg.sum(axis=1)/cg.shape[1]
    rotation, translation, RMSD = multifit(cg, ref)
    rot = np.transpose(rotation, axes=(0,2,1))
    centered = a - COM[:,None,:]
    rotated = np.einsum('...ij,...jk->...ik',centered,rot)
    fitted = rotated + COM[:,None, :]
    translated = fitted - translation[:,None,:]
    return translated, RMSD

############
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('cg_npy')
parser.add_argument('aa_npy')
parser.add_argument('--libcg', nargs='+', help="list of cg monomer libraries")
parser.add_argument('--libaa', nargs='+', help="list of aa monomer libraries")
parser.add_argument('--noext', help="do not replace termini", action="store_true")
args = parser.parse_args()
############
chains = np.load(args.cg_npy)
libraries_aa = [ np.load(lib) for lib in args.libaa ]
libraries_cg = [ np.load(lib) for lib in args.libcg ]

Nres = len(libraries_aa)
Nchains = len(chains)
Nat_aa = sum([lib.shape[1] for lib in libraries_aa])

if len(chains.shape) == 2:
    if chains.shape[1] == 3:
        chains = chains.reshape((1,chains.shape[0],3))
        Nchains = 1
    else:
        chains = npy2to3(chains)

atoms_cg = [0] # first atom in each cg residue
for lib in libraries_cg:
    atoms_cg.append(atoms_cg[-1]+lib.shape[1])
atoms_cg.append(atoms_cg[-1]+lib.shape[-1])

atoms_aa = [0] # first atom in each cg residue
for lib in libraries_aa:
    atoms_aa.append(atoms_aa[-1]+lib.shape[1])
atoms_aa.append(atoms_aa[-1]+lib.shape[-1])


print chains.shape
monomers_chains = [ chains[:, atoms_cg[i]:atoms_cg[i+1] ,:] for i in range(Nres) ]

r = range(Nres)
if args.noext:
    r = r[1:-1]

chains_aa = np.zeros( (Nchains, Nat_aa, 3) )
print chains.shape, chains_aa.shape

#j=0 ### changed 18/04/18
for j in range(len(chains)):
    chain = chains[j]
    for i in r:
        monomer = chain[atoms_cg[i]:atoms_cg[i+1],:]
        fitted, RMSD = fit_another_multi_npy(libraries_aa[i], libraries_cg[i], monomer)
        best = fitted[ RMSD.argmin(axis=0)]
        chains_aa[j, atoms_aa[i]:atoms_aa[i+1]] = best
        #for k in range(atoms[i], atoms[i+1]):
        #    chains_aa[j, k,:] = best[k-atoms[i],:]

print chains.shape
np.save(args.aa_npy, chains_aa)

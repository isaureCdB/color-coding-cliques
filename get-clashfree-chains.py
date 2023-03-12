#!/usr/bin/env python3

import argparse, sys
import numpy as np
from math import *

# Select non-clashing chains of fragments
# chains are in format:
# #header <mean (root-mean-sq) ligand rmsd> <mean (geometric mean) rank>  <rms-overlap-rmsd> <ranks> <ligand rmsds>

#######################################
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('chains', help="chains of pose indices, in txt") #AAA-5frag-1.8A-100000.chains
parser.add_argument('clashes', help="clashing dist in A") #AAA-5frag-1.8A-clash1.5.npy
parser.add_argument('--chains', help="print clashfree chains")
parser.add_argument('--indices', help="print indices clashfree chains, from 1")

args = parser.parse_args()
#######################################

lines = [ l for l in open(args.chains).readlines()]
nfrag = int(0.5 * (len(lines[1].split()) - 3))

chains = [ [int(i)-1 for i in l.split()[3:3+nfrag]] for l in lines[1:]]
clashes = np.load(args.clashes)

sel = []
# to print not clashing chains:
for nc, chain in enumerate(chains):
    clash = False
    for i in range(nfrag):
        for j in range(i+1, nfrag):
            if [chain[i], chain[j]] in clashes:
                break
    if not clash:
        sel.append(nc)

if args.chains:
    print(lines[0], end="")
    for i in sel:
        print(lines[i+1], end="")

if args.indices:
    for i in sel:
        print(i+1)


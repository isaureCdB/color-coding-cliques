#!/usr/bin/env python3

'''
Get list of pairs of incompatible poses 
because only seen at one same position in chains
'''

import numpy as np
import argparse

#######################################
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('graph', help="graph of poses connectivity, in npz")
parser.add_argument('--npz', help="output in npz format")
parser.add_argument('--npy', help="output in npy format")
parser.add_argument('--txt', help="output in txt format")
args = parser.parse_args()
#######################################

def get_single_position_pairs(pools):
    # Get list of pairs of poses non-connectable
    # because only seen at one same position in chain
    all_pairs = set()
    nfrag = len(pools)
    #compare each pool to union of other pool
    for i in range(nfrag):
        pool = pools[i]
        for j in range(nfrag):
            if i != j:
                pool = [ i for i in pool if i not in pools[j]]
        pairs = set([ (a, b) for a in pool for b in pool if a<b ])
    all_pairs.update(pairs)
    return sorted(list(all_pairs))

def get_pools(gaph):
    nfrag = graph["nfrags"]
    # get the pool of possible poses at each position in chain
    pools = []
    for f in range(nfrag-1):
        f_pool = set(graph["interactions-%i"%f][:,0])
        pools.append(f_pool)
    f_pool = set(graph["interactions-%i"%(nfrag-2)][:,1])
    pools.append(f_pool)
    return pools

graph = np.load(args.graph)
pools = get_pools(graph)
outp = get_single_position_pairs(pools)

if args.npz is not None:
    np.savez(args.npz, outp)

if args.npy is not None:
    np.save(args.npy, outp)

if args.txt is not None:
    f = open(args.txt, "w")
    for o in outp:
        print("%i %i"%(o[0], o[1]), file=f)
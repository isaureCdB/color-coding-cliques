#!/usr/bin/env python3

'''
get list of non-clashing overlapping pairs of poses:
 - pairs overlapping by 2 nucleotides (frag1_n2n3 frag2_n1n2)
 - pairs overlapping by 1 nucl (frag1_n3 frag3_n1) and not clashing (frag1_n1 frag3_n3)
'''

import numpy as np
import argparse, sys
import random
import opt_einsum
from scipy.spatial import KDTree

#######################################
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('connection_graph', help="graph of poses connectivity, in npz")
parser.add_argument('coor', help="atomic coordinates in numpy format")
parser.add_argument('--atoms_n1', help="atom indices of nucl 1", nargs='+')
parser.add_argument('--atoms_n3', help="atom indices of nucl 3", nargs='+')
parser.add_argument('--dist', help="clashing distance in A", default=3.5)
parser.add_argument('--mem', help="mem usage (bytes) ", default=10**8)
parser.add_argument('--npz', help="list of pairs in npz format")
parser.add_argument('--npy', help="list of pairs in npy format")
parser.add_argument('--txt', help="list of pairs in txt format")
parser.add_argument('--clashes', help="list of clashes in npy format")

args = parser.parse_args()
#######################################

CLASH_THRESHOLD = float(args.dist)
MAX_CHUNK_MEM = int(args.mem)
MAX_CHUNK_1 = int(np.sqrt(MAX_CHUNK_MEM / 8 / 3))
MAX_CHUNK_2 = int(MAX_CHUNK_MEM / 8 / 3 / 9)

# debug: smaller chunk sizes
###MAX_CHUNK_1, MAX_CHUNK_2, MAX_CHUNK_3 = 50, 1000, 1000 ###
#print(MAX_CHUNK_1, MAX_CHUNK_2, MAX_CHUNK_3, file=sys.stderr)

def detect_clashes_all(coor1, coor2, pairs):
    pairs = np.array(pairs, dtype=int)
    print(pairs)
    print(len(coor1), len(coor2))
    c1 = coor1[pairs[:, 0]]
    c2 = coor2[pairs[:, 1]] 
    dif0 = c1[:, :, None, :] - c2[:, None, :, :]
    distsq = np.einsum("ijkl,ijkl->ijk", dif0, dif0)  # = (dif0*dif0).sum(axis=3)
    min_distsq = distsq.reshape(len(distsq), -1).min(axis=1)
    clashing = (min_distsq < CLASH_THRESHOLD**2)
    clashes = pairs[clashing]
    #print("CLASH!", len(clashes), len(pairs), file=sys.stderr)
    if not len(clashes):
        return None
    #print(clashes)
    return clashes

def detect_clashes_COM(coor1, coor2, pairs_to_eval, homo=False):
    #retain pairs that have close enough COMs to clash

    coor1_com = coor1.mean(axis=1)
    dif1_com = coor1 - coor1_com[:, None, :]
    dif1_com = dif1_com.reshape(-1, 3)

    coor2_com = coor2.mean(axis=1)
    dif2_com = coor2 - coor2_com[:, None, :]
    dif2_com = dif2_com.reshape(-1, 3)

    dist1_com = (dif1_com * dif1_com).sum(axis=1)
    dist2_com = (dif2_com * dif2_com).sum(axis=1)
    max_dist_com = np.sqrt(max(dist1_com.max(), dist2_com.max()))
    #print(max_dist_com, file=sys.stderr)
    com_threshold = 2 * max_dist_com + CLASH_THRESHOLD

    tree1 = KDTree(coor1_com)
    tree2 = KDTree(coor2_com)
    pairs_lists = tree1.query_ball_tree(tree2, com_threshold)
    p1len = [len(l) for l in pairs_lists]
    p1 = np.repeat(np.arange(len(coor1_com)), p1len)
    p2 = np.concatenate(pairs_lists)
    pairs = np.stack((p1, p2), axis=1)
    # debug: keep all pairs
    ###import itertools; pairs = np.array(list(itertools.product(range(len(coor1_com)), range(len(coor2_com))))) ###

    # filter pairs based on eval list
    assert pairs_to_eval[:, 0].min() >= 0, pairs_to_eval[:, 0].min()
    assert pairs_to_eval[:, 1].min() >= 0, pairs_to_eval[:, 1].min()
    _ = coor1_com[pairs_to_eval[:, 0]]
    _ = coor2_com[pairs_to_eval[:, 1]]
    pairs_key = len(coor1_com) * pairs[:, 1] + pairs[:, 0]
    pairs_key = np.sort(pairs_key)
    eval_key = len(coor1_com) * pairs_to_eval[:, 1] + pairs_to_eval[:, 0]
    eval_key = np.sort(eval_key)
    pairs_to_keep = np.intersect1d(pairs_key, eval_key) #, assume_unique=True)
    print("PAIRS BEFORE EVAL", pairs, pairs.shape, len(coor1_com))
    pairs_col1 = pairs_to_keep % len(coor1_com)
    pairs_col2 = pairs_to_keep // len(coor1_com)
    pairs = np.stack((pairs_col1, pairs_col2), axis=1)
    print("PAIRS AFTER EVAL", pairs, pairs.shape)

    clashes = []
    for n in range(0, len(pairs), MAX_CHUNK_2):
        pairs_chunk = pairs[n:n+MAX_CHUNK_2]
        print("2 CHUNK", n, n+len(pairs_chunk), len(pairs), file=sys.stderr)
        #if filtering by nucleotide - nucleotide COMs distances:
        #clashes_chunks = detect_clashes2(coor1, coor2, pairs_chunk)
        clashes_chunk = detect_clashes_all(coor1, coor2, pairs_chunk)
        if clashes_chunk is not None:
            clashes.append(clashes_chunk)
    if len(clashes):
        result = np.concatenate(clashes)
        return result
    else:
        return None

def detect_clashes(coor1, coor2, pairs):
    clashes0 = []
    print("detect cash %i vs %i"%(len(coor1), len(coor2)), file=sys.stderr)
    for n1 in range(0, len(coor1), MAX_CHUNK_1):
        coor_chunk1 = coor1[n1:n1+MAX_CHUNK_1]
        chunk_pairs1 = pairs.copy()
        column = chunk_pairs1[:,0]
        column -= n1
        chunk_pairs1_keep = (column >= 0) & (column < len(coor_chunk1))
        chunk_pairs1 = chunk_pairs1[chunk_pairs1_keep]        
        for n2 in range(0, len(coor2), MAX_CHUNK_1):
            coor_chunk2 = coor2[n2:n2+MAX_CHUNK_1]
            chunk_pairs2 = chunk_pairs1.copy()
            column = chunk_pairs2[:, 1]
            column -= n2
            chunk_pairs2_keep = (column >= 0) & (column < len(coor_chunk2))
            chunk_pairs2 = chunk_pairs2[chunk_pairs2_keep]
            pairs_to_eval = chunk_pairs2
            print("1 CHUNK", n1, n1+len(coor_chunk1), n2, n2+len(coor_chunk2), file=sys.stderr)
            clashes_chunk = detect_clashes_COM(coor_chunk1, coor_chunk2, pairs_to_eval, homo=False)
            if clashes_chunk is not None:
                print("CHUNK!", clashes_chunk.shape, file=sys.stderr)
                clashes_chunk2 = clashes_chunk + (n1, n2)
                if len(clashes_chunk2):
                    #print("CHUNK!!", clashes_chunk2.shape, file=sys.stderr)
                    clashes0.append(clashes_chunk2)
        ###
        pc = len(coor1) * len(coor2) // (n1+MAX_CHUNK_1)
        print("Done %.0f %%"%(100*pc), file=sys.stderr)
        ###

    clashes = []
    if len(clashes0):
        for chunk in clashes0:            
            clashes.extend(chunk.tolist())
    return clashes

def connection_1nucl(g):
    print("get connection_1nucl")
    connect1nucl = []
    nfrag = g["nfrags"]
    maps = []
    for n in range(nfrag-1):
        interactions = g['interactions-%i'%n] 
        pool = np.unique(interactions[:,0])
        curr_map = {pose:set() for pose in pool}
        for pose1, pose2 in interactions:
            curr_map[pose1].add(pose2)
        maps.append(curr_map)
    
    for n in range(nfrag-2):
        print("get connection_1nucl frag%i - frag%i"%(n+1, n+3))
        map1, map2 = maps[n], maps[n+1]        
        for pose1, poses2 in map1.items():
            curr_connect1nucl = set()
            for pose2 in poses2:
                poses3 = map2[pose2]
                curr_connect1nucl.update(poses3)
            connect1nucl.extend([(pose1, pose3) for pose3 in curr_connect1nucl])
 
    return np.array(connect1nucl)

g = np.load(args.connection_graph) #npz file from connect.py
nfrag = g["nfrags"]

# get pairs connected by 2 nucl
connect2nucl = [ (i[0], i[1]) for n in range(nfrag-2) for i in g['interactions-%i'%n]]
print("%i pairs connected by 2 nucl"%len(connect2nucl))

# get pairs connected by 1 nucl
connect1nucl = connection_1nucl(g)
print("%i pairs connected by 1 nucl"%len(connect1nucl))

## extract the coordinates of the terminal
# nucleotides that might clash (frag1_n1 and frag3_n3).
coordinates = np.load(args.coor) 
coor = coordinates.reshape((len(coordinates), -1, 3))

at1 = [int(i)-1 for i in args.atoms_n1]
at3 = [int(i)-1 for i in args.atoms_n3]

coor1 = coor[:, at1]
coor2 = coor[:, at3]

## detect clashes (frag1_n1 and frag3_n3).
clashes_map = np.array(detect_clashes(coor1, coor2, connect1nucl))

## remove poses clashing with themselves
print(clashes_map.shape)
mask = np.where(clashes_map[:, 0] != clashes_map[:, 1])
clashes_map_diff = clashes_map[mask]

## re-map clashes
clashes = np.stack((connect1nucl[clashes_map_diff[:,0], 0], connect1nucl[clashes_map_diff[:,1], 1]), axis=1)
print("%i pairs connected by 1 nucl that clash"%len(clashes))

## get pairs connected by 1 nucl and NOT clashing

### np.setdiff1d needs 1D arrays
### convert 2D arrays into 1D arrays by creating unique keys
### ex: [[52, 23], [9, 66] => [52023, 90066]
nmax = max(clashes_map[:, 1].max(), connect1nucl[:, 1].max())
clashes_key = nmax * clashes_map[:, 0] + clashes_map[:, 1]
connect1nucl_key = nmax * connect1nucl[:, 0] + connect1nucl[:, 1]

noclash_map_key = np.setdiff1d(connect1nucl_key, clashes_key, assume_unique = True) 

noclash_map_col0 = noclash_map_key // nmax
noclash_map_col1 = noclash_map_key % nmax
 
noclash = np.stack((connect1nucl[noclash_map_col0, 0], connect1nucl[noclash_map_col1, 1]), axis=1)

# Merge connected non-clashing pairs
connect = np.unique(np.concatenate((noclash, connect2nucl), axis=0), axis=0)

if args.npz is not None:
    np.savez(args.npz, connect)

if args.npy is not None:
    np.save(args.npy, connect)

if args.txt is not None:
    f = open(args.txt, "w")
    for p0, p1 in connect:
        print("%i %i"%(p0, p1), file=f)

if args.clashes is not None:
    print("%i clashes"%len(clashes))
    np.save(args.clashes, clashes)

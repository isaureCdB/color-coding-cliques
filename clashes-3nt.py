#!/usr/bin/env python3

'''
Brute force detection of clashes in trinucleotide poses
NOTE: for now, assumes homo trinucleotide poses!
(real version would take an extra three numbers which is the number of (all!) atoms per nucleotide)

In addition, a fixed radius for all heavy atoms is assumed, which isn't quite correct:
  The vdw radii are ~1.5 A for nitrogen/oxygen and ~1.75 A for carbon/phosphorus
Here, a clash threshold of 3.5 A is assumed
'''
import numpy as np
import argparse
import sys
import random
import opt_einsum
from scipy.spatial import KDTree

#######################################
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('coor', help="atomic coordinates in numpy format")
parser.add_argument('--ignore', help="list of pairs to ignore, npy format", nargs='+')
parser.add_argument('--dist', help="clashing distance in A", default=3.5)
parser.add_argument('--mem', help="mem usage (bytes) ", default=10**8)
parser.add_argument('--npz', help="output in npz format")
parser.add_argument('--npy', help="output in npy format")
parser.add_argument('--txt', help="output in txt format")
args = parser.parse_args()
#######################################

CLASH_THRESHOLD = float(args.dist)
MAX_CHUNK_MEM = int(args.mem)

MAX_CHUNK_1 = int(np.sqrt(MAX_CHUNK_MEM / 8 / 3))
MAX_CHUNK_2 = int(MAX_CHUNK_MEM / 8 / 3 / 9)

coordinates_aa = np.load(args.coor) 
coor = coordinates_aa.reshape((len(coordinates_aa), -1, 3))

# Add the symetric pairs to the list of pairs to ignore
if args.ignore:
    ignore_lists = [np.load(i) for i in args.ignore if len(np.load(i))>0]
    ignore_asym = np.concatenate( ignore_lists , axis=0)
    ignore_inverted = np.stack((ignore_asym[:,1], ignore_asym[:,0]), axis=1)
    ignore = np.concatenate( (ignore_asym, ignore_inverted) , axis=0)
else:
    ignore = None

trinucsize = coor.shape[1]
MAX_CHUNK_3 = int(MAX_CHUNK_MEM / 8 / 3 / trinucsize**2)

# debug: smaller chunk sizes
###MAX_CHUNK_1, MAX_CHUNK_2, MAX_CHUNK_3 = 50, 1000, 1000 ###
#print(MAX_CHUNK_1, MAX_CHUNK_2, MAX_CHUNK_3, file=sys.stderr)

assert trinucsize % 3 == 0  # must be homo-trinucleotide!
nucsize = trinucsize // 3

coor_per_nuc = coor.reshape(len(coor), 3, nucsize, 3)
coor_com_nuc = coor_per_nuc.mean(axis=2)
dif_com_nuc = (coor_per_nuc - coor_com_nuc[:, :, None, :])
nucleotide_radius = np.sqrt((dif_com_nuc**2).sum(axis=3).max())

def detect_clashes3(coor1, coor2, pairs):
    c1, c2 = coor1[pairs[:, 0]], coor2[pairs[:, 1]]  
    dif0 = c1[:, :, None, :] - c2[:, None, :, :]
    distsq = np.einsum("ijkl,ijkl->ijk", dif0, dif0)  # = (dif0*dif0).sum(axis=3)
    min_distsq = distsq.reshape(len(distsq), -1).min(axis=1)
    clashing = (min_distsq < CLASH_THRESHOLD**2)
    clashes = pairs[clashing]
    #print("CLASH!", len(clashes), len(pairs), file=sys.stderr)
    if not len(clashes):
        return None
    return clashes

def detect_clashes2(coor1, coor2, pairs):
    #retain pairs with nucleotides that have close enough COMs to clash
    coor1_per_nuc = coor1.reshape(len(coor1), 3, nucsize, 3)
    coor1_com_nuc = coor1_per_nuc.mean(axis=2)
    coor2_per_nuc = coor2.reshape(len(coor2), 3, nucsize, 3)
    coor2_com_nuc = coor2_per_nuc.mean(axis=2)
    c1, c2 = coor1_com_nuc[pairs[:, 0]], coor2_com_nuc[pairs[:, 1]]  
    dif0 = c1[:, :, None, :] - c2[:, None, :, :]
    distsq = np.einsum("ijkl,ijkl->ijk", dif0, dif0)  # = (dif0*dif0).sum(axis=3)
    min_distsq = distsq.reshape(-1, 9).min(axis=1)
    clash_threshold_nuc = CLASH_THRESHOLD + 2 * nucleotide_radius
    pairs_to_keep = pairs[(min_distsq < clash_threshold_nuc**2)]
    # debug: keep all pairs
    ###pairs_to_keep = pairs ###
    clashes = []
    for n in range(0, len(pairs_to_keep), MAX_CHUNK_3):
        pairs_chunk = pairs_to_keep[n:n+MAX_CHUNK_3]
        #print("3 CHUNK", n, n+len(pairs_chunk), file=sys.stderr)
        clashes_chunk = detect_clashes3(coor1, coor2, pairs_chunk)
        if clashes_chunk is not None:
            clashes.append(clashes_chunk)
    return clashes

def detect_clashes1(coor1, coor2, homo, pairs_to_ignore=None):
    #retain pairs that have close enough COMs to clash
    #pairs to ignore: counting from zero! 

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

    #if coor1 and coor2 are the same chunk
    if homo:
        keep = (pairs[:, 1] - pairs[:, 0] > 0)
        pairs = pairs[keep]

    # filter pairs based on ignore list
    if pairs_to_ignore is not None:
        assert pairs_to_ignore[:, 0].min() >= 0, pairs_to_ignore[:, 0].min()
        assert pairs_to_ignore[:, 1].min() >= 0, pairs_to_ignore[:, 1].min()
        assert pairs_to_ignore[:, 0].max() < len(coor1_com), (pairs_to_ignore[:, 0].max(), len(coor1_com))
        assert pairs_to_ignore[:, 1].max() < len(coor2_com), (pairs_to_ignore[:, 1].max(), len(coor2_com))

        pairs_key = len(coor1_com) * pairs[:, 1] + pairs[:, 0]
        pairs_key = np.sort(pairs_key)
        ignore_key = len(coor1_com) * pairs_to_ignore[:, 1] + pairs_to_ignore[:, 0]
        ignore_key = np.sort(ignore_key)
        pairs_to_keep = np.setdiff1d(pairs_key, ignore_key) #, assume_unique=True)
        print("PAIRS BEFORE IGNORE", pairs, pairs.shape, len(coor1_com))
        pairs_col1 = pairs_to_keep % len(coor1_com)
        pairs_col2 = pairs_to_keep // len(coor1_com)
        pairs_old = pairs ###
        pairs = np.stack((pairs_col1, pairs_col2), axis=1)
        print("PAIRS AFTER IGNORE", pairs, pairs.shape)

    clashes = []
    for n in range(0, len(pairs), MAX_CHUNK_2):
        pairs_chunk = pairs[n:n+MAX_CHUNK_2]
        print("2 CHUNK", n, n+len(pairs_chunk), len(pairs), file=sys.stderr)
        clashes_chunks = detect_clashes2(coor1, coor2, pairs_chunk)
        if clashes_chunks:
            clashes.extend(clashes_chunks)
    if len(clashes):
        return np.concatenate(clashes)
    else:
        return None
   

def detect_clashes(coor, ignore=None):
    clashes0 = []
    chunk_ignore = None
    for n1 in range(0, len(coor), MAX_CHUNK_1):
        coor_chunk1 = coor[n1:n1+MAX_CHUNK_1]
        chunk_ignore1 = None
        if ignore is not None:
            chunk_ignore1 = ignore.copy()
            column = chunk_ignore1[:,0]
            column -= n1
            chunk_ignore1_keep = (column >= 0) & (column < len(coor_chunk1))
            chunk_ignore1 = chunk_ignore1[chunk_ignore1_keep]
        for n2 in range(n1, len(coor), MAX_CHUNK_1):
            coor_chunk2 = coor[n2:n2+MAX_CHUNK_1]
            chunk_ignore2 = None
            if ignore is not None:
                chunk_ignore2 = chunk_ignore1.copy()
                column = chunk_ignore2[:, 1]
                column -= n2
                chunk_ignore2_keep = (column >= 0) & (column < len(coor_chunk2))
                chunk_ignore2 = chunk_ignore2[chunk_ignore2_keep]
                pairs_to_ignore = chunk_ignore2
            else:
                pairs_to_ignore = None
                
            print("1 CHUNK", n1, n1+len(coor_chunk1), n2, n2+len(coor_chunk2), file=sys.stderr)
            clashes_chunk = detect_clashes1(coor_chunk1, coor_chunk2, homo=(n1==n2), pairs_to_ignore=pairs_to_ignore)
            if clashes_chunk is not None:
                #print("CHUNK!", clashes_chunk.shape, file=sys.stderr)
                clashes_chunk2 = clashes_chunk + (n1, n2)
                if len(clashes_chunk2):
                    #print("CHUNK!!", clashes_chunk2.shape, file=sys.stderr)
                    clashes0.append(clashes_chunk2)
    clashes = []
    if len(clashes0):
        for chunk in clashes0:            
            clashes.extend(chunk.tolist())
    return clashes

# positive control
random.seed(0)
sample = random.sample(range(len(coor)), min(100,len(coor)) )
coor0 = coor[sample]
tree = KDTree(coor0.reshape(-1,3))
p = tree.query_pairs(CLASH_THRESHOLD,output_type="ndarray")
p //= trinucsize
p = p[p[:, 0] - p[:, 1] != 0]
p2 = 1000 * p[:,0] + p[:,1]
p2 = p2.astype(np.int64) 
truepos0 = np.unique(p2)
truepos = np.stack([truepos0//1000, truepos0 % 1000],axis=1)

print("Positive control (sample of 100)", file=sys.stderr)
print("True clashes:", len(truepos), file=sys.stderr)
test = detect_clashes(coor0)
print("Calculated clashes:", len(test), file=sys.stderr)

clashes = detect_clashes(coor, ignore)

ind1 = np.argsort(clashes[:, 1])
sorted = clashes[ind1]
ind0 = np.argsort(sorted[:, 0], kind="stable")
clashes = sorted[ind0]

print(clashes[0])

if args.npz:
    np.savez(args.npz, clashes)

if args.npy:
    np.save(args.npy, clashes)

if args.txt:
    clashes1 = np.array(clashes) + 1
    f=open(args.txt,"w")
    for p1, p2 in clashes1: 
        print("%i %i"%(p1,p2), file = f)

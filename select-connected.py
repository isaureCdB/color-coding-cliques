#!/usr/bin/env python3

import numpy as np
import argparse

'''
extract the list of unique connected poses
'''

#######################################
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('connection_graph', help="graph of poses connectivity, in npz")
parser.add_argument('--start', help="start indices")
args = parser.parse_args()
#######################################

npz = np.load(args.connection_graph)
nfrags = npz["nfrags"]
poses = set()
for n in range(nfrags-1):
    inter = npz["interactions-%d"%n]
    poses.update(inter[:,0])

inter = npz["interactions-%d"%(nfrags-2)]
poses.update(inter[:,1])

if args.start:
    start = int(args.start)

for p in sorted(list(poses)):
    if start:
        p += start
    print(p)

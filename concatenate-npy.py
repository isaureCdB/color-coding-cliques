#!/usr/bin/env python3

import numpy as np
import argparse

#######################################
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('outp', help="output in npy format")
parser.add_argument('inp', help="inputs in npy format", nargs='+')
parser.add_argument('--axis', default=0)
parser.add_argument('--sort', action="store_true")
parser.add_argument('--unique', action="store_true")
args = parser.parse_args()
#######################################

def sortnp(a): 
    ind1 = np.argsort(a[:, 1]) 
    s = a[ind1] 
    ind0 = np.argsort(s[:, 0], kind="stable") 
    return s[ind0] 


merge_list = [np.load(i) for i in args.inp if len(np.load(i)) > 0]
for m in merge_list:
    print(np.shape(m))

merged = np.concatenate(merge_list , axis=args.axis)
print("merged")
print(np.shape(merged))

result = merged
if args.sort:
    result = sortnp(merged)

if args.unique:
    unique = np.unique(result, axis=0)
    result = unique

    print("unique")
    print(np.shape(unique))

np.save(args.outp, result)
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

merge_lists = [np.load(i) for i in args.inp if len(np.load(i)) > 0]
merged = np.concatenate( merge_lists , axis=args.axis)

if args.sort:
    merged = np.sort(merged, axis=0)

if args.unique:
    merged = np.unique(merged, axis=0)

np.save(args.outp, merged)
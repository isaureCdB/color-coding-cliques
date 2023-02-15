#!/usr/bin/env python3

import numpy as np
import argparse

#######################################
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('outp', help="output in npy format")
parser.add_argument('inp', help="inputs in npy format", nargs='+')
parser.add_argument('--axis', help="array axis", default=0)
args = parser.parse_args()
#######################################

merge_lists = [np.load(i) for i in args.inp if len(np.load(i)) > 0]
merged = np.concatenate( merge_lists , axis=args.axis)
np.save(args.outp, merged)
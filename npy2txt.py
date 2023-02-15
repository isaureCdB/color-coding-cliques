#!/usr/bin/env python3

import numpy as np
import argparse

#######################################
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('inp', help="inputs in npy format")
parser.add_argument('--offset')
args = parser.parse_args()
#######################################

inp = np.load(args.inp)
for i in inp:
    if args.offset:
        i += args.offset
    print(i)
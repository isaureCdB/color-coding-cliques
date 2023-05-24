#!/usr/bin/env python3

import numpy as np
from typing import Tuple

def _select_element(matrix, kt, random_number):
    sums = matrix.sum(axis=1)
    # only consider sums > 0 and above -20 kt from the maximum
    threshold = max(sums.max() - 20 * kt, 0)
    mask = (sums > threshold)
    if mask.sum() == 0:
        return None
    probs = np.exp((sums[mask]-threshold)/kt)
    probs /= probs.sum()
    probs /= probs.sum()
    cum_probs = np.cumsum(probs)
    element0 = np.searchsorted(cum_probs, random_number, "left")
    return np.nonzero(mask)[0][element0]
       
def largest_clique(matrix:np.ndarray, kt:float) -> Tuple[int]:
    """finds the largest clique heuristically
    kt: temperature.
      two elements with number-of-connections X and X+kt have 
       a factor e (2.71) difference in probability to be selected
    """
    assert matrix.ndim == 2 and matrix.shape[0] == matrix.shape[1]
    size = len(matrix)
    random_numbers = np.random.random(size)

    m = matrix.copy()
    clique = []
    for n in range(size):
        element = _select_element(m, kt, random_numbers[n])
        if element is None:
            break
        element_sum = m[element].sum()
        #print(n, "%2d" % element, "%2d" % element_sum)
        eliminated_partners = (~m[element]).copy()
        m[eliminated_partners] = 0
        m[:, eliminated_partners] = 0
        m[element] = 0
        clique.append(element)
    return tuple(sorted(clique))
    
if __name__ == "__main__":
    import sys
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("connection_list", help="connection list in Numpy format")
    parser.add_argument("--cliques", type=int, required=True, help="Number of cliques to find")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--kt", type=float, required=True, help="""Temperature.
    two elements with number-of-connections X and X+kt have 
    a factor e (2.71) difference in probability to be selected"""
    )                
    parser.add_argument("--trials", type=int, required=True, help="Number of probabilistic trials")
    args = parser.parse_args()
    np.random.seed(args.seed)
    connection_list = np.load(args.connection_list)
    assert connection_list.ndim == 2, connection_list.shape
    assert connection_list.shape[1] == 2, connection_list.shape
    assert connection_list.dtype == int, connection_list.dtype
    assert connection_list.min() >= 0, connection_list.min()

    size = connection_list.max() + 1
    
    matrix = np.eye(size, size, dtype=bool) #force diagonal to be positive

    matrix[connection_list[:, 0], connection_list[:, 1]] = 1

    # force matrix to be symmetric
    triu_indices = np.triu_indices(size, 1, size)
    matrix.T[triu_indices] = matrix[triu_indices]

    print("Number of poses: {}".format(size), file=sys.stderr)
    print("Number of clashes: {}".format(matrix.sum()), file=sys.stderr)

    m = matrix.copy()
    for clique in range(args.cliques):
        maxcliquesize = 0
        for n in range(args.trials):
            curr_maxclique = largest_clique(m, args.kt)
            if len(curr_maxclique) > maxcliquesize:
                maxclique = curr_maxclique
                maxcliquesize = len(maxclique)
            print(clique +1, n+1, "Probabilistic maximum clique", len(curr_maxclique), file=sys.stderr)
        if maxcliquesize == 0:
            break
        print(maxclique)
        m[list(maxclique)] = 0
        m[:, list(maxclique)] = 0
# color-coding-cliques
creates cliques for color-coding assembly of docking poses of ssRNA 3-nucleotides

The main script is cliques-of-clashing-poses.sh

-------------- Variables: \
\
graph=$1   # graph of poses connexions in npz format\
d=$2       # clashing distance (heavy atoms, in A)\
coor=$3    # all-atom (heavy) coordinates of the connected poses, in npy format shaped as (pose, atom, coor)\
nat1=$4    # number of heavy atoms in the 1st nucleotide\
nat2=$5    # number of heavy atoms in the 2nd nucleotide\
nat3=$6    # number of heavy atoms in the 3rd nucleotide\


-------------- required files: \
\

- graph of poses connections, in npz format (numpy dictionary), with indices of poses corresponding to their docking rank from 0

- coordinates of all (and only) the connected poses in all-atom, in npy format shaped as (pose, atom, coor). Poses must be sorted by their docking rank.



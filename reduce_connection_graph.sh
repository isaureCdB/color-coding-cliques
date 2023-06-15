# create graph and coordinates with only connected poses


graph=$1      # graph of poses connexions in npz format
coor=$2				# coarse-grained coordinate of all poses, in npy format (pose, atom, coor)
o=$3          # overlapping max dist (in A) of compatible poses. Recommanded in range [1.0-2.0]]
#nposes=$4 		# max number of connected poses per fragments(for memory saving). Recommanded <100000
motif=$5      # 3-nt sequence of the fragments
nfrag=$6      # Number of assembled fragments in the graph
nat1=$7    #number of pseudo-atoms in the 1st nucleotide
nat2=$8    #number of pseudo-atoms in the 2nd nucleotide
nat3=$9    #number of pseudo-atoms in the 3rd nucleotide

name=${graph%%.npz}

## !!!!! for $name.npz indices from 1 !!!!!
## if $name.npz indices are from 1, change to --offset 0
$COLCO/select-connected.py $name.npz --offset 1 > $name-connected.list

# $motif-e6-aa.npy = all-atom coordinates of the docked poses
$COLCO/select-npy.py $coor $name-connected-aa.npy --f $name-connected.list

$COLCO/select-npy.py $coor preatoms.npy  --f $name-connected.list --atom `seq $(($nat1+1)) $(($nat1+$nat2+$nat3))` # coarse-grained coor of 2 last nucl
$COLCO/select-npy.py $coor postatoms.npy --f $name-connected.list --atom `seq 1 $(($nat1+$nat2))` # coarse-grained coor of 2 first nucl

# create graph with only those poses
$COLCO/connect-npz.py 2 $o 1000000 100 2frag-connected.npz preatoms.npy preatoms.npy postatoms.npy postatoms.npy
$COLCO/connect-homo-npz.py 2frag-connected.npz $nfrag $name-connected.npz >& $name-connected-NPZ.log

# count total number of chains
count_chains_npz.py  $name.npz 
count_chains_npz.py  $name.npz 
echo "check that you got twice the same number of chains"

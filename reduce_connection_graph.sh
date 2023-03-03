# create graph and coordinates with only connected poses


s="./color-coding-cliques" #scripts repo

o=1.8         # overlapping RMSD (for overlapping 3-nt)
nposes=100000 # max number of connected poses (for memory saving)
motif=AAA     # 3-nt sequence of the fragments
nfrag=5       # Number of assemblde fragments in the graph

graph=$motif-${nfrag}frag-${o}A

## !!!!! for $graph.npz indices from 1 !!!!!
## if $graph.npz indices are from 1, change to --offset 0
$s/select-connected.py $graph.npz --offset 1 > $graph-connected.list

# $motif-e6-aa.npy = all-atom coordinates of the docked poses
$s/select-npy.py ../$motif-e6-aa.npy $graph-connected-aa.npy --f $graph-connected.list

$s/select-npy.py ../$motif-e6-preatoms.npy preatoms.npy   --f $graph-connected.list # coarse-grained coor for 2 last nucl
$s/select-npy.py ../$motif-e6-postatoms.npy postatoms.npy --f $graph-connected.list # coarse-grained coor for 2 first nucl

# create graph with only those poses
$s/connect-npz.py 2 $o $nposes 100 2frag-connected.npz preatoms.npy preatoms.npy postatoms.npy postatoms.npy
$s/connect-homo-npz.py 2frag-connected.npz $nfrag $graph-connectedB.npz >& $graph-connected-NPZ.log

# count total number of chains
count_chains_npz.py  $graph.npz 
count_chains_npz.py  $graph.npz 
echo "check that you got twice the same number of chains"

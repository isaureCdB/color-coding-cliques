# !!! numbering from 0 in npy files, and from 1 in txt files !!!

s="./color-coding-cliques" #scripts repo
d=3.5
o=1.8
nposes=100000
motif=AAA
nfrag=5
graph=$motif-${nfrag}frag-${o}A

# create graph and coordinates with only connected poses

# !!!!! for $graph.npz indices from 1 !!!!!
count_chains_npz.py  $graph.npz 
$s/select-connected.py $graph.npz --offset 1 > $graph-connected.list

# !!!!! if $graph.npz indices from 1 : !!!!!
#$s/select-connected.py $graph.npz --offset 0 > $graph-connected.list

#$motif-e6-aa.npy = coordinates of the docked poses
#$s/select-npy.py $motif-e6-aa.npy $graph-connected-aa.npy --f $graph-connected.list

$s/select-npy.py ../$motif-e6-preatoms.npy preatoms.npy   --f $graph-connected.list
$s/select-npy.py ../$motif-e6-postatoms.npy postatoms.npy --f $graph-connected.list

$s/connect-npz.py 2 $o $nposes 100 2frag-connected.npz preatoms.npy preatoms.npy postatoms.npy postatoms.npy
$s/connect-homo-npz.py 2frag-connected.npz $nfrag $graph-connectedB.npz >& $graph-connected-NPZ.log

# get list "samepos" of pairs of poses seen only at one same position in chain
# (those can't be together in a chain)
$s/same-position.py $graph-connected.npz --npy $graph-samepos.npy 
$s/npy2txt.py $graph-samepos.npy --offset 1 > $graph-samepos.txt

# get list "noclash" of non-clashing overlapping pairs of poses:
# - pairs overlapping by 2 nucleotides (frag1_n2n3 frag2_n1n2)
# - pairs overlapping by 1 nucl (frag1_n3 frag3_n1) and not clashing (frag1_n1 frag3_n3)
$s/overlapping-noclash.py $graph-connected.npz $graph-connected-aa.npy --atoms_n1  `seq 1 22` --atoms_n3 `seq 45 66` \
    --npy $graph-noclash$d.npy  --clashes $graph-clash$d-n1n3.npy --dist $d

$s/npy2txt.py $graph-noclash$d.npy --offset 1 > $graph-noclash$d.txt
$s/npy2txt.py $graph-clash$d-n1n3.npy --offset 1 > $graph-clash$d-n1n3.txt

#get list L3 of other clashing pairs of poses (excluding "noclash" and "samepos")
$s/clashes-3nt.py $graph-connected-aa.npy --dist $d --ignore $graph-samepos.npy $graph-noclash$d.npy \
    --npy $graph-clash$d.npy

$s/npy2txt.py  $graph-clash$d.npy --offset 1 > $graph-clash$d.txt

# get list of incompatible poses:
#   - Poses in samepos are incompatible
#   - Poses in L3 are incompatible
$s/concatenate-npy.py $graph-clash$d-incompatible.npy $graph-samepos.npy $graph-clash$d.npy --sort --unique

$s/npy2txt.py  $graph-clash$d.npy --offset 1 > $graph-clash$d.txt

#create cliques of incompatible poses:
#TODO

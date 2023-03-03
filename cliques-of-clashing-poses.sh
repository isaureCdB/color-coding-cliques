# !!! numbering from 0 in npy files, and from 1 in txt files !!!

s="./color-coding-cliques" #scripts repo

d=1.5         # clashing distance (A, all-atom)
o=1.8         # overlapping RMSD (for overlapping 3-nt)
nposes=100000 # max number of connected poses (for memory saving)
motif=AAA     # 3-nt sequence of the fragments
nfrag=5       # Number of assemblde fragments in the graph

graph=$motif-${nfrag}frag-${o}A

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

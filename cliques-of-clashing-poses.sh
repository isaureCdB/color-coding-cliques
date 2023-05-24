# !!! numbering from 0 in npy files, and from 1 in txt files !!!

#export COLCO={path to scripts repo}

d=1.5         # clashing distance (A, all-atom)
o=1.8         # overlapping RMSD (for overlapping 3-nt)
nposes=100000 # max number of connected poses (for memory saving)
motif=AAA     # 3-nt sequence of the fragments
nfrag=5       # Number of assembled fragments in the graph

graph=$motif-${nfrag}frag-${o}A

# get list "samepos" of pairs of poses seen only at one same position in chain
# (those can't be together in a chain)
$COLCO/same-position.py $graph-connected.npz --npy $graph-samepos.npy 
$COLCO/npy2txt.py $graph-samepos.npy --offset 1 > $graph-samepos.txt

# get list "noclash" of non-clashing overlapping pairs of poses:
# - pairs overlapping by 2 nucleotides (frag1_n2n3 frag2_n1n2)
# - pairs overlapping by 1 nucl (frag1_n3 frag3_n1) and not clashing (frag1_n1 frag3_n3)
$COLCO/overlapping-noclash.py $graph-connected.npz $graph-connected-aa.npy --atoms_n1  `seq 1 22` --atoms_n3 `seq 45 66` \
    --npy $graph-noclash$d.npy  --clashes $graph-clash$d-n1n3.npy --dist $d

$COLCO/npy2txt.py $graph-noclash$d.npy --offset 1 > $graph-noclash$d.txt
$COLCO/npy2txt.py $graph-clash$d-n1n3.npy --offset 1 > $graph-clash$d-n1n3.txt

#get list L3 of other clashing pairs of poses (excluding "noclash" and "samepos")
$COLCO/clashes-3nt.py $graph-connected-aa.npy --dist $d --ignore $graph-samepos.npy $graph-noclash$d.npy \
    --npy $graph-clash$d.npy

$COLCO/npy2txt.py  $graph-clash$d.npy --offset 1 > $graph-clash$d.txt

# get list of incompatible poses:
#   - Poses in samepos are incompatible
#   - Poses in L3 are incompatible
$COLCO/concatenate-npy.py $graph-clash$d-incompatible.npy $graph-samepos.npy $graph-clash$d.npy --sort --unique

$COLCO/npy2txt.py  $graph-clash$d.npy --offset 1 > $graph-clash$d.txt

#create cliques of incompatible poses:
python3 clique.py $graph-clash$d-incompatible.npy --seed 6 --kt 1 --trials 10 --cliques 200 > $graph-clash$d-cliques.txt

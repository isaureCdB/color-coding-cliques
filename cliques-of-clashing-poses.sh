# !!! numbering from 0 in npy files, and from 1 in txt files !!!

#export COLCO={path to scripts repo}

graph=$1  #graph of poses connexions in npz format
d=$2      #clashing distance (heavy atoms)
coor=$3    #all-atom (heavy) coordinates of the connected poses, in npy format shaped as (pose, atom, coor).
nat1=$4     #number of heavy atoms in the 1st nucleotide
nat2=$5     #number of heavy atoms in the 2nd nucleotide
nat3=$6     #number of heavy atoms in the 3rd nucleotide

name=${graph%%.npz}

# get list "samepos" of pairs of poses seen only at one same position in chain
# (those can't be together in a chain)
$COLCO/same-position.py $name.npz --npy $name-samepos.npy 
$COLCO/npy2txt.py $name-samepos.npy --offset 1 > $name-samepos.txt

# get list "noclash" of non-clashing overlapping pairs of poses:
# - pairs overlapping by 2 nucleotides (frag1_n2n3 frag2_n1n2)
# - pairs overlapping by 1 nucl (frag1_n3 frag3_n1) and not clashing (frag1_n1 frag3_n3)
$COLCO/overlapping-noclash.py $name.npz $coor --atoms_n1  `seq 1 $nat1` --atoms_n3 `seq $(($nat1+$nat2+1)) $(($nat1+$nat2+$nat3))` \
    --npy $name-noclash$d.npy  --clashes $name-clash$d-n1n3.npy --dist $d

$COLCO/npy2txt.py $name-noclash$d.npy --offset 1 > $name-noclash$d.txt
$COLCO/npy2txt.py $name-clash$d-n1n3.npy --offset 1 > $name-clash$d-n1n3.txt

#get list L3 of other clashing pairs of poses (excluding "noclash" and "samepos")
$COLCO/clashes-3nt.py $coor --dist $d --ignore $name-samepos.npy $name-noclash$d.npy \
    --npy $name-clash$d.npy

$COLCO/npy2txt.py  $name-clash$d.npy --offset 1 > $name-clash$d.txt

# get list of incompatible poses:
#   - Poses in samepos are incompatible
#   - Poses in L3 are incompatible
$COLCO/concatenate-npy.py $name-clash$d-incompatible.npy $name-samepos.npy $name-clash$d.npy --sort --unique

$COLCO/npy2txt.py  $name-clash$d.npy --offset 1 > $name-clash$d.txt

#create cliques of incompatible poses:
python3 $COLCO/clique.py $name-clash$d-incompatible.npy --seed 6 --kt 1 --trials 10 --cliques 200 > $name-clash$d-cliques.txt

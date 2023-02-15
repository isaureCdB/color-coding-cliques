# ${motif} -5frag -5A=> ${motif}

# !!! Any printed/saved numbering starts at 1 !!!

d=5.0
cutoff=7
nposes=10
motif=AAA
nfrag=5
graph=${motif}-${nfrag}frag-${cutoff}A

################# to test on subset
python3 ~/Scripts/connect-npz.py 2 $cutoff $nposes 100 ${motif}-2frag-${cutoff}A.npz \
    ../../$motif-e6-preatoms.npy ../../$motif-e6-preatoms.npy ../../$motif-e6-postatoms.npy ../../$motif-e6-postatoms.npy

python3 ~/Scripts/connect-homo-npz.py ${motif}-2frag-${cutoff}A.npz $nfrag $graph.npz >& $graph-NPZ.log
#################

# create graph and coordinates with only connected poses
~/Scripts/select_connected.py $graph.npz --start 1 > $graph-connected.list
~/Scripts/select-npy.py ${motif}.npy  $graph-connected.npy --f $graph-connected.list

~/Scripts/select-npy.py ../../$motif-e6-preatoms.npy preatoms.npy   --f $graph-connected.list
~/Scripts/select-npy.py ../../$motif-e6-postatoms.npy postatoms.npy --f $graph-connected.list

~/Scripts/connect-npz.py 2 $cutoff $nposes 100 2frag-connected.npz preatoms.npy preatoms.npy postatoms.npy postatoms.npy
~/Scripts/connect-homo-npz.py 2frag-connected.npz $nfrag $graph-connected.npz >& $graph-connected-NPZ.log

# get list L1 of pairs of poses seen only at one same position in chain
# (those can't be together in a chain)
./scripts/same-position.py $graph-connected.npz --npy $graph-L1.npy --txt $graph-L1.txt

# get list L2 of non-clashing overlapping pairs of poses:
# - pairs overlapping by 2 nucleotides (frag1_n2n3 frag2_n1n2)
# - pairs overlapping by 1 nucl (frag1_n3 frag3_n1) and not clashing (frag1_n1 frag3_n3)
./scripts/overlapping-noclash.py $graph-connected.npz $graph-connected.npy --atoms_n1  `seq 1 7` --atoms_n3 `seq 15 21` \
    --npy $graph-clash$d-L2.npy --txt $graph-clash$d-L2.txt --clashes $graph-clash$d-n1n3.npy

#get list L3 of other clashing pairs of poses (excluding L2 and L1)
./scripts/clashes-3nt.py $graph-connected.npy --dist $d --ignore $graph-L1.npy $graph-clash$d-L2.npy \
    --npy $graph-clash$d-L3.npy --txt $graph-clash$d-L3.txt

# get list of incompatible poses:
#   - Poses in L1 are incompatible
#   - Poses in L3 are incompatible
./scripts/concatenate-npy.py $graph-clash$d-incompatible.npy $graph-L1.npy $graph-clash$d-L3.npy 
    
#create cliques of incompatible poses:
#TODO

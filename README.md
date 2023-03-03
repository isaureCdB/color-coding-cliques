# color-coding-cliques
creates cliques for color-coding assembly of ssRNA 3nt

The main script is cliques-of-clashing-poses.sh

-------------- Variables:
d=1.5         # clashing distance (A, all-atom)
o=1.8         # overlapping RMSD (for overlapping 3-nt)
nposes=100000 # max number of connected poses (for memory saving)
motif=AAA     # 3-nt sequence of the fragments
nfrag=5       # Number of assemblde fragments in the graph

graph=$motif-${nfrag}frag-${o}A

-------------- Pre-existing files:
$graph-connected.npz = graph of poses connections, in npz format (numpy dictionary) 
graph-connected-aa.npy = coordinates of the connected poses in all-atom



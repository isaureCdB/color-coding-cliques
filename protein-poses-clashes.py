import sys
import numpy as np
from scipy.spatial import KDTree

protein = np.load(sys.argv[1])
assert protein.ndim == 2 and protein.shape[1] == 3, protein.shape

poses = np.load(sys.argv[2])
assert poses.ndim == 3 and poses.shape[2] == 3, poses.shape

clash_radius = float(sys.argv[3])

tree_protein = KDTree(protein)
tree_poses = KDTree(poses.reshape(-1, 3))

clashes = tree_poses.query_ball_tree(tree_protein, r=clash_radius)
nclashes = np.array([len(lis) for lis in clashes]).reshape(poses.shape[:2]).sum(axis=1)

for c in nclashes: print(c)
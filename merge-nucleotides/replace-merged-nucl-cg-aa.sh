convert-cg-aa.py mergednucl.npy mergednucl-aa.npy  \
    --libcg Ar.npy Ar.npy Ar.npy Ar.npy Ar.npy Ar.npy Ar.npy Ar.npy Ar.npy Ar.npy Ar.npy \
    --libaa A.npy A.npy A.npy A.npy A.npy A.npy A.npy A.npy A.npy A.npy A.npy


npy2pdb.py mergednucl-aa.npy A8-aa.pdb > mergednucl-aa.pdb

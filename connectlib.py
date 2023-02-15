copyright 2014-2017 Sjoerd de Vries, TUM

import sys, numpy as np, weakref, itertools

#variables for computation parallelisation
MAXCHUNK = 4000000

#CLUSTERING = list of RMSD cutoffs for hierarchical clustering
#second-to-last is deredundant criterion (only one of two poses at this distance will be kept,
# as they are considered as identical)

CLUSTERING = [10, 8, 6, 5, 4, 3.5, 3, 2.5, 2, 1.7, 1.5, 1.0, 0.5, 0.1, 0]
####CLUSTERING = [0.1, 0]

MAX_CLUSTERING = len(CLUSTERING) - 1

# CLUST_MARGIN * clust_radius is the margin to consider if one member of a cluster
# could overlap with a member of another cluster based on the RMSD of the cluster representatives.
# 2 is the theoretical worse case; change to 9999 to disable all effects of clustering
CLUST_MARGIN = 2
MINCHILD = 100

###################################################
try:
    import cffi
    import _get_msd
    ffi = _get_msd.ffi
    def get_msd(c1, c2):
        #actually computes square-deviation (sd), not msd
        def npdata(a):
            return a.__array_interface__["data"][0]

        nc1 = len(c1)
        nc2 = len(c2)
        natom = c1.shape[1]
        assert c2.shape[1] == natom
        msd = np.empty((nc1, nc2), dtype=float)
        _get_msd.lib.get_msd(
          nc1, nc2, natom,
          ffi.cast("double *", npdata(c1)),
          ffi.cast("double *", npdata(c2)),
          ffi.cast("double *", npdata(msd)),
        )
        return msd
except ImportError:
    print("Cannot find cffi, you will lose some speed", file=sys.stderr)
    def get_msd(c1, c2):
        #actually computes square-deviation (sd), not msd
        d = c1[:, np.newaxis, :, :] - c2[np.newaxis, :, :, :]
        # dsq = d * d
        # msd = dsq.sum(axis=3).sum(axis=2)
        msd = np.einsum("...ijk,...ijk->...i", d,d)
        return msd

# Depth-first decomposition
def decompose(clusters, clusnr, max_rmsd):
    best_conlevel = None # connection-level = priority of the cluster to be decomposed
                         # (here, correspond to depth)
    bestc = None # best cluster = cluster with highest priority
    clus = clusters[clusnr]
    for c in clus:
        if not len(c.children): continue
        if best_conlevel is None or c.conlevel > best_conlevel:
            best_conlevel = c.conlevel
            bestc = c
    if bestc is None: return False
    c = bestc
    clusters[clusnr].remove(c)
    if (clusnr % 2): # if preatoms
        c.decompose(max_rmsd, fwd=True)        # check connections with postatoms of forward/downstream frag
        c.decompose_intra(fwd=False) # check connections with postatoms of same frag
    else: # if postatoms
        c.decompose_intra(fwd=True) # check connections with preatoms of same frag
        c.decompose(max_rmsd, fwd=False)      # check connections with preatoms of backward/upstream frag
    c.reparent_children()
    for child in c.children:
        # additional ways to ensure depth search.
        conlevel = 0
        for con in itertools.chain(child.back_connections,child.connections):
            v = con.clusterlevel
            if v is not None and v > conlevel:
                conlevel = v
        child.conlevel = conlevel * child.clusterlevel
    clusters[clusnr].update(set(c.children))
    for cc in c.children:
        cc.check_deletion()
    return True

###################################################
class Cluster(object):
    __slots__ = ("clustid", "_splittable", "clusterlevel", "coors", "ranks", "all_ranks", "children", "nodes", "totstruc", "parent", "connections", \
      "back_connections", "_splitting", "_checking_delete", "conlevel", "_dead", "clusters")
    def __init__(self, clusters, clustid, clusterlevel, coors, ranks):
        """
        A cluster can be in three forms:
        - leaf form. The cluster does not have any child clusters (yet).
            .coors contains all coordinates of the cluster (the first one is the representative).
            .ranks contains all ATTRACT ranks of the cluster.
        - node form. The cluster contains child clusters.
            .coors contains the coordinates of each representative of each child cluster.
            .ranks contains the ATTRACT rank of each representative of each child cluster.
        - singleton form. The cluster contains only a single member.
        There is a fourth form with coors = None, but this is not currently used in the code.
        The cluster is initially in leaf form (or in singleton form).
        """
        self.clusters = clusters
        self._dead = False
        self.clustid = clustid # tuple
        # clustid of top-level clust of frag1-postatoms is (1,)
        #            2nd-level ........................ are (1,1), (1,2), ...
        #            top-level clust of frag1-preatoms  is (1001,)
        #            2nd-level ........................ are (1001,1), (1001,2)
        #            top-level ........ frag2-postatoms is (2,)
        #            top-level ........ frag2-preatoms  is (1002,)
        # ...
        self._splittable = True
        if coors is not None and len(coors) == 1: #singleton
            self.clusterlevel = MAX_CLUSTERING
        else:
            assert clusterlevel is None or clusterlevel < MAX_CLUSTERING
            self.clusterlevel = clusterlevel #contains the clusterlevel of the cluster itself, not of the cluster children!
        if self.clusterlevel == MAX_CLUSTERING:
            self._splittable = False
        self.coors = coors #coordinates of the representative of the children clusters (or all coors if this clust is a leaf)
        self.ranks = ranks #ranks of the representative of the children clusters (or all ranks if this clust is a leaf)
        self.all_ranks = set(ranks) #ranks of the representative of the children clusters
        self.children = []
        self.nodes = 1 # nb of children
        if coors is not None:
            self.totstruc = coors.shape[0]
        self.parent = None
        self.connections = []
        self.back_connections = []
        self._splitting = False
        self._checking_delete = False
        self.conlevel = self.clusterlevel
        if self.clusterlevel is None:
            self.conlevel = -1
        if coors is not None and clusterlevel == MAX_CLUSTERING - 1: # deredundant level (cf line 17)
            r = self
            for cnr in range(len(self.coors)):
                c = Cluster(self.clusters, self.clustid + (cnr,), MAX_CLUSTERING, coors[cnr:cnr+1], ranks[cnr:cnr+1])
                c.parent = r
                self.children.append(c)
            self.nodes = len(self.children)
            self._splittable = False

    def _cluster(self, clusterlevel):
        '''subdivide cluster'''
        assert not self.children
        c = self.coors
        #clus: coordinates of the representative of each cluster
        clus = c[:1]
        #clus_indices: the coors indices of the structures of each cluster
        clus_indices = [[0]]
        chunksize = 20

        radius = CLUSTERING[clusterlevel]
        max_sd = radius * radius * c.shape[1]

        #This variable keeps track, for every structure in the chunk, into which new cluster it is sorted
        which_new_clust = np.zeros(chunksize, dtype=int)

        clustid = self.clustid
        if clustid is None: clustid = ()
        for n in range(1, len(c), chunksize):
            chunk = c[n:n+chunksize]

            #intra-chunk msd
            d = chunk[:, np.newaxis, :, :] - chunk[np.newaxis, :, :, :]
            intra_msd = np.einsum("...ijk,...ijk->...i", d,d)

            #chunk-cluster msd: compare to existing reference structures
            d = chunk[:, np.newaxis, :, :] - clus[np.newaxis, :, :, :]
            inter_msd = np.einsum("...ijk,...ijk->...i", d,d)

            for nn in range(len(chunk)):
                sort_clust = None
                # identify sets of structures that should go into the same new cluster
                # if they do not go in an existing cluster
                # => avoid creating to many references
                close_intra_clusts = (intra_msd[nn] < max_sd).nonzero()[0]
                intra_new_clusts = [which_new_clust[k] for k in close_intra_clusts if k < nn and which_new_clust[k] != -1]
                if len(intra_new_clusts):
                    sort_clust = min(intra_new_clusts)
                close_inter_clusts = (inter_msd[nn] < max_sd).nonzero()[0]
                if len(close_inter_clusts):
                    sort_clust2 = min(close_inter_clusts)
                    if sort_clust is None or sort_clust > sort_clust2:
                        sort_clust = sort_clust2
                if sort_clust is None:
                    #new cluster
                    which_new_clust[nn] = len(clus)
                    clus = np.append(clus, chunk[nn][np.newaxis,:,:], axis=0)
                    clus_indices.append([n+nn])
                else:
                    clus_indices[sort_clust].append(n+nn)
                    which_new_clust[nn] = -1

        indices = [i[0] for i in clus_indices]

        # Re-cluster to the closest reference
        clus_indices = [i[:1] for i in clus_indices]
        for n in range(0, len(c), chunksize):
            chunk = c[n:n+chunksize]
            d = chunk[:, np.newaxis, :, :] - clus[np.newaxis, :, :, :]
            inter_msd = np.einsum("...ijk,...ijk->...i", d,d)
            sort_clusts = np.argmin(inter_msd, axis=1)
            for nn in range(len(chunk)):
                if (n+nn) in indices: continue
                sort_clust = sort_clusts[nn]
                clus_indices[sort_clust].append(n+nn)

        for cnr,c in enumerate(clus_indices):
            ind = clus_indices[cnr]
            coors = self.coors[ind]
            ranks = self.ranks[ind]
            c = Cluster(self.clusters, clustid+(cnr+1,), clusterlevel, coors, ranks)
            c.parent = self
            self.children.append(c)
        self.nodes = len(self.children)
        self.totstruc = sum([c.totstruc for c in self.children])
        return clus, indices

    def cluster(self, clusterlevel):
        """Converts a cluster in leaf form to a cluster in node form."""
        assert clusterlevel < MAX_CLUSTERING
        clus, indices = self._cluster(clusterlevel)
        self.coors = clus
        self.ranks = self.ranks[indices]
        r = self
        for c in self.children:
            c.parent = r
        for c in self.children:
            if c._splittable:
                break
        else:
            self._splittable = False

    def dissolve(self, clusterlevel):
        """Dissolve all direct child clusters, clustering them and linking all of their children to us"""
        self.ranks = np.concatenate([c.ranks for c in self.children])
        newchildren = []
        coors = []
        while len(self.children):
            child = self.children.pop()
            child.cluster(clusterlevel)
            coors.append(child.coors)
            newchildren += child.children
            self.totstruc = child.totstruc
        self.coors = np.concatenate(coors, axis=0)
        self.children = newchildren
        r = self
        for c in self.children:
            c.parent = r
        self.nodes = len(newchildren)
        for c in self.children:
            if c._splittable:
                break
        else:
            self._splittable = False

    def reorganize(self):
        """Prune superfluous levels of clustering by dissolving child clusters with less than MINCHILD children
        Evokes reorganize() also on our children """
        for c in self.children:
            c.reorganize()
        if not len(self.children): return
        if len(self.children) >= MINCHILD: return
        oldchildren = [c for c in self.children if len(c.children)]
        if not len(oldchildren): return
        while len(self.children) < MINCHILD and len(oldchildren):
            child = oldchildren.pop(0)
            self.children.remove(child)
            self.children += child.children
            oldchildren += [c for c in child.children if len(c.children)]
        r = self
        coors = [c.coors[0] for c in self.children]
        self.coors = np.array(coors)
        for c in self.children:
            c.parent = r
        self.nodes = sum([c.nodes for c in self.children])
        for c in self.children:
            if c._splittable:
                break
        else:
            self._splittable = False

    def split(self):
        """Sub-cluster ourselves and then our children, until we are a singleton or have more than 1 child
        Returns whether or not we are still splittable further"""
        assert self.parent is not None
        if not self.children:
            assert self.clusterlevel is not None
            if self.clusterlevel == MAX_CLUSTERING: return False
            self.clusterlevel += 1
            if self.clusterlevel == MAX_CLUSTERING:
                self._splittable = False
                return False
            self.cluster(self.clusterlevel)
            r = self
            while len(self.children) == 1:
                child = self.children[0]
                self.clusterlevel = child.clusterlevel
                if self.clusterlevel >= MAX_CLUSTERING - 1:
                    self._splittable = False
                    break
                self.clusterlevel += 1
                child.cluster(self.clusterlevel)
                children = child.children
                for c in children: c.parent = r
                self.coors = child.coors
                self.ranks = child.ranks
                self.children = children
                self.nodes = len(children)
                self.totstruc = sum([c.totstruc for c in children])

            self.parent.add_nodes(self.nodes - 1)
            for c in self.children:
                if c._splittable:
                    break
            else:
                self._splittable = False
            return True
        else:
            self._splitting = True
            oldnodes = self.nodes
            ok = False
            for c in self.children:
                c_splittable = c._splittable
                has_split = c.split()
                if has_split: ok = True
            newnodes = self.nodes
            self._splitting = False
            if self.parent is not None and newnodes > oldnodes:
                self.parent.add_nodes(newnodes-oldnodes)
            for c in self.children:
                if c._splittable:
                    break
            else:
                self._splittable = False
            return ok

    def add_nodes(self, nodes):
        self.nodes += nodes
        if self.parent is not None and not self._splitting:
            self.parent.add_nodes(nodes)

    def check_deletion(self):
        """Check if cluster is dead because of lack of connections
        If so, remove the cluster and propagate check_deletion along the tree
        """
        if self._dead: return
        if self._checking_delete: return
        self._checking_delete = True
        while 1:
            has_c1, has_c2 = len(self.back_connections), len(self.connections)
            if has_c1 and has_c2: break # Forward and backward connections, cluster is alive

            # If not, check position in the chain, to know if we need
            # both forward and backward connections
            frag0 = self.clustid[0]
            if frag0 > 1000:
                pos = 2 * (frag0-1001) + 1
            else:
                pos = 2 * (frag0-1)
            if pos == 0 and has_c2: break # we are at the beginning of chain
                                          # and we do have fwd connections
            if pos == len(self.clusters) - 1 and has_c1: break # we are at the end of chain
                                                          # and we do have bwd connections

            #Cluster is dead. Remove it, remove all connections, and invoke check_deletion on the connected clusters
            self._dead = True
            if self.parent is None:
                self.clusters[pos].remove(self)
            else:
                self.parent.children.remove(self)
                self.parent.check_deletion()
            # remove all our connection,
            # check if the thereby disconnected clusters are dead
            for o in list(self.back_connections):
                o.connections.remove(self)
                o.check_deletion()
            for o in list(self.connections):
                o.back_connections.remove(self)
                o.check_deletion()
            break
        self._checking_delete = False

    def reparent_children(self):
        for c in self.children:
            assert c.parent is self
            c.parent = self.parent
            if self.parent is None:
                frag0 = self.clustid[0]
                if frag0 > 1000:
                    pos = 2 * (frag0-1001) + 1
                else:
                    pos = 2 * (frag0-1)
                self.clusters[pos].add(c)

    def decompose(self, max_rmsd, fwd):
        '''
        When a node-form parent cluster is subdivided into children clusters,
        propagate from parent to children the connections with downstream(fwd=True) or upstream(fwd=False) clusters that are in overlap range
        then delete parent and connect children to the grand-parent (which is normally root)
        '''
        c1 = self.coors   # coordinates of the representatives of the children clusters at next level.
        if fwd:
            if not self.connections: return
            others = list(self.connections)
            # c2 = coord of the representative of the connected clusters of fwd frag
            c2 = np.concatenate([c.coors[0] for c in self.connections]).reshape(len(self.connections), len(self.connections[0].coors[0]), 3)
        else:
            if not self.back_connections: return
            others = list(self.back_connections)
            # c2 = coord of the representative of the connected clusters of bwd frag
            c2 = np.concatenate([c.coors[0] for c in self.back_connections]).reshape(len(self.back_connections), len(self.back_connections[0].coors[0]), 3)

        # For each of the representative of the connected clusters of fwd/bwd fragment (=c2),
        # check connectivity with each of the representative of the children clusters (= c1)
        # Divide c2 into chunks for memory saving
        chunksize = MAXCHUNK//len(c1)
        for chunkpos in range(0, len(c2), chunksize):
            # Account for clustering cutoff to check if any pose in cluster1 could
            # overlap with any pose in cluster2
            c_max_rmsd = [CLUSTERING[child.clusterlevel] * CLUST_MARGIN for child in self.children]
            max_rmsd0 = np.fromiter(c_max_rmsd,count=len(c_max_rmsd),dtype=float)[:, np.newaxis]
            max_rmsd2 = max_rmsd0 + max_rmsd # max_rmsd = cutoff given by user
            max_sd = (max_rmsd2**2) * c1.shape[1]

            c2_chunk = c2[chunkpos:chunkpos+chunksize]
            others_chunk = others[chunkpos:chunkpos+chunksize]
            if fwd:
                ocon = [o.back_connections for o in others_chunk]
                childcon = [c.connections for c in self.children]
            else:
                ocon = [o.connections for o in others_chunk]
                childcon = [c.back_connections for c in self.children]

            for o in ocon:
                o.remove(self)

            msd = get_msd(c1, c2_chunk)

            ch = self.children
            old_childnr = None

            for childnr, onr in zip(*np.where(msd < max_sd)):
                if childnr != old_childnr:
                    c_child = ch[childnr]
                    c_childcon = childcon[childnr]
                    old_childnr = childnr
                ocon[onr].append(c_child)
                c_childcon.append(others_chunk[onr])

            for o in others:
                o.check_deletion()

    def decompose_intra(self, fwd):
        if fwd:
            if not self.connections: return
            others = list(self.connections)
        else:
            if not self.back_connections: return
            others = list(self.back_connections)

        if fwd:
            ocon = [o.back_connections for o in others]
            childcon = [c.connections for c in self.children]
        else:
            ocon = [o.connections for o in others]
            childcon = [c.back_connections for c in self.children]
        for o in ocon:
            o.remove(self)

        for childnr, child in enumerate(self.children):
            cranks = child.all_ranks
            for onr, o in enumerate(others):
                if cranks.intersection(o.all_ranks):
                    ocon[onr].append(child)
                    childcon[childnr].append(o)

        for o in others:
            o.check_deletion()

    def verify(self, max_rmsd):
        ''' check that all kept connections have indeed low overlap rmsd'''
        if len(self.children): return
        if self.totstruc > 1: return
        cons = [con for con in self.connections if not len(con.children) and con.totstruc == 1]
        if not len(cons): return
        concoors = np.concatenate([con.coors[:1] for con in cons])
        c2 = self.coors[0]
        c1 = concoors
        d = c1 - c2
        msd = np.einsum("...ijk,ijk->...i", d,d)
        msd_low = (msd<(max_rmsd**2*c2.shape[0]))
        for n in range(len(cons)):
            assert msd_low[n]
        return

    def all_children(self):
        if not len(self.children):
            yield self
        else:
            for cc in self.children:
                for v in cc.all_children():
                    yield v

    def check_parentage(self, parent):
        parent_id = parent.clustid if parent is not None else None
        self_parent_id = self.parent.clustid if self.parent is not None else None
        assert self.parent == parent, (self.clustid, parent_id, self_parent_id)
        if self.parent is None:
            frag0 = self.clustid[0]
            if frag0 > 1000:
                pos = 2 * (frag0-1001) + 1
            else:
                pos = 2 * (frag0-1)
            assert self in self.clusters[pos], self.clustid
        for c in self.children:
            c.check_parentage(self)

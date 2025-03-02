from MDAnalysis.lib.distances import capped_distance,self_capped_distance
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
import MDAnalysis as mda
import numpy as np
import scipy.sparse
from scipy.spatial import KDTree, Voronoi
from scipy.linalg import eigh
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.distances import capped_distance
import warnings
import pandas as pd


class ClusterArea(AnalysisBase):
    def __init__(self, universe, residueGroup, selResidueGroup, cutoff=12.0):
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = [i for i in residueGroup]

        self.cutoff = cutoff

        self.headAtoms = self.u.atoms[[]]
        self.headSp = {}

        for i in range(len(self.residues)):
            sp = self.residues[i]
            self.headSp[sp] = ' '.join(residueGroup[sp])
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]))

        self._n_residues = self.headAtoms.n_residues
        self.resids = self.headAtoms.resindices

        self.selResiduesResid = self.u.select_atoms('resname %s' % ' '.join(selResidueGroup)).resindices
        self._residuesMask = np.in1d(self.headAtoms.residues.resindices, self.selResiduesResid)

    def _single_frame(self):
        pairs = capped_distance(
            self.headAtoms.positions,
            self.headAtoms.positions,
            max_cutoff=self.cutoff,
            box=self._ts.dimensions,
            return_distances=False,
        )
        ori, expan = np.unique(self.resids[pairs], axis=0).T

        dif = ori != expan
        ori = ori[dif]
        expan = expan[dif]

        one = np.ones_like(ori)
        neighborId = scipy.sparse.csr_matrix(
            (one, (ori, expan)),
            dtype=np.int8,
            shape=(self._n_residues, self._n_residues),
        )

        neighborIdFilter = neighborId[self._residuesMask][:, self._residuesMask]

        _, com_labels = scipy.sparse.csgraph.connected_components(neighborIdFilter)

        unique_com_labels, counts = np.unique(com_labels, return_counts=True)
        largest_label = unique_com_labels[np.argmax(counts)]
        # frame_resindices = self.headAtoms.residues.resindices[self._residuesMask]
        # largest_cluster_resindices = frame_resindices[com_labels == largest_label]
        area_arr = np.zeros([self._n_residues])
        centerPos = self.headAtoms.center_of_geometry(compound='residues')
        kdTree = KDTree(centerPos)
        for i, point in enumerate(centerPos):
            _, idxs = kdTree.query(point, 13)
            nearest_neighbors = np.vstack((centerPos[idxs[1:]], point))
            area_arr[i] = self.get_voronoi_area(nearest_neighbors)
        areaSel = np.sum(area_arr[self._residuesMask][com_labels == largest_label])
        print(areaSel)
    def get_voronoi_area(self, points):
        mean = np.mean(points, axis=0)
        centered_points = points - mean
        eigenvalues, eigenvectors = eigh(np.cov(centered_points.T))
        normal = eigenvectors[:, np.argmin(eigenvalues)]
        projected_p = points - np.dot(centered_points, normal)[:, np.newaxis] * normal
        reference_point = projected_p[-1]
        x_axis = (reference_point - mean) / np.linalg.norm(reference_point - mean)
        y_axis = np.cross(normal, x_axis)
        y_axis = y_axis / np.linalg.norm(y_axis)
        local_coordinates = projected_p - mean
        local_coordinates = np.column_stack((np.dot(local_coordinates, x_axis), np.dot(local_coordinates, y_axis)))
        vor = Voronoi(local_coordinates)
        region_index = vor.point_region[len(points) - 1]
        region_vertices = vor.regions[region_index]
        if -1 in region_vertices:
            return 0
        finite_vertices_indices = [index for index in region_vertices if index != -1]
        vertices = vor.vertices[finite_vertices_indices]
        return self.polygon_area(vertices)

    @staticmethod
    def polygon_area(vertices):
        n = vertices.shape[0]
        area = 0.0
        for i in range(n):
            x1, y1 = vertices[i]
            x2, y2 = vertices[(i+1) % n]
            area += x1 * y2 - y1 * x2
        return 0.5 * abs(area)


if __name__ == '__main__':
    import MDAnalysis as mda

    u = mda.Universe('E:/awork/lnb/june/01/B/ach.gro')
    cls1 = ClusterArea(u, {'DPPC': ['GL1', 'GL2'], 'D3PC': ['GL1', 'GL2'], 'CHOL': ['ROH']}, ['DPPC', 'CHOL'])
    cls1.run()
# u = mda.Universe("E:/awork/lnb/june/01/B/ach.gro")
# atm1pos = u.select_atoms('name PO4 ROH').positions
# atm1 = u.select_atoms('name PO4 ROH')
# heihei = u.select_atoms('resid 1 and name PO4')
# ls1 = AtomNeighborSearch(atm1)
# center = atm1.center_of_mass()
# print(np.linalg.norm(heihei.positions - center))
# atomresid = u.select_atoms('name PO4 ROH').resids
#
#
# #
# # z1 = np.zeros([atm1.n_residues, atm1.n_residues])
# # z1[ref, neigh] = 1
# # print(np.where(z1[:,]!=0))
#
# selResiudes = u.select_atoms('resname DPPC CHOL')
# residues = u.select_atoms('resname DPPC D3PC CHOL')
#
# centerPos = residues.center_of_mass(compound='residues')
#
# kdTree = KDTree(centerPos)
#
# area_arr = np.zeros([residues.n_residues])
# for i, point in enumerate(centerPos):
#     _, idxs = kdTree.query(point, 13)
#     nearest_neighbors = np.vstack((centerPos[idxs[1:]], point))
#     area_arr[i] = get_voronoi_area(nearest_neighbors)
#
# pairs = capped_distance(
#             atm1pos,
#             atm1pos,
#             max_cutoff=50,
#             box=u.trajectory.ts.dimensions,
#             return_distances=False,
#         )
#
#
# ref, neigh = np.unique(pairs, axis=0).T
#
# deleteSelf = ref != neigh
# ref = ref[deleteSelf]
# neigh = neigh[deleteSelf]
# one = np.zeros_like(ref)
#
# neighborId = scipy.sparse.csr_matrix(
#             (one, (ref, neigh)),
#             dtype=np.int8,
#             shape=(residues.n_residues, residues.n_residues),
#         )
# filter = np.in1d(residues.residues.resindices, selResiudes.residues.resindices)
# neighborFilter = neighborId[filter][:, filter]
# _, com_labels = scipy.sparse.csgraph.connected_components(neighborFilter)
# unique_com_labels, counts = np.unique(com_labels, return_counts=True)
# largest_label = unique_com_labels[np.argmax(counts)]
# frame_resindices = residues.residues.resindices[filter]
# largest_cluster_resindices = frame_resindices[com_labels == largest_label]
# largest_cluster = max(counts)
# print(largest_cluster)
# print(largest_cluster_resindices)
#
# areaSel = np.sum(area_arr[filter])
# print(areaSel)
# print(np.sum(area_arr))



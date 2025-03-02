"""
2024.6.17
未修改完成
1.无法修复点的位置在边界的情况
2.没有修改搜索k近邻问题

"""

import pandas as pd
try:
    from .analysis_base import AnalysisBase
except:
    from analysis_base import AnalysisBase
import MDAnalysis as mda
import numpy as np
from scipy.spatial import KDTree, Voronoi
from scipy.linalg import eigh
import warnings
import time
warnings.filterwarnings('ignore')
# suppress some MDAnalysis warnings when writing PDB/GRO files
__all__ = ['Area']


class Area(AnalysisBase):
    def __init__(self, universe, residueGroup: dict, k: int, file_path=None):
        super().__init__(universe.trajectory)
        self.u = universe
        self.k = k + 1
        self.filePath = file_path
        self.residues = [i for i in residueGroup]
        self.headSp = {
            sp: residueGroup[sp][0] for sp in self.residues
        }

        self.headAtoms = self.u.atoms[[]]
        for i in range(len(self.residues)):
            sp = self.residues[i]
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]), updating=False)

        self._n_residues = self.headAtoms.n_residues

        self.resnames = self.headAtoms.resnames
        self.resids = self.headAtoms.resids

        self.results.Area = None
        self.current = 0

    @property
    def Area(self):
        return self.results.Area

    def _prepare(self):
        self.results.Area = np.full([self._n_residues, self.n_frames],
                                    fill_value=np.NaN)

    def _single_frame(self):
        headPos = self.headAtoms.positions
        kdtree = KDTree(headPos)
        area_arr = np.zeros([self._n_residues])
        # for i, point in enumerate(headPos):
        _, idxs = kdtree.query(headPos, self.k)
        # nearest_neighbors = headPos[idxs]
        # area_arr[i] = get_voronoi_area(nearest_neighbors)
        # self.results.Area[:, self._frame_index] = area_arr / 100

    def _conclude(self):
        if self.filePath:
            try:
                from .analysis_base import WriteExcelLipids
            except:
                from analysis_base import WriteExcelLipids
            lipids_ratio = {sp: self.u.select_atoms(f'resname {sp}').n_residues for sp in self.residues}
            dict_parameter = {'step': self.step, 'n_frames': self.n_frames, 'resids': self.resids,
                              'resnames': self.resnames,
                              'positions': self.headAtoms.positions, 'results': self.results.Area,
                              'file_path': self.filePath, 'description': 'Area(nm^2)',
                              'value_divition': 1, 'lipids_type': lipids_ratio}
            WriteExcelLipids(**dict_parameter).run()


def get_voronoi_area(points, point_index=0):
    if not isinstance(points, np.ndarray):
        points = np.array(points)
    if points.shape[1] != 3:
        raise ValueError("输入点必须是三维的。")

    mean = np.mean(points, axis=0)
    centered_points = points - mean

    cov_matrix = np.cov(centered_points.T)

    _, eigenvectors = eigh(cov_matrix)

    x_axis = eigenvectors[:, 2]
    y_axis = eigenvectors[:, 1]

    local_x = np.dot(centered_points, x_axis)
    local_y = np.dot(centered_points, y_axis)
    local_coords_2d = np.column_stack((local_x, local_y))

    vor = Voronoi(local_coords_2d)

    region_index = vor.point_region[point_index]
    region_vertices = vor.regions[region_index]

    if -1 in region_vertices:
        return np.NaN
    finite_vertices = vor.vertices[region_vertices]
    return polygon_area(finite_vertices)


def polygon_area(vertices):
    x = vertices[:, 0]
    y = vertices[:, 1]
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))




if __name__ == "__main__":
    import time
    # u = mda.Universe("E:/awork/lnb/june/01/B/ach.gro", "E:/awork/lnb/june/01/B/lnb_nojump.xtc")
    # u = mda.Universe("E:/1130.gro")
    u = mda.Universe("E:/ach.gro")
    cls2 = Area(u, {'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}, 21, file_path='E:/excel/area_k.xlsx')
    print(cls2.residues)
    t1 = time.time()
    cls2.run()
    t2 = time.time()
    print('time', t2 - t1)


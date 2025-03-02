"""
2024.6.17
未修改完成
1.无法修复点的位置在边界的情况
2.没有修改搜索k近邻问题

"""
from itertools import combinations

import pandas as pd
from shapely import Polygon

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
        for i, point in enumerate(headPos):
            _, idxs = kdtree.query(point, self.k)
            nearest_neighbors = headPos[idxs]
            area_arr[i] = get_voronoi_area(nearest_neighbors)
        self.results.Area[:, self._frame_index] = area_arr * 0.01

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


def get_voronoi_area(points):
    # 输入验证
    if not isinstance(points, np.ndarray):
        points = np.array(points)

    if points.shape[1] != 3:
        raise ValueError("输入点必须是三维的。")

    centered_points = points - points[0]

    cov_matrix = np.cov(centered_points.T)

    _, eigenvectors = eigh(cov_matrix)
    # 选择前两个特征向量作为局部坐标系的x和y轴
    x_axis = eigenvectors[:, 2]
    y_axis = eigenvectors[:, 1]
    # 转换到局部二维坐标系
    local_x = np.dot(centered_points, x_axis)
    local_y = np.dot(centered_points, y_axis)
    local_coords_2d = np.column_stack((local_x, local_y))
    area = voronoi_area(local_coords_2d[0], local_coords_2d)
    return area
    # vor = Voronoi(local_coords_2d)
    # # 获取指定点的Voronoi区域
    # region_index = vor.point_region[0]
    # region_vertices = vor.regions[region_index]
    # # 检查区域是否为有限区域
    # if -1 in region_vertices:
    #     return np.NaN
    # finite_vertices = vor.vertices[region_vertices]
    # return polygon_area(finite_vertices)

def circumcenter(A, B, C):
    D = 2 * (A[0]*(B[1]-C[1]) + B[0]*(C[1]-A[1]) + C[0]*(A[1]-B[1]))
    if D == 0:
        return None  # Collinear points
    Px = ((A[0]**2 + A[1]**2)*(B[1]-C[1]) + (B[0]**2 + B[1]**2)*(C[1]-A[1]) + (C[0]**2 + C[1]**2)*(A[1]-B[1])) / D
    Py = ((A[0]**2 + A[1]**2)*(C[0]-B[0]) + (B[0]**2 + B[1]**2)*(A[0]-C[0]) + (C[0]**2 + C[1]**2)*(B[0]-A[0])) / D
    return (Px, Py)

def cross(o, a, b):
    return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])

def convex_hull(points):
    points = sorted(points, key=lambda p: (p[0], p[1]))
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]

def polygon_area(vertices):
    area = 0.0
    n = len(vertices)
    for i in range(n):
        j = (i + 1) % n
        area += vertices[i][0] * vertices[j][1] - vertices[j][0] * vertices[i][1]
    return abs(area) / 2.0

def voronoi_area(ref, nei):
    circumcenters = []
    for p1, p2 in combinations(nei, 2):
        cc = circumcenter(ref, p1, p2)
        if cc is not None:
            circumcenters.append(cc)
    if len(circumcenters) < 3:
        return 0.0  # Not enough points to form a polygon
    hull = convex_hull(circumcenters)
    if len(hull) < 3:
        return 0.0  # Not a polygon
    area = polygon_area(hull)
    return area
if __name__ == "__main__":
    import time
    # u = mda.Universe("E:/awork/lnb/june/01/B/ach.gro", "E:/awork/lnb/june/01/B/lnb_nojump.xtc")
    u = mda.Universe("E:/out_bicell_new.gro")
    # u = mda.Universe("E:/ach.gro")
    cls2 = Area(u, {'DPPC': ['PO4']}, 21, file_path='E:/excel/area_21_fuck.csv')
    print(cls2.residues)
    t1 = time.time()
    cls2.run()
    t2 = time.time()
    print('time', t2 - t1)

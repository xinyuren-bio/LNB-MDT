import warnings
warnings.filterwarnings('ignore')

import MDAnalysis as mda
import numpy as np
from scipy.spatial import KDTree, Voronoi
from scipy.linalg import eigh

try:
    from .analysis_base import *
except:
    from analysis_base import *
# suppress some MDAnalysis warnings when writing PDB/GRO files
__all__ = ['Area']


class Area(AnalysisBase):
    def __init__(self, universe, residueGroup: dict, k: int, file_path=None):
        super().__init__(universe.trajectory)
        self.u = universe
        self.k = k + 1
        self.filePath = file_path
        self.residues = list(residueGroup)
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
        n_points = len(headPos)
        k = self.k

        # ================== 向量化邻近点查询 ==================
        kd_tree = KDTree(headPos)
        _, all_idxs = kd_tree.query(headPos, k=k)  # 批量查询所有点

        # 构建扩展的邻近点数组 (n_points, k, 3)
        expanded_neighbors = headPos[all_idxs]

        # ================== 批量计算局部坐标系 ==================
        # 计算每个组的中心点 (n_points, 3)
        means = np.mean(expanded_neighbors, axis=1, keepdims=True)
        centered_points = expanded_neighbors - means  # (n_points, k, 3)

        # 计算协方差矩阵 (n_points, 3, 3)
        cov_matrices = np.einsum('nki,nkj->nij', centered_points, centered_points) / k

        # 批量特征分解
        _, eigenvectors = np.linalg.eigh(cov_matrices)

        # 提取局部坐标轴 (n_points, 3, 3)
        x_axes = eigenvectors[:, :, 2]  # 第一主方向
        y_axes = eigenvectors[:, :, 1]  # 第二主方向

        # 转换到局部二维坐标系 (n_points, k, 2)
        local_x = np.einsum('nkj,nj->nk', centered_points, x_axes)  # (n_points, k)
        local_y = np.einsum('nkj,nj->nk', centered_points, y_axes)  # (n_points, k)
        local_coords_2d = np.stack((local_x, local_y), axis=-1)  # (n_points, k, 2)

        # ================== 批量计算Voronoi面积 ==================
        def compute_voronoi_area(i):
            points = local_coords_2d[i]  # (k, 2)
            try:
                vor = Voronoi(points)
                region_index = vor.point_region[0]  # 第一个点是中心点
                region_vertices = vor.regions[region_index]
                if -1 in region_vertices:  # 检查是否为有限区域
                    return np.NaN
                finite_vertices = vor.vertices[region_vertices]
                return polygon_area(finite_vertices)
            except Exception:
                return np.NaN

        # 逐点计算Voronoi面积
        area_arr = np.array([compute_voronoi_area(i) for i in range(n_points)])

        # 替换NaN值为平均值
        mean_area = np.nanmean(area_arr)
        area_arr[np.isnan(area_arr)] = mean_area

        # 存储结果
        self.results.Area[:, self._frame_index] = area_arr * 0.01

    def _conclude(self):
        if self.filePath:
            lipids_ratio = {sp: self.u.select_atoms(f'resname {sp}').n_residues for sp in self.residues}
            dict_parameter = {'step': self.step, 'n_frames': self.n_frames, 'resids': self.resids,
                              'resnames': self.resnames,
                              'positions': self.headAtoms.positions, 'results': self.results.Area,
                              'file_path': self.filePath, 'description': 'Area(nm^2)',
                              'value_divition': 1, 'lipids_type': lipids_ratio}
            WriteExcelLipids(**dict_parameter).run()


def polygon_area(vertices):
    # 计算多边形面积（顺时针或逆时针）
    x = vertices[:, 0]
    y = vertices[:, 1]
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


if __name__ == "__main__":
    import time
    u = mda.Universe("E:/ach.gro", 'E:/ach.xtc')
    cls2 = Area(u, {'DPPC': ['PO4'], 'DAPC':['PO4'], 'CHOL':['ROH']}, 18, file_path='E:untitled2.csv')
    t1 = time.time()
    cls2.run(1, 100,verbose=True)
    t2 = time.time()
    print('time', t2 - t1)

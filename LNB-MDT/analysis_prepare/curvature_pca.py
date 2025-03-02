try:
    from .analysis_base import AnalysisBase
except:
    from analysis_base import AnalysisBase
import MDAnalysis as mda
import numpy as np
from numpy.linalg import lstsq, eigh
from scipy.spatial import KDTree
import warnings

warnings.filterwarnings('ignore')

__all__ = ['Curvature']


class Curvature(AnalysisBase):
    def __init__(self, u, residueGroup: dict, k: int, path=None, method='mean'):
        super().__init__(u.trajectory)
        # self.residues_dict = residueGroup
        self.u = u
        self.k = k
        self.filePath = path
        self.method = method

        self.residues = [i for i in residueGroup]

        self.headAtoms = self.u.atoms[[]]
        self.headSp = {}

        for i in range(len(self.residues)):
            sp = self.residues[i]

            self.headSp[sp] = residueGroup[sp][0]

            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]), updating=False)

        self._n_residues = self.headAtoms.n_residues
        self.resnames = self.headAtoms.resnames
        self.resids = self.headAtoms.resids
        self._headMask = {
            sp: self.resnames == sp for sp in self.residues
        }

        self.results.MeanCurvature = None
        self.results.GussainCurvature = None
        self.results.Normal = None

    @property
    def MeanCurvature(self):
        return self.results.MeanCurvature

    @property
    def GussainCurvature(self):
        return self.results.GussainCurvature

    def _prepare(self):
        self.results.MeanCurvature = np.full([self._n_residues, self.n_frames]
                                             , fill_value=np.NaN)

        self.results.GussainCurvature = np.full([self._n_residues, self.n_frames]
                                                , fill_value=np.NaN)

        self.results.Normal = np.full([self._n_residues, self.n_frames, 3],
                                      fill_value=np.NaN)

    def _single_frame(self):
        head_positions = self.headAtoms.positions
        kd_tree = KDTree(head_positions)
        for i, point in enumerate(head_positions):
            _, idxs = kd_tree.query(point, self.k)
            nearest_neighbors = head_positions[idxs]
            lambda_1, lambda_2 = get_eigenvalues(nearest_neighbors)
            if self.method == 'mean':
                k_mean = (lambda_1 + lambda_2) * 0.5
                self.results.MeanCurvature[i, self._frame_index] = k_mean

            elif self.method == 'gussain':
                k_gussain = lambda_2 * lambda_1
                self.results.GussainCurvature[i, self._frame_index] = k_gussain

    def _conclude(self):
        if self.filePath:
            try:
                from .analysis_base import WriteExcelLipids
            except:
                from analysis_base import WriteExcelLipids
            lipids_ratio = {sp: self.u.select_atoms(f'resname {sp}').n_residues for sp in self.residues}
            dict_parameter = {'step': self.step, 'n_frames': self.n_frames, 'resids': self.resids,
                              'resnames': self.resnames,
                              'positions': self.headAtoms.positions, 'results': self.results.MeanCurvature,
                              'file_path': self.filePath, 'description': 'Mean Curvature(nm -1)',
                              'value_divition': 1, 'lipids_type': lipids_ratio}
            WriteExcelLipids(**dict_parameter).run()


def get_eigenvalues(points):
    mean = np.mean(points, axis=0)
    centered_points = points - mean
    covariance_matrix = np.cov(centered_points.T)
    eigenvalues, _ = eigh(covariance_matrix)
    # print(eigenvalues)
    return eigenvalues[2], eigenvalues[1]


if __name__ == "__main__":
    # 导入结构文件、轨迹文件（可选）
    u = mda.Universe("E:/ach.gro")
    # 上传参数：
    #           1. 需要分析的残基及其相关原子
    #           2. K值
    #           3. 结果保存路径
    cls = Curvature(u, {'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}, 20, path='E:/excel/cur_k_pca.xlsx')
    # 执行计算
    cls.run()


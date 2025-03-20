import warnings
warnings.filterwarnings('ignore')

import pandas as pd
from scipy.linalg import eigh
import numpy as np

try:
    from .analysis_base import *
except:
    from analysis_base import *

__all__ = ['Anisotropy']


class Anisotropy(AnalysisBase):

    def __init__(self, u, residuesGroup, filePath):
        super().__init__(u.trajectory)
        self.u = u
        self.residues = list(residuesGroup)
        self.filePath = filePath

        self.headSp = {
            sp: residuesGroup[sp][0] for sp in self.residues
        }
        print(self.headSp)
        self.headAtoms = self.u.atoms[[]]

        for i in range(len(self.residues)):
            sp = self.residues[i]
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]), updating=False)

        self._n_residues = self.headAtoms.n_residues

        self.parameters = str(residuesGroup)

    @property
    def Anisotropy(self):
        return self.results.Anisotropy

    def _prepare(self):
        self.results.Anisotropy = np.full([self.n_frames], fill_value=np.NaN)

    def _single_frame(self):
        atomsPos = self.headAtoms.positions
        ovariance_matrix = np.cov(atomsPos.T)
        # 计算特征值
        eigenvalues, _ = eigh(ovariance_matrix)
        # print(q, eigenvalues)
        # 计算椭圆偏心率
        anisotropy = np.sqrt(1 - eigenvalues[0] ** 2 / eigenvalues[2] ** 2)
        self.results.Anisotropy[self._frame_index] = anisotropy

    def _conclude(self):
        if self.filePath:
            dict_parameter = {'step': self.step, 'n_frames': self.n_frames,
                             'results': self.results.Anisotropy,
                              'file_path': self.filePath, 'description': 'Anisotropy',
                              'parameters': self.parameters}
            WriteExcelBubble(**dict_parameter).run()


if __name__ == "__main__":

    import MDAnalysis as mda
    u = mda.Universe('E:/ach.gro', 'E:/ach.xtc')
    cls1 = Anisotropy(u, {'DPPC':['PO4'], 'DAPC':['PO4'], 'CHOL':['ROH']}, filePath='E:/untitled1.csv')
    cls1.run(0, 100)


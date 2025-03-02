import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
from scipy.linalg import eig

try:
    from .analysis_base import *
except:
    from analysis_base import *

__all__ = ['PCA']


class PCA(AnalysisBase):
    def __init__(self, u, residuesGroup, filePath):
        super().__init__(u.trajectory)
        self.u = u
        self.residues = list(residuesGroup)
        self.filePath = filePath

        self.headSp = {
            sp: residuesGroup[sp][0] for sp in self.residues
        }
        self.headAtoms = self.u.atoms[[]]

        for i in range(len(self.residues)):
            sp = self.residues[i]
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]), updating=False)

        self.n_residues = self.headAtoms.n_residues

    def _prepare(self):
        self.results.PCA = np.full([self.n_frames], fill_value=np.NaN)

    def _single_frame(self):
        head_pos = self.headAtoms.positions
        center = self.headAtoms.center_of_mass()
        head2center = head_pos - center
        covMatix = np.cov(head2center, rowvar=False)
        eigenValues, _ = eig(covMatix)
        sorted_indices = np.argsort(eigenValues)[::-1]
        eigenvalues = eigenValues[sorted_indices]
        varianceRatio = eigenvalues[-1] / eigenvalues[0]
        self.results.PCA[self._frame_index] = varianceRatio

    def _conclude(self):
        if self.filePath:
            dict_parameter = {'step': self.step, 'n_frames': self.n_frames,
                              'results': self.results.PCA,
                              'file_path': self.filePath, 'description': 'PCA',
                              'value_divition': 1}
            WriteExcelBubble(**dict_parameter).run()


if __name__ == '__main__':
    import MDAnalysis as mda
    u = mda.Universe("E:/ach.gro", "E:/ach.xtc")
    cls2 = PCA(u, ['DPPC','D3PC','CHOL'], path='E:/excel')
    cls2.run(1,100,1,verbose=True)
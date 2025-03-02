import warnings
warnings.filterwarnings('ignore')

import MDAnalysis as mda
import numpy as np
from scipy.spatial import KDTree
from scipy.linalg import eigh

try:
    from .analysis_base import *
except:
    from analysis_base import *


__all__ = ['Height']


class Height(AnalysisBase):

    def __init__(self, universe, residuesGroup: dict, k: int = None, file_path: str = None):
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residuesGroup)
        self.k = k
        self.file_path = file_path

        self.headSp = {sp: ' '.join(residuesGroup[sp][0]) for sp in residuesGroup}
        self.tailSp = {sp: ' '.join(residuesGroup[sp][-1]) for sp in residuesGroup}

        self.headAtoms = self.u.atoms[[]]
        self.tailAtoms = self.u.atoms[[]]

        for i in range(len(self.residues)):
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (self.residues[i], self.headSp[self.residues[i]]), updating=False)

            self.tailAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (self.residues[i], self.tailSp[self.residues[i]]), updating=False)

        self._n_residues = self.headAtoms.n_residues
        self.resids = self.headAtoms.resids
        self.resnames = self.headAtoms.resnames
        self.resArrange = np.argsort(np.argsort(self.headAtoms.resindices))
        self.results.Height = None

    def _prepare(self):
        self.results.Height = np.full([self._n_residues, self.n_frames],
                                      fill_value=np.NaN)

    def _single_frame(self):
        centerHeadSp = self.headAtoms.positions
        normals = get_normals(self.k, centerHeadSp)
        centerTailSp = self.tailAtoms.center_of_geometry(compound='residues')
        head_to_tail = centerHeadSp - centerTailSp
        distance = (np.abs(np.einsum('ij,ij->i', normals, head_to_tail)))[self.resArrange]
        self.results.Height[:, self._frame_index] = distance * 0.1  # A° to nm

    @property
    def Height(self):
        return self.results.Height

    def _conclude(self):
        if self.file_path:
            lipids_ratio = {sp: self.u.select_atoms(f'resname {sp}').n_residues for sp in self.residues}
            dict_parameter = {'step':self.step, 'n_frames': self.n_frames, 'resids':self.resids, 'resnames':self.resnames,
             'positions':self.headAtoms.positions, 'results':self.results.Height, 'file_path':self.file_path,'description':'Height(nm)' ,
                              'value_divition':1, 'lipids_type':lipids_ratio}
            WriteExcelLipids(**dict_parameter).run()


if __name__ == "__main__":
    import time
    t1 = time.time()
    u = mda.Universe("E:/ach.gro", 'E:/ach.xtc')
    dict_residue = {'DPPC': (['PO4'], ['C4B', 'C4A']), 'DAPC':(['PO4'], ['C5A', 'C5B']), 'CHOL':(['ROH'], ['R5'])}
    cls2 = Height(u, dict_residue, k=21, file_path='E:/untitled4.csv')
    cls2.run(1, 100)
    t2 = time.time()
    print(t2 - t1)
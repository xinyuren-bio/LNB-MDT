import warnings
warnings.filterwarnings('ignore')

import numpy as np

try:
    from .analysis_base import *
except:
    from analysis_base import *

__all__ = ['Gyration_py']


class Gyration_py(AnalysisBase):
    def __init__(self, u, residuesGroup, filePath):
        super().__init__(u.trajectory)
        self.u = u
        self.residues = list(residuesGroup)
        self.filePath = filePath
        print(self.residues)
        self.headSp = {
            sp: ' '.join(residuesGroup[sp]) for sp in self.residues
        }
        print(self.headSp)
        self.headAtoms = self.u.atoms[[]]

        for i in range(len(self.residues)):
            sp = self.residues[i]
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]), updating=False)

        self._n_residues = self.headAtoms.n_residues
        self.results.Gyration = None

    def _prepare(self):
        self.results.Gyration = np.zeros([self.n_frames])

    def _single_frame(self):
        atomsPos = self.headAtoms.center_of_mass(compound='residues') / 10
        center = self.headAtoms.center_of_mass() / 10
        repos = atomsPos - center
        masses = self.headAtoms.residues.masses

        self.results.Gyration[self._frame_index] = get_radius_of_gyration(repos, masses)

    def _conclude(self):
        if self.filePath:
            try:
                from .analysis_base import WriteExcelBubble
            except:
                from analysis_base import WriteExcelBubble
            dict_parameter = {'step': self.step, 'n_frames': self.n_frames,
                              'results': self.results.Gyration,
                              'file_path': self.filePath, 'description': 'Gyration(nm)',
                              'value_divition': 1}
            WriteExcelBubble(**dict_parameter).run()


def get_radius_of_gyration(positions, masses):
    total_masses = np.sum(masses)
    n = positions.shape[0]
    view_norm = np.zeros([n])
    total = 0
    for i in range(n):
        for j in range(3):
            view_norm[i] += (positions[i, j] * positions[i, j])
        total += view_norm[i] * masses[i]
    heihei = total/ total_masses
    return heihei ** 0.5


if __name__ == "__main__":

    import MDAnalysis as mda
    u = mda.Universe('E:/ach.gro', 'E:/ach.xtc')
    cls1 = Gyration_py(u, {'DPPC':['PO4'], 'DAPC':['PO4'], 'CHOL':['ROH']}, filePath='E:/untitled1.csv')
    cls1.run(0, 100)


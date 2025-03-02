import logging
import warnings
from collections import OrderedDict

import numpy as np
import pandas as pd
from ..analysis.analysis_base import AnalysisBase

warnings.filterwarnings('ignore')


class RDF(AnalysisBase):
    def __init__(self
                 , u
                 , center_residues: tuple
                 , rdf_residues: dict[str, list]
                 , n_circle: int = 50
                 , ranges: tuple = None
                 , save_path: str = None) -> None:

        super().__init__(u.trajectory)
        self.u = u
        self.n_circle: int = n_circle

        if ranges is None or ranges[0] >= ranges[1] or len(ranges) != 2:
            self.ranges = np.min(self.u.dimensions[:3] / 2)
        else:
            self.ranges = ranges

        self.save_path = save_path

        self.center_residues: str = ' '.join(center_residues)
        self.rdf_residues: list[str] = [i for i in rdf_residues]

        self.rdf_atoms = self.u.atoms[[]]
        self._num_residues = [0]

        for i, sp in enumerate(self.rdf_residues):
            self.rdf_atoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, ' '.join(rdf_residues[sp])))
            self._num_residues.append(self.rdf_atoms.n_residues)

        self._sum_volume: float = 0

        self.results.RDF_all_frame = None
        self.results.RDF = None

    @property
    def RDF(self):
        return self.results.RDF

    def _prepare(self):
        self.results.RDF_all_frame = np.zeros([len(self.rdf_residues), self.n_circle])
        self.results.RDF = np.zeros([len(self.rdf_residues), self.n_circle])

    def _single_frame(self):
        self._sum_volume += self._ts.volume

        center_position = self.u.select_atoms(self.center_residues).center_of_mass()

        residues_positions = self.rdf_atoms.center_of_mass(compound='residues')

        distances_to_center = np.linalg.norm(residues_positions - center_position, axis=1)
        for i in range(len(self.rdf_residues)):
            count, _ = np.histogram(distances_to_center[self._num_residues[i]:self._num_residues[i+1]]
                                    , bins=self.n_circle
                                    , range=self.ranges)

            self.results.RDF_all_frame[i] += count

    def _conclude(self):
        _, edges = np.histogram([-1], self.n_circle, self.ranges)
        ball_volume = 4 / 3 * np.power(edges, 3) * np.pi
        shell_volume = np.diff(ball_volume)
        for i in range(len(self.rdf_residues)):
            self.results.RDF[self._num_residues[i]: self._num_residues[i+1]] \
                = (self.results.RDF_all_frame[self._num_residues[i]: self._num_residues[i+1]]
                   / shell_volume / (int(self._num_residues[i + 1] - self._num_residues[i]) / self._sum_volume * self.n_frames))

    # def writeExcel(self):
    #     arr1 = self.r.reshape([-1, 1]).T
    #     combined_array = np.row_stack((arr1, self.results.Rad))
    #     columnHead = ['Distance'] + [sp for sp in self.residues]
    #     df = pd.DataFrame(combined_array.T, columns=columnHead)
    #     explanation_df = pd.DataFrame(['Radial Distribution'])
    #     with pd.ExcelWriter(self.filePath, engine='openpyxl') as writer:
    #         explanation_df.to_excel(writer, sheet_name='Sheet1', index=False, header=False)
    #         df.to_excel(writer, sheet_name='Sheet1', index=False, startrow=1, header=True)


if __name__ == "__main__":
    import MDAnalysis as mda

    u = mda.Universe("C:/Users/59563/1.gro")
    # u = mda.Universe("E:/ach.gro")
    cls = RDF(u, {'N2': ['N2'], 'DPPC': ['NC3'], 'DAPC': ['NC3'], 'CHOL': ['ROH']}, ncircle=50,
                 filePath='E:/excel/rad.xlsx')
    # cls = CalRad(u, {'DPPC':['NC3'],'D3PC':['NC3'],'CHOL':['ROH']}, ncircle=50, filePath='E:/excel/rad1.xlsx')
    cls.run()

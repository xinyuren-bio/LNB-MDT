import warnings
import os
import sys
warnings.filterwarnings('ignore')

import pandas as pd
from scipy.linalg import eigh
import numpy as np
import MDAnalysis as mda
from joblib import Parallel, delayed

if __name__ == '__main__':
    current_file_path = os.path.abspath(__file__)
    current_dir = os.path.dirname(current_file_path)
    package_root = os.path.abspath(os.path.join(current_dir, '..'))
    if package_root not in sys.path:
        sys.path.insert(0, package_root)
from analysis.analysis_base import *

__all__ = ['Anisotropy']


class Anisotropy(AnalysisBase):
    def __init__(self, universe: mda.Universe, residueGroup: dict, file_path: str = None, parallel: bool = False, n_jobs: int = -1):
        super().__init__(universe.trajectory)
        self.u = universe
        self.file_path = file_path
        self.parallel = parallel
        self.n_jobs = n_jobs
        
        self.residues = list(residueGroup)
        self.headSp = {sp: ' '.join(residueGroup[sp]) for sp in residueGroup}
        print("Head atoms:", self.headSp)

        self.headAtoms = self.u.atoms[[]]

        # 为了保持原子顺序，循环选择指定脂质类型的头部原子并添加到 AtomGroup
        for i in range(len(self.residues)):
            self.headAtoms += self.u.select_atoms(
                f'resname {self.residues[i]} and name {self.headSp[self.residues[i]]}', 
                updating=False
            )
        
        if self.headAtoms.n_atoms == 0:
            raise ValueError("Atom selection is empty. Please check your `residueGroup` dictionary and atomic names.")

        self._n_residues = self.headAtoms.n_residues
        self.resids = self.headAtoms.resids
        self.resnames = self.headAtoms.resnames
        self.results.Anisotropy = None

        self.parameters = str(residueGroup)

    @property
    def Anisotropy(self):
        return self.results.Anisotropy

    def _prepare(self):
        self.results.Anisotropy = np.full(self.n_frames, fill_value=np.nan)

    @staticmethod
    def _calculate_anisotropy(positions: np.ndarray) -> float: 
        covariance_matrix = np.cov(positions.T)
        eigenvalues, _ = eigh(covariance_matrix)
        anisotropy = np.sqrt(1 - eigenvalues[0]**2 / eigenvalues[2]**2)
        return anisotropy

    def _single_frame(self):
        positions = self.headAtoms.positions
        self.results.Anisotropy[self._frame_index] = self._calculate_anisotropy(positions)

    def _get_positions_generator(self):
        for ts in self.u.trajectory[self.start:self.stop:self.step]:
            yield self.headAtoms.positions

    def run(self, start=None, stop=None, step=None, verbose=None):
        if self.parallel:
            self.start = start if start is not None else 0
            self.stop = stop if stop is not None and stop < self._trajectory.n_frames else self._trajectory.n_frames
            self.step = step if step is not None else 1
            self.n_frames = len(range(self.start, self.stop, self.step))
            self._prepare()
            print(f"Running in parallel on {self.n_jobs} jobs...")
            positions_generator = self._get_positions_generator()
            verbose_level = 10 if verbose else 0
            results_list = Parallel(n_jobs=self.n_jobs, verbose=verbose_level)(
                delayed(self._calculate_anisotropy)(pos) for pos in positions_generator
            )
            self.results.Anisotropy = np.array(results_list)
            self._conclude()
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose)
    def _conclude(self):
        if self.file_path:
            dict_parameter = {
                'frames': list(range(self.start, self.stop, self.step)),
                'results': self.results.Anisotropy,
                'file_path': self.file_path,
                'description': 'Anisotropy',
                'parameters': self.parameters
            }
            WriteExcelBubble(**dict_parameter).run()
            print(f"Results saved to {self.file_path}")


if __name__ == "__main__":
    gro_file = "cases/lnb.gro"
    xtc_file = "cases/md.xtc"
    csv_file_serial = "cases/csv/anisotropy_serial.csv"
    csv_file_parallel = "cases/csv/anisotropy_parallel.csv"

    u = mda.Universe(gro_file, xtc_file)
    
    residue_group = {'DPPC':['PO4'], 'DUPC':['PO4'], 'CHOL':['ROH']}

    # 1. 串行
    analysis_serial = Anisotropy(u, residue_group, file_path=csv_file_serial, parallel=False)
    analysis_serial.run(start=0, stop=100, verbose=True)

    # 2. 并行
    analysis_parallel = Anisotropy(u, residue_group, file_path=csv_file_parallel, parallel=True, n_jobs=2)
    analysis_parallel.run(start=0, stop=500 , verbose=True)


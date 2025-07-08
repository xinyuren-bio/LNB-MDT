import warnings
import os
import sys
warnings.filterwarnings('ignore')

import numpy as np
import MDAnalysis as mda
from scipy.linalg import eig
from joblib import Parallel, delayed

if __name__ == '__main__':
    current_file_path = os.path.abspath(__file__)
    current_dir = os.path.dirname(current_file_path)
    package_root = os.path.abspath(os.path.join(current_dir, '..'))
    if package_root not in sys.path:
        sys.path.insert(0, package_root)
from analysis.analysis_base import * 

__all__ = ['PCA']


class PCA(AnalysisBase):
    """
    A class for performing Principal Component Analysis on a group of atoms.
    """
    def __init__(self, universe, residues_group: dict, file_path: str = None, 
                 parallel: bool = False, n_jobs: int = -1):
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residues_group)
        self.file_path = file_path
        self.parallel = parallel
        self.n_jobs = n_jobs

        self.atomSp = {sp: ' '.join(residues_group[sp]) for sp in self.residues}
        print("Selected Atoms:", self.atomSp)

        self.headAtoms = self.u.atoms[[]]
        for sp in self.residues:
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.atomSp[sp]), updating=False)

        self.n_atoms = self.headAtoms.n_atoms
        self.results.PCA = None
        self.parameters = str(residues_group)

    def _prepare(self):
        self.results.PCA = np.full(self.n_frames, fill_value=np.nan)

    @staticmethod
    def _calculate_pca_ratio_for_frame(positions):
        if positions.shape[0] < 3:
            return np.nan
        
        center = np.mean(positions, axis=0)
        centered_pos = positions - center
        cov_matrix = np.cov(centered_pos, rowvar=False)
        
        eigenvalues, _ = eig(cov_matrix)
        eigenvalues = np.real(eigenvalues)
        eigenvalues.sort()

        if eigenvalues[2] == 0:
            return np.nan
            
        variance_ratio = eigenvalues[0] / eigenvalues[2]
        return variance_ratio

    def _single_frame(self):
        positions = self.headAtoms.positions
        ratio = PCA._calculate_pca_ratio_for_frame(positions)
        self.results.PCA[self._frame_index] = ratio

    def _get_inputs_generator(self):
        for _ in self.u.trajectory[self.start:self.stop:self.step]:
            yield (self.headAtoms.positions.copy(),)

    def run(self, start=None, stop=None, step=None, verbose=None):
        self.start = start if start is not None else 0
        self.stop = stop if stop is not None and stop < self._trajectory.n_frames else self._trajectory.n_frames
        self.step = step if step is not None else 1
        self.n_frames = len(range(self.start, self.stop, self.step))
        self._prepare()

        if self.parallel:
            print(f"Running in parallel on {self.n_jobs} jobs...")
            verbose_level = 10 if verbose else 0
            inputs_generator = self._get_inputs_generator()
            results_list = Parallel(n_jobs=self.n_jobs, verbose=verbose_level)(
                delayed(PCA._calculate_pca_ratio_for_frame)(*inputs) for inputs in inputs_generator
            )
            if results_list:
                self.results.PCA = np.array(results_list)
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose)

        self._conclude()

    def _conclude(self):
        if self.file_path:
            dict_parameter = {
                'frames': list(range(self.start, self.stop, self.step)),
                'results': self.results.PCA,
                'file_path': self.file_path,
                'description': 'PCA Variance Ratio (min/max)',
                'parameters': self.parameters
            }
            # Assuming WriteExcelBubble is defined and available
            WriteExcelBubble(**dict_parameter).run()
            print(f"Analysis complete. Results will be saved to {self.file_path}")


if __name__ == '__main__':
    gro_file = "cases/lnb.gro"
    xtc_file = "cases/md.xtc"
    csv_file_serial = "cases/csv/pca_serial.csv"
    csv_file_parallel = "cases/csv/pca_parallel.csv"

    u = mda.Universe(gro_file, xtc_file)
    residues_group = {'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}

    #1. 串行
    print("--- Running Serial Analysis ---")
    analysis_serial = PCA(u, residues_group, file_path=csv_file_serial, parallel=False)
    analysis_serial.run(verbose=True)

    #2. 并行
    print("--- Running Parallel Analysis ---")
    analysis_parallel = PCA(u, residues_group, file_path=csv_file_parallel, parallel=True, n_jobs=-1)
    analysis_parallel.run(verbose=True)
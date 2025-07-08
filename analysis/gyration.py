import warnings
import os
import sys
warnings.filterwarnings('ignore')

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

__all__ = ['Gyration']


def get_radius_of_gyration(positions, masses):
    total_mass = np.sum(masses)
    if total_mass == 0:
        return np.nan
    
    squared_distances = np.sum(positions**2, axis=1)
    mass_weighted_sq_dist = np.dot(squared_distances, masses)
    
    return np.sqrt(mass_weighted_sq_dist / total_mass)


class Gyration(AnalysisBase):
    """
    A class for calculating the radius of gyration for a selection of atoms.
    """
    def __init__(self, universe, residues_group: dict, file_path: str = None, 
                 parallel: bool = False, n_jobs: int = -1):
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residues_group)
        self.file_path = file_path
        self.parallel = parallel
        self.n_jobs = n_jobs

        self.headSp = {sp: ' '.join(residues_group[sp]) for sp in self.residues}
        print("Selected Atoms:", self.headSp)
        
        self.headAtoms = self.u.atoms[[]]
        for sp in self.residues:
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]), updating=False)
        
        self._n_residues = self.headAtoms.n_residues
        self.results.Gyration = None
        self.parameters = str(residues_group)

    def _prepare(self):
        self.results.Gyration = np.full(self.n_frames, fill_value=np.nan)

    @staticmethod
    def _calculate_gyration_for_frame(residue_com, group_com, residue_masses):
        relative_positions = residue_com - group_com
        return get_radius_of_gyration(relative_positions, residue_masses)

    def _single_frame(self):
        residue_com = self.headAtoms.center_of_mass(compound='residues') / 10
        group_com = self.headAtoms.center_of_mass() / 10
        residue_masses = self.headAtoms.residues.masses

        gyration_radius = Gyration._calculate_gyration_for_frame(residue_com, group_com, residue_masses)
        self.results.Gyration[self._frame_index] = gyration_radius

    def _get_inputs_generator(self):
        residue_masses = self.headAtoms.residues.masses
        for _ in self.u.trajectory[self.start:self.stop:self.step]:
            residue_com = self.headAtoms.center_of_mass(compound='residues') / 10
            group_com = self.headAtoms.center_of_mass() / 10
            yield (residue_com, group_com, residue_masses)

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
                delayed(Gyration._calculate_gyration_for_frame)(*inputs) for inputs in inputs_generator
            )
            if results_list:
                self.results.Gyration = np.array(results_list)
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose)

        self._conclude()

    def _conclude(self):
        if self.file_path:
            try:
                from .analysis_base import WriteExcelBubble
            except ImportError:
                from analysis_base import WriteExcelBubble
            
            dict_parameter = {
                'frames': list(range(self.start, self.stop, self.step)),
                'results': self.results.Gyration,
                'file_path': self.file_path,
                'description': 'Gyration (nm)',
                'parameters': self.parameters
            }
            WriteExcelBubble(**dict_parameter).run()
            print(f"Analysis complete. Results will be saved to {self.file_path}")


if __name__ == "__main__":
    import time
    gro_file = "cases/lnb.gro"
    xtc_file = "cases/md.xtc"
    csv_file_serial = "cases/csv/gyration_serial.csv"
    csv_file_parallel = "cases/csv/gyration_parallel.csv"

    u = mda.Universe(gro_file, xtc_file)
    residues_group = {'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}

    # 1. 串行
    print("--- Running Serial Analysis ---")
    t1 = time.time()
    analysis_serial = Gyration(u, residues_group, file_path=csv_file_serial, parallel=False)
    analysis_serial.run(verbose=True)
    t2 = time.time()
    print(f"Serial execution time: {t2 - t1:.2f} seconds")
    
    # 2. 并行
    print("--- Running Parallel Analysis ---")
    t1 = time.time()
    analysis_parallel = Gyration(u, residues_group, file_path=csv_file_parallel, parallel=True, n_jobs=-1)
    analysis_parallel.run()
    t2 = time.time()
    print(f"Parallel execution time: {t2 - t1:.2f} seconds")
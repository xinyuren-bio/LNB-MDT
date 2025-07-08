import warnings
import os
import sys
warnings.filterwarnings('ignore')

import MDAnalysis as mda
import numpy as np
from scipy.spatial import KDTree
from scipy.linalg import eigh
from joblib import Parallel, delayed

if __name__ == '__main__':
    current_file_path = os.path.abspath(__file__)
    current_dir = os.path.dirname(current_file_path)
    package_root = os.path.abspath(os.path.join(current_dir, '..'))
    if package_root not in sys.path:
        sys.path.insert(0, package_root)
from analysis.analysis_base import *

__all__ = ['Height']

def get_normals(k, points):
    n_points = points.shape[0]
    if n_points < k:
        return np.full_like(points, np.nan)

    kdtree = KDTree(points)
    _, idxs = kdtree.query(points, k=k)
    neighborhoods = points[idxs]

    means = np.mean(neighborhoods, axis=1, keepdims=True)
    centered = neighborhoods - means
    cov_matrices = np.einsum('nki,nkj->nij', centered, centered)

    eigenvalues, eigenvectors = eigh(cov_matrices)
    normals = eigenvectors[:, :, 0]

    geo_center = np.mean(points, axis=0)
    center_vectors = geo_center - points
    dots = np.einsum('ij,ij->i', normals, center_vectors)
    
    flip_mask = dots < 0
    normals[flip_mask] *= -1

    return normals


class Height(AnalysisBase):
    """
    A class for calculating the height of lipid molecules in LNB system.
    """

    def __init__(self, universe, residuesGroup: dict, k: int = 20, file_path: str = None, 
                 parallel: bool = False, n_jobs: int = -1):
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residuesGroup)
        self.k = k
        self.file_path = file_path
        self.parallel = parallel
        self.n_jobs = n_jobs

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

        self.parameters = str(residuesGroup) + 'K:' + str(self.k)

    def _prepare(self):
        self.results.Height = np.full([self._n_residues, self.n_frames],
                                      fill_value=np.nan)

    @staticmethod
    def _calculate_height_for_frame(head_pos, tail_pos, k, res_arrange):
        if head_pos.shape[0] == 0:
            return np.array([])
        normals = get_normals(k, head_pos)
        head_to_tail = head_pos - tail_pos
        distance = (np.abs(np.einsum('ij,ij->i', normals, head_to_tail)))[res_arrange]
        return distance * 0.1

    def _single_frame(self):
        head_pos = self.headAtoms.positions
        tail_pos = self.tailAtoms.center_of_geometry(compound='residues')
        heights = Height._calculate_height_for_frame(head_pos, tail_pos, self.k, self.resArrange)
        self.results.Height[:, self._frame_index] = heights

    @property
    def Height(self):
        return self.results.Height

    def _get_inputs_generator(self):
        for _ in self.u.trajectory[self.start:self.stop:self.step]:
            head_pos = self.headAtoms.positions.copy()
            tail_pos = self.tailAtoms.center_of_geometry(compound='residues')
            yield (head_pos, tail_pos, self.k, self.resArrange)

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
                delayed(Height._calculate_height_for_frame)(*inputs) for inputs in inputs_generator
            )
            if results_list:
                self.results.Height = np.array(results_list).T
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose)

        self._conclude()

    def _conclude(self):
        if self.file_path:
            lipids_ratio = {sp: self.u.select_atoms(f'resname {sp}').n_residues for sp in self.residues}
            dict_parameter = {
                'frames': list(range(self.start, self.stop, self.step)),
                'resids': self.resids,
                'resnames': self.resnames,
                'positions': self.headAtoms.positions,
                'results': self.results.Height,
                'file_path': self.file_path,
                'description': 'Height (nm)',
                'parameters': self.parameters,
                'lipids_type': lipids_ratio
            }
            # Assuming WriteExcelLipids is defined and available
            WriteExcelLipids(**dict_parameter).run()
            print(f"Analysis complete. Results saved to {self.file_path}")
        else:
            return self.results.Height


if __name__ == "__main__":
    import time
    gro_file = "cases/lnb.gro"
    xtc_file = "cases/md.xtc"
    csv_file_serial = "cases/csv/height_serial.csv"
    csv_file_parallel = "cases/csv/height_parallel.csv"

    if not (os.path.exists(gro_file) and os.path.exists(xtc_file)):
         print(f"Error: Input files not found. Please ensure '{gro_file}' and '{xtc_file}' exist.")
    else:
        u = mda.Universe(gro_file, xtc_file)
        residues_group = {'DPPC': (['PO4'], ['C4B', 'C4A']), 'DUPC':(['PO4'], ['C3A', 'C4B']), 'CHOL':(['ROH'], ['R5'])}
        
        print("--- Running Serial Analysis ---")
        t1 = time.time()
        analysis_serial = Height(u, residues_group, k=21, file_path=csv_file_serial, parallel=False)
        analysis_serial.run(verbose=True)
        t2 = time.time()
        print(f"Serial execution time: {t2 - t1:.2f} seconds")

        print("\n" + "="*50 + "\n")

        print("--- Running Parallel Analysis ---")
        t1 = time.time()
        analysis_parallel = Height(u, residues_group, k=21, file_path=csv_file_parallel, parallel=True, n_jobs=-1)
        analysis_parallel.run(verbose=True)
        t2 = time.time()
        print(f"Parallel execution time: {t2 - t1:.2f} seconds")
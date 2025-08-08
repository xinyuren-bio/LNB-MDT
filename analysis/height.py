import warnings
import os
import sys
import argparse  # 添加argparse导入
import ast       # 用于安全地解析residues_group字典字符串
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


# --- Command-line Argument Parsing ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform Height analysis on molecular dynamics trajectories."
    )

    parser.add_argument(
        "--gro-file",
        type=str,
        default="cases/lnb.gro",
        help="Path to the GRO file (topology file)."
    )
    parser.add_argument(
        "--xtc-file",
        type=str,
        default="cases/md.xtc",
        help="Path to the XTC file (trajectory file)."
    )
    parser.add_argument(
        "--output-csv",
        type=str,
        default="cases/csv/height_results.csv",
        help="Path to the output CSV file for height results."
    )
    parser.add_argument(
        "--residues",
        type=str,
        default="{'DPPC': (['PO4'], ['C4B', 'C4A']), 'DUPC':(['PO4'], ['C3A', 'C4B']), 'CHOL':(['ROH'], ['R5'])}",
        help="A dictionary string defining residue groups for analysis. E.g., \"{'DPPC': (['PO4'], ['C4B', 'C4A']), 'CHOL':(['ROH'], ['R5'])}\""
    )
    parser.add_argument(
        "--k-value",
        type=int,
        default=20,
        help="K value for height calculation."
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Enable parallel processing for height calculation."
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=-1,
        help="Number of jobs to run in parallel. -1 means using all available CPU cores."
    )
    parser.add_argument(
        "--start-frame",
        type=int,
        default=0,
        help="Starting frame for analysis (0-indexed)."
    )
    parser.add_argument(
        "--stop-frame",
        type=int,
        help="Stopping frame for analysis (exclusive). Defaults to end of trajectory."
    )
    parser.add_argument(
        "--step-frame",
        type=int,
        default=1,
        help="Step size for frames during analysis."
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output during analysis."
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Parse residues_group from string
    try:
        residues_group_parsed = ast.literal_eval(args.residues)
        if not isinstance(residues_group_parsed, dict):
            raise ValueError("Residues argument must be a dictionary string.")
    except (ValueError, SyntaxError) as e:
        print(f"Error: Could not parse residues argument: {e}")
        print("Please ensure it's a valid dictionary string, e.g., \"{'DPPC': (['PO4'], ['C4B', 'C4A']), 'CHOL':(['ROH'], ['R5'])}\"")
        sys.exit(1)

    print("\n--- Initializing MDAnalysis Universe ---")
    try:
        u = mda.Universe(args.gro_file, args.xtc_file)
    except Exception as e:
        print(f"Error loading MDAnalysis Universe: {e}")
        print("Please check if GRO/XTC files exist and are valid.")
        sys.exit(1)

    print("\n--- Running Height Analysis ---")
    height_analysis = Height(
        u,
        residues_group_parsed,
        k=args.k_value,
        file_path=args.output_csv,
        parallel=args.parallel,
        n_jobs=args.n_jobs
    )
    height_analysis.run(
        start=args.start_frame,
        stop=args.stop_frame,
        step=args.step_frame,
        verbose=args.verbose
    )

    print("\n--- Analysis Finished ---")
    print(f"Results saved to: {args.output_csv}")
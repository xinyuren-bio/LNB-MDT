import warnings
import os
import sys
import argparse 
import ast      
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


# --- Command-line Argument Parsing ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform Anisotropy analysis on molecular dynamics trajectories."
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
        default="cases/csv/anisotropy_results.csv",
        help="Path to the output CSV file for anisotropy results."
    )
    parser.add_argument(
        "--residues",
        type=str,
        default="{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}",
        help="A dictionary string defining residue groups for analysis. E.g., \"{'DPPC': ['PO4'], 'CHOL': ['ROH']}\""
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Enable parallel processing for anisotropy calculation."
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
        print("Please ensure it's a valid dictionary string, e.g., \"{'DPPC': ['PO4'], 'CHOL': ['ROH']}\"")
        sys.exit(1)

    print("\n--- Initializing MDAnalysis Universe ---")
    try:
        u = mda.Universe(args.gro_file, args.xtc_file)
    except Exception as e:
        print(f"Error loading MDAnalysis Universe: {e}")
        print("Please check if GRO/XTC files exist and are valid.")
        sys.exit(1)

    print("\n--- Running Anisotropy Analysis ---")
    anisotropy_analysis = Anisotropy(
        u,
        residues_group_parsed,
        file_path=args.output_csv,
        parallel=args.parallel,
        n_jobs=args.n_jobs
    )
    anisotropy_analysis.run(
        start=args.start_frame,
        stop=args.stop_frame,
        step=args.step_frame,
        verbose=args.verbose
    )

    print("\n--- Analysis Finished ---")
    print(f"Results saved to: {args.output_csv}")


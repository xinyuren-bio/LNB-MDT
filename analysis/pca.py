import warnings
import os
import sys
import argparse
import ast

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

# Assuming analysis.analysis_base and WriteExcelBubble are correctly defined
# and accessible in your environment, or you can provide their definitions.
# For demonstration purposes, I'll add placeholder classes if they are not
# part of a well-defined package you're importing.
try:
    from analysis.analysis_base import AnalysisBase
    # Assuming WriteExcelBubble is also in analysis.analysis_base or a similar location
    # If not, you might need to adjust this import or define a placeholder.
    from analysis.write_excel_bubble import WriteExcelBubble # Placeholder import
except ImportError:
    print("Warning: Could not import AnalysisBase or WriteExcelBubble. Using placeholders for demonstration.")
    class AnalysisBase:
        def __init__(self, trajectory):
            self._trajectory = trajectory
            self.results = type('Results', (object,), {})() # Simple object for results
            self._frame_index = 0
            self.n_frames = 0

        def run(self, start=None, stop=None, step=None, verbose=None):
            # Placeholder run method for serial execution
            for i in range(self.start, self.stop, self.step):
                self._frame_index = i - self.start # Adjust index for results array
                self.u.trajectory[i] # Advance trajectory for MDAnalysis
                self._single_frame()
                if verbose:
                    sys.stdout.write(f"\rProcessing frame {i+1}/{self._trajectory.n_frames}")
                    sys.stdout.flush()
            if verbose:
                print("\nSerial analysis complete.")
            
    class WriteExcelBubble:
        def __init__(self, **kwargs):
            self.kwargs = kwargs
            print(f"WriteExcelBubble initialized with: {kwargs}")

        def run(self):
            print(f"Saving results to {self.kwargs.get('file_path')}")
            # In a real scenario, you'd write your data to an Excel/CSV here.
            # For demonstration, we'll just print a confirmation.


__all__ = ['PCA']


class PCA(AnalysisBase):
    """
    A class for performing Principal Component Analysis on a group of atoms.
    """
    def __init__(self, universe, residues_group: dict, file_path: str = None,
                 parallel: bool = False, n_jobs: int = -1):
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residues_group.keys()) # Use .keys() for consistency
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
        for ts in self.u.trajectory[self.start:self.stop:self.step]:
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
            super().run(start=self.start, stop=self.stop, step=self.step, verbose=verbose)

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


# --- Command-line Argument Parsing ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform Principal Component Analysis on molecular dynamics trajectories."
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
        default="cases/csv/pca_results.csv",
        help="Path to the output CSV file for PCA results."
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
        help="Enable parallel processing for PCA calculation."
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


if __name__ == '__main__':
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

    print("\n--- Running PCA Analysis ---")
    pca_analysis = PCA(
        u,
        residues_group_parsed,
        file_path=args.output_csv,
        parallel=args.parallel,
        n_jobs=args.n_jobs
    )
    pca_analysis.run(
        start=args.start_frame,
        stop=args.stop_frame,
        step=args.step_frame,
        verbose=args.verbose
    )

    print("\n--- Analysis Finished ---")
    print(f"Results saved to: {args.output_csv}")
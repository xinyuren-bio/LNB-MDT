import os
import sys
import argparse
import ast
import time
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import MDAnalysis as mda
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

__all__ = ['Density']

class Density(AnalysisBase):
    def __init__(self, universe, ResiudeGroup: dict, GasGroup: dict, radius: float = 50, filePath: str = None, parallel: bool = False, n_jobs: int = -1):
        self.radius = radius
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(ResiudeGroup)
        self.gas_group = list(GasGroup)
        self.file_path = filePath
        self.parallel = parallel
        self.n_jobs = n_jobs
        self.headSp = {sp: ' '.join(ResiudeGroup[sp]) for sp in ResiudeGroup}
        self.gas_sp = {sp: ' '.join(GasGroup[sp]) for sp in GasGroup}

        self.headAtoms = self.u.atoms[[]]
        self.gas_atoms = self.u.atoms[[]]
        for sp in self.residues:
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]), updating=False)
        for sp in self.gas_group:
            self.gas_atoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.gas_sp[sp]), updating=False)

        self._n_residues = self.headAtoms.n_residues
        self.results.DensityWithtimes = None
        self.parameters = str(GasGroup) + f'Radius:{self.radius}'

    def _prepare(self):
        self.results.DensityWithtimes = np.full([self.n_frames], fill_value=np.nan)

    def _single_frame(self):
        self.results.DensityWithtimes[self._frame_index] = self._calculate_density_with_times()

    def _calculate_density_with_times(self):
        bubble_center = self.headAtoms.center_of_mass()
        gas = self.gas_atoms.select_atoms(f'point {bubble_center[0]} {bubble_center[1]} {bubble_center[2]} {self.radius}').n_residues
        bubble_volume = 4 / 3 * np.pi * (self.radius ** 3)
        return 14 / bubble_volume * gas

    def run(self, start=None, stop=None, step=None, verbose=None, callBack=None):
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
                delayed(Density._calculate_density_with_times)() for inputs in inputs_generator
            )
            if results_list:
                self.results.DensityWithtimes = np.array(results_list).T
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose, callBack=callBack)

        self._conclude()

    def _conclude(self):
        if self.file_path:
            dict_parameter = {
                'frames': list(range(self.start, self.stop, self.step)),
                'results': self.results.DensityWithtimes,
                'file_path': self.file_path,
                'description': 'Density With Times',
                'parameters': self.parameters
            }
            # Assuming WriteExcelBubble is defined and available
            WriteExcelBubble(**dict_parameter).run()


# --- Command-line Argument Parsing ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform Density Analysis on molecular dynamics trajectories."
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
        default="cases/csv/density_results.csv",
        help="Path to the output CSV file for Density results."
    )

    parser.add_argument(
        "--residues",
        type=str,
        default="{'DPPC': ['PO4'], 'DUPC': ['PO4'], 'CHOL': ['ROH']}",
        help="A dictionary string defining residue groups for analysis. E.g., \"{'DPPC': ['PO4'], 'CHOL': ['ROH']}\""
    )

    parser.add_argument(
        "--gas-group",
        type=str,
        default="{'N2': ['N2']}",
        help="A dictionary string defining residue groups for analysis. E.g., \"{'DPPC': ['PO4'], 'CHOL': ['ROH']}\""
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=50,
        help="Radius for Density calculation."
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Enable parallel processing for Density calculation."
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=2,
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
    start_time = time.time()

    # Parse residues_group from string
    try:
        residues_group_parsed = ast.literal_eval(args.residues)
        gas_group_parsed = ast.literal_eval(args.gas_group)
        if not isinstance(residues_group_parsed, dict):
            raise ValueError("Residues argument must be a dictionary string.")
        if not isinstance(gas_group_parsed, dict):
            raise ValueError("Gas group argument must be a dictionary string.")
    except (ValueError, SyntaxError) as e:
        print(f"Error: Could not parse gas_group argument: {e}")
        print("Please ensure it's a valid dictionary string, e.g., \"{'N2': ['N2']}\"")
        sys.exit(1)

    # Display selected residues and atoms
    print("Selected Gas Group:")
    for residue, atoms in residues_group_parsed.items():
        print(f"  {residue}: {atoms}")

    try:
        u = mda.Universe(args.gro_file, args.xtc_file)
    except Exception as e:
        print(f"Error loading MDAnalysis Universe: {e}")
        print("Please check if GRO/XTC files exist and are valid.")
        sys.exit(1)

    print(f"\nRunning Density Analysis on {u.trajectory.n_frames} frames...")
    density_analysis = Density(
        u,
        ResiudeGroup=residues_group_parsed,
        GasGroup=gas_group_parsed,
        filePath=args.output_csv,
        parallel=args.parallel,
        n_jobs=args.n_jobs,
        radius=args.radius
    )
    density_analysis.run(
        start=args.start_frame,
        stop=args.stop_frame,
        step=args.step_frame,
        verbose=args.verbose
    )

    end_time = time.time()
    elapsed_time = end_time - start_time
    
    print(f"\nAnalysis completed in {elapsed_time:.2f} seconds")
    print(f"Results saved to: {args.output_csv}")
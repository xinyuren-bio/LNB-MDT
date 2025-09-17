import warnings
import os
import sys
import argparse
import ast
warnings.filterwarnings('ignore')

import numpy as np
import scipy.sparse
import scipy.stats
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
from joblib import Parallel, delayed

if __name__ == '__main__':
    current_file_path = os.path.abspath(__file__)
    current_dir = os.path.dirname(current_file_path)
    package_root = os.path.abspath(os.path.join(current_dir, '..'))
    if package_root not in sys.path:
        sys.path.insert(0, package_root)
from analysis.analysis_base import * 

__all__ = ['NCluster']


def make_data(residues_group):
    ls1 = []
    for residue, atoms in residues_group.items():
        str_atom = ' '.join(atoms)
        ls1.append(f'resname {residue} and name {str_atom}')
    return ls1


class NCluster(AnalysisBase):
    """
    A class for counting the number of lipid clusters above a size threshold.
    """
    def __init__(self, universe, residues_group: dict, filePath: str = None, N_cutoff: int = 10, cutoff: float = 12.0,
                 parallel: bool = False, n_jobs: int = -1):
        super().__init__(universe.trajectory)
        self.u = universe
        self.cutoff = cutoff
        self.N_cutoff = N_cutoff
        self.file_path = filePath
        self.parallel = parallel
        self.n_jobs = n_jobs

        sel_residues_and_atoms = make_data(residues_group)
        self.atoms_search = self.u.atoms[[]]

        for str_sel in sel_residues_and_atoms:
            self.atoms_search += self.u.select_atoms(str_sel, updating=False)

        if self.atoms_search.n_atoms == 0:
            raise ValueError('No atoms found in the selection. Please check your residues_group.')

        # Create a dense, 0-based index for each residue in the selection
        self.atom_residue_indices = (
            scipy.stats.rankdata(
                self.atoms_search.resindices,
                method="dense",
            ) - 1
        )
        self.n_selection_residues = self.atoms_search.n_residues
        
        self.results.NCluster = None
        self.parameters = f"{residues_group}, N_cutoff:{self.N_cutoff}, Cutoff:{self.cutoff}"

    @property
    def NCluster(self):
        return self.results.NCluster

    def _prepare(self):
        self.results.NCluster = np.full(self.n_frames, 0, dtype=int)

    @staticmethod
    def _calculate_ncluster_for_frame(positions, box, cutoff, N_cutoff, atom_residue_indices, n_selection_residues):
        if n_selection_residues <= 1:
            return 1 if N_cutoff <= 1 and n_selection_residues > 0 else 0

        pairs = capped_distance(
            positions,
            positions,
            max_cutoff=cutoff,
            box=box,
            return_distances=False
        )
        
        if pairs.shape[0] == 0:
            return n_selection_residues if N_cutoff <= 1 else 0

        residue_index_pairs = np.unique(atom_residue_indices[pairs], axis=0)
        
        if residue_index_pairs.shape[0] == 0:
            return n_selection_residues if N_cutoff <=1 else 0

        ref, nei = residue_index_pairs[residue_index_pairs[:, 0] != residue_index_pairs[:, 1]].T
        data = np.ones_like(ref, dtype=np.int8)

        neighbours_frame = scipy.sparse.csr_matrix(
            (data, (ref, nei)),
            shape=(n_selection_residues, n_selection_residues)
        )

        _, com_labels = scipy.sparse.csgraph.connected_components(neighbours_frame, directed=False)
        _, counts = np.unique(com_labels, return_counts=True)
        
        clusters_above_cutoff = counts[counts > N_cutoff]
        return len(clusters_above_cutoff)

    def _single_frame(self):
        n_clusters = NCluster._calculate_ncluster_for_frame(
            self.atoms_search.positions,
            self._ts.dimensions,
            self.cutoff,
            self.N_cutoff,
            self.atom_residue_indices,
            self.n_selection_residues
        )
        self.results.NCluster[self._frame_index] = n_clusters

    def _get_inputs_generator(self):
        for ts in self.u.trajectory[self.start:self.stop:self.step]:
            yield (
                self.atoms_search.positions.copy(),
                ts.dimensions,
                self.cutoff,
                self.N_cutoff,
                self.atom_residue_indices,
                self.n_selection_residues
            )

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
                delayed(NCluster._calculate_ncluster_for_frame)(*inputs) for inputs in inputs_generator
            )
            if results_list:
                self.results.NCluster = np.array(results_list)
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose, callBack=callBack)

        self._conclude()

    def _conclude(self):
        if self.file_path:
            dict_parameter = {
                'frames': list(range(self.start, self.stop, self.step)),
                'results': self.results.NCluster,
                'file_path': self.file_path,
                'description': 'Number of Clusters',
                'parameters': self.parameters
            }
            # Assuming WriteExcelBubble is defined and available
            WriteExcelBubble(**dict_parameter).run()
            print(f"Analysis complete. Results will be saved to {self.file_path}")


# --- Command-line Argument Parsing ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform N-Cluster analysis on molecular dynamics trajectories."
    )

    parser.add_argument(
        "--gro-file", "-g",
        type=str,
        default="cases/lnb.gro",
        help="Path to the GRO file (topology file)."
    )
    parser.add_argument(
        "--xtc-file", "-x",
        type=str,
        default="cases/md.xtc",
        help="Path to the XTC file (trajectory file)."
    )
    parser.add_argument(
        "--output-csv", "-o",
        type=str,
        default="cases/csv/ncluster_results.csv",
        help="Path to the output CSV file for n-cluster results."
    )
    parser.add_argument(
        "--residues", "-r",
        type=str,
        default="{'DAPC': ['GL1', 'GL2'], 'DPPC': ['PO4']}",
        help="A dictionary string defining residue groups for analysis. E.g., \"{'DAPC': ['GL1', 'GL2'], 'DPPC': ['PO4']}\""
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=12.0,
        help="Cutoff distance for clustering (in Angstroms)."
    )
    parser.add_argument(
        "--n-cutoff",
        type=int,
        default=10,
        help="Minimum cluster size threshold."
    )
    parser.add_argument(
        "--parallel", "-p",
        action="store_true",
        help="Enable parallel processing for n-cluster calculation."
    )
    parser.add_argument(
        "--n-jobs", "-j",
        type=int,
        default=-1,
        help="Number of jobs to run in parallel. -1 means using all available CPU cores."
    )
    parser.add_argument(
        "--start-frame", "-s",
        type=int,
        default=0,
        help="Starting frame for analysis (0-indexed)."
    )
    parser.add_argument(
        "--stop-frame", "-e",
        type=int,
        help="Stopping frame for analysis (exclusive). Defaults to end of trajectory."
    )
    parser.add_argument(
        "--step-frame", "-t",
        type=int,
        default=1,
        help="Step size for frames during analysis."
    )
    parser.add_argument(
        "--verbose", "-v",
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
        print("Please ensure it's a valid dictionary string, e.g., \"{'DAPC': ['GL1', 'GL2'], 'DPPC': ['PO4']}\"")
        sys.exit(1)

    print("\n--- Initializing MDAnalysis Universe ---")
    try:
        u = mda.Universe(args.gro_file, args.xtc_file)
    except Exception as e:
        print(f"Error loading MDAnalysis Universe: {e}")
        print("Please check if GRO/XTC files exist and are valid.")
        sys.exit(1)

    print("\n--- Running N-Cluster Analysis ---")
    ncluster_analysis = NCluster(
        u,
        residues_group_parsed,
        cutoff=args.cutoff,
        N_cutoff=args.n_cutoff,
        file_path=args.output_csv,
        parallel=args.parallel,
        n_jobs=args.n_jobs
    )
    ncluster_analysis.run(
        start=args.start_frame,
        stop=args.stop_frame,
        step=args.step_frame,
        verbose=args.verbose
    )

    print("\n--- Analysis Finished ---")
    print(f"Results saved to: {args.output_csv}")
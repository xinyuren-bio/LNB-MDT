import warnings
import os
import sys
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
    def __init__(self, universe, residues_group: dict, file_path: str = None, N_cutoff: int = 10, cutoff: float = 12.0,
                 parallel: bool = False, n_jobs: int = -1):
        super().__init__(universe.trajectory)
        self.u = universe
        self.cutoff = cutoff
        self.N_cutoff = N_cutoff
        self.file_path = file_path
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
                delayed(NCluster._calculate_ncluster_for_frame)(*inputs) for inputs in inputs_generator
            )
            if results_list:
                self.results.NCluster = np.array(results_list)
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose)

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


if __name__ == "__main__":
    gro_file = "cases/lnb.gro"
    xtc_file = "cases/md.xtc"
    csv_file_serial = "cases/csv/ncluster_serial.csv"
    csv_file_parallel = "cases/csv/ncluster_parallel.csv"

    u = mda.Universe(gro_file, xtc_file)
    residues_group = {'DAPC': ['GL1', 'GL2'], 'DPPC': ['PO4']}

    #1. 串行
    analysis_serial = NCluster(
        u, 
        residues_group, 
        cutoff=12.0, 
        N_cutoff=10,
        file_path=csv_file_serial, 
        parallel=False
    )
    analysis_serial.run(verbose=True)

    #2. 并行
    analysis_parallel = NCluster(
        u, 
        residues_group, 
        cutoff=12.0,
        N_cutoff=10,
        file_path=csv_file_parallel, 
        parallel=True, 
        n_jobs=-1
    )
    analysis_parallel.run(verbose=True)
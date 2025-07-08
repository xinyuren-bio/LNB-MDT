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

__all__ = ['Cluster']


class Cluster(AnalysisBase):
    """
    A class for analyzing lipid clustering in a bilayer system.
    """

    def __init__(self, universe, residues_group: dict, cutoff: float = 8.0, file_path: str = None,
                 parallel: bool = False, n_jobs: int = -1):
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residues_group)
        self.cutoff = cutoff
        self.file_path = file_path
        self.parallel = parallel
        self.n_jobs = n_jobs

        self.atomSp = {sp: ' '.join(residues_group[sp]) for sp in self.residues}
        print("Atoms for clustering:", self.atomSp)

        self.headAtoms = self.u.atoms[[]]
        for sp in self.residues:
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.atomSp[sp]), updating=False)

        self.n_residues = self.headAtoms.n_residues
        self.atom_resindices = self.headAtoms.resindices
        self.residue_names = self.headAtoms.residues.resnames

        self.results.LargestClusterSize = None
        self.results.LargestClusterResidues = None
        self.parameters = str(residues_group) + 'Cutoff:' + str(self.cutoff)

    @property
    def LargestClusterSize(self):
        return self.results.LargestClusterSize

    @property
    def LargestClusterResidues(self):
        return self.results.LargestClusterResidues

    def _prepare(self):
        self.results.LargestClusterSize = np.full(self.n_frames, 0, dtype=int)
        self.results.LargestClusterResidues = np.full(self.n_frames, fill_value=None, dtype=object)

    @staticmethod
    def _calculate_cluster_for_frame(positions, box, cutoff, atom_resindices, residue_names, n_residues):
        pairs = capped_distance(
            positions,
            positions,
            max_cutoff=cutoff,
            box=box,
            return_distances=False
        )

        if pairs.shape[0] == 0:
            return 1 if n_residues > 0 else 0, np.array([]) if n_residues > 0 else np.array([])

        residue_indices_pairs = np.unique(atom_resindices[pairs], axis=0)

        if residue_indices_pairs.shape[0] == 0:
             return 1 if n_residues > 0 else 0, np.array([]) if n_residues > 0 else np.array([])

        ref, nei = residue_indices_pairs[residue_indices_pairs[:, 0] != residue_indices_pairs[:, 1]].T
        data = np.ones_like(ref, dtype=np.int8)

        neighbours_frame = scipy.sparse.csr_matrix(
            (data, (ref, nei)),
            shape=(n_residues, n_residues)
        )

        n_components, com_labels = scipy.sparse.csgraph.connected_components(neighbours_frame, directed=False)

        if n_components == 0:
            return 0, np.array([])
        
        unique_com_labels, counts = np.unique(com_labels, return_counts=True)
        largest_cluster_index = np.argmax(counts)
        largest_label = unique_com_labels[largest_cluster_index]
        max_size = counts[largest_cluster_index]

        residues_in_cluster = residue_names[com_labels == largest_label]

        return max_size, residues_in_cluster

    def _single_frame(self):
        size, residues = Cluster._calculate_cluster_for_frame(
            self.headAtoms.positions,
            self._ts.dimensions,
            self.cutoff,
            self.atom_resindices,
            self.residue_names,
            self.n_residues
        )
        self.results.LargestClusterSize[self._frame_index] = size
        self.results.LargestClusterResidues[self._frame_index] = residues

    def _get_inputs_generator(self):
        for ts in self.u.trajectory[self.start:self.stop:self.step]:
            yield (
                self.headAtoms.positions.copy(),
                ts.dimensions,
                self.cutoff,
                self.atom_resindices,
                self.residue_names,
                self.n_residues
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
                delayed(Cluster._calculate_cluster_for_frame)(*inputs) for inputs in inputs_generator
            )
            if results_list:
                sizes, residues_list = zip(*results_list)
                self.results.LargestClusterSize = np.array(sizes)
                self.results.LargestClusterResidues = np.array(residues_list, dtype=object)
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose)

        self._conclude()

    def _conclude(self):
        if self.file_path:
            dict_parameter = {
                'frames': list(range(self.start, self.stop, self.step)),
                'results': self.results.LargestClusterSize,
                'file_path': self.file_path,
                'description': 'Largest Cluster Size',
                'parameters': self.parameters
            }
            WriteExcelBubble(**dict_parameter).run()
            print(f"Analysis complete. Results will be saved to {self.file_path}")


if __name__ == "__main__":
    gro_file = "cases/lnb.gro"
    xtc_file = "cases/md.xtc"
    csv_file_serial = "cases/csv/cluster_serial.csv"
    csv_file_parallel = "cases/csv/cluster_parallel.csv"

    u = mda.Universe(gro_file, xtc_file)
    residues_group = {'DPPC': ['PO4'], 'DUPC': ['PO4'], 'CHOL': ['ROH']}

    #1. 串行
    print("--- Running Serial Analysis ---")
    analysis_serial = Cluster(u, residues_group, cutoff=10.0, file_path=csv_file_serial, parallel=False)
    analysis_serial.run(verbose=True)
    
    #2. 并行
    print("--- Running Parallel Analysis ---")
    analysis_parallel = Cluster(u, residues_group, cutoff=10.0, file_path=csv_file_parallel, parallel=True, n_jobs=-1)
    analysis_parallel.run(verbose=True)
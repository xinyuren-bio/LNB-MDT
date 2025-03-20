import warnings
warnings.filterwarnings('ignore')

from MDAnalysis.lib.distances import capped_distance
import numpy as np
import scipy.sparse
import scipy.stats

try:
    from .analysis_base import *
except:
    from analysis_base import *

__all__ = ['NCluster']

# suppress some MDAnalysis warnings when writing PDB/GRO files

class NCluster(AnalysisBase):
    def __init__(self, universe, residues_group, file_path=None, N_cutoff=10, cutoff=12):
        super().__init__(universe.trajectory)

        self.u = universe

        self.cutoff = cutoff
        self.N_cutoff = N_cutoff

        self.file_path = file_path

        sel_residues_and_atoms = make_data(residues_group)
        self.atoms_search = self.u.atoms[[]]

        for str_sel in sel_residues_and_atoms:
            self.atoms_search += self.u.select_atoms(str_sel, updating=False)

        if self.atoms_search.n_residues <= 1:
            raise ValueError('No atoms found in the selection')

        self.resids_search_ori = self.atoms_search.resindices
        self.resids_residues_ori = self.atoms_search.residues.resids
        self.resids_search_sorted = (
            scipy.stats.rankdata(
                self.resids_search_ori,
                method="dense",
            )
            - 1
        )
        self.results.NCluster = None

    @property
    def NCluster(self):
        return self.results.NCluster


    def _prepare(self):
        self.results.NCluster = np.zeros([self.n_frames])

    def _single_frame(self):
        positions_single_frame = self.atoms_search.positions
        pairs = capped_distance(
            positions_single_frame
            , positions_single_frame
            , max_cutoff=self.cutoff
            , box=self._ts.dimensions
            , return_distances=False
        )

        ref, nei = np.unique(self.resids_search_sorted[pairs], axis=0).T

        index_needs = ref != nei

        ref = ref[index_needs]
        nei = nei[index_needs]

        data = np.ones_like(ref)

        neighbours_frame = scipy.sparse.csr_matrix(
            (data, (ref, nei))
            , dtype=np.int8
            , shape=(self.atoms_search.n_residues, self.atoms_search.n_residues)
        )

        _, com_labels = scipy.sparse.csgraph.connected_components(neighbours_frame)

        unique_com_labels, counts = np.unique(com_labels, return_counts=True)
        clusters_above_cutoff = counts[counts > self.N_cutoff]
        n_clusters = len(clusters_above_cutoff)

        self.results.NCluster[self._frame_index] = n_clusters


    def _conclude(self):
        if self.file_path:
            dict_parameter = {'step': self.step, 'n_frames': self.n_frames,
                             'results': self.results.NCluster,
                              'file_path': self.file_path, 'description': 'Number of Cluster',
                              'value_divition': 1}
            WriteExcelBubble(**dict_parameter).run()


def make_data(residues_group):
    ls1 = []
    for residue, atoms in residues_group.items():
        str_atom = ' '.join(atoms)
        ls1.append(f'resname {residue} and name {str_atom}')
    return ls1


if __name__ == "__main__":
    import MDAnalysis as mda
    u = mda.Universe('E:/ach.gro', 'E:/ach.xtc')
    cls1 = NCluster(u, {'DAPC': ['GL1', 'GL2']}, cutoff=12, file_path='E:/untitled2.csv', N_cutoff=10)
    cls1.run(0, 100)

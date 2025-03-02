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

__all__ = ['Cluster']


# suppress some MDAnalysis warnings when writing PDB/GRO files


class Cluster(AnalysisBase):
    def __init__(self, universe, residues_group, file_path=None, cutoff=12):
        super().__init__(universe.trajectory)

        self.u = universe

        self.cutoff = cutoff

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
        self.results.Cluster = None
        self.results.Resindices = None

    @property
    def Cluster(self):
        return self.results.Cluster

    @property
    def Resindices(self):
        return self.results.Resindices

    def _prepare(self):
        self.results.Cluster = np.zeros([self.n_frames])
        self.results.Resindices = np.full(self.n_frames, fill_value=0, dtype=object)

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
        largest_label = unique_com_labels[np.argmax(counts)]
        self.results.Cluster[self._frame_index] = max(counts)
        self.results.Resindices[self._frame_index] = self.resids_residues_ori[com_labels == largest_label]

    def _conclude(self):
        if self.file_path:
            dict_parameter = {'step': self.step, 'n_frames': self.n_frames,
                             'results': self.results.Cluster,
                              'file_path': self.file_path, 'description': 'Cluster',
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
    cls1 = Cluster(u, {'DAPC': ['NC3']}, cutoff=20, file_path='E:/untitled1.csv')
    cls1.run(0, 100)

import warnings
warnings.filterwarnings('ignore')

import MDAnalysis as mda
import numpy as np
from scipy.spatial import KDTree
from scipy.linalg import eigh

try:
    from .analysis_base import *
except:
    from analysis_base import *


__all__ = ['Height']


class Height(AnalysisBase):
    """
    A class for calculating the height of lipid molecules in LNB system.
    
    This class analyzes the vertical distance between the head and tail groups of lipid molecules
    in an MD trajectory. It uses local normal vectors to account for membrane
    curvature and lipid tilt.
    """

    def __init__(self, universe, residuesGroup: dict, k: int = None, file_path: str = None):
        """
        Initialize the Height analysis class.
        
        Parameters
        ----------
        universe : MDAnalysis.Universe
            The MDAnalysis Universe object containing the molecular dynamics trajectory.
            This should include both the structure file (e.g., .gro) and trajectory file (e.g., .xtc).
            
        residuesGroup : dict
            A dictionary specifying the head and tail atoms for each lipid type.
            Format: {'lipid_name': (['head_atom_names'], ['tail_atom_names'])}
            Example: {
                'DPPC': (['PO4'], ['C4A', 'C4B']),
                'CHOL': (['ROH'], ['R5'])
            }
            
        k : int, optional
            The number of nearest neighbors to use for calculating local normal vectors.
            A larger k value results in smoother normal vectors but requires more computation.
            If None, a default value will be used.
            
        file_path : str, optional
            The path where the analysis results will be saved as a CSV file.
            If None, results will not be saved to disk.
            
        Attributes
        ----------
        headAtoms : MDAnalysis.AtomGroup
            The selected head atoms of all lipid molecules.
            
        tailAtoms : MDAnalysis.AtomGroup
            The selected tail atoms of all lipid molecules.
            
        _n_residues : int
            The total number of lipid molecules in the system.
            
        resids : numpy.ndarray
            The residue IDs of all lipid molecules.
            
        resnames : numpy.ndarray
            The residue names of all lipid molecules.
            
        resArrange : numpy.ndarray
            An array used to maintain the original order of atoms.
            
        results.Height : numpy.ndarray
            A 2D array storing the height values for all lipid molecules across all frames.
            Shape: (n_residues, n_frames)
        """
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residuesGroup)
        self.k = k
        self.file_path = file_path

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

        self.parameters = str(residuesGroup)+ 'K:' + str(self.k)
    def _prepare(self):
        self.results.Height = np.full([self._n_residues, self.n_frames],
                                      fill_value=np.nan)

    def _single_frame(self):
        centerHeadSp = self.headAtoms.positions
        normals = get_normals(self.k, centerHeadSp)
        centerTailSp = self.tailAtoms.center_of_geometry(compound='residues')
        head_to_tail = centerHeadSp - centerTailSp
        distance = (np.abs(np.einsum('ij,ij->i', normals, head_to_tail)))[self.resArrange]
        self.results.Height[:, self._frame_index] = distance * 0.1  # A° to nm

    @property
    def Height(self):
        return self.results.Height

    def _conclude(self):
        if self.file_path:
            lipids_ratio = {sp: self.u.select_atoms(f'resname {sp}').n_residues for sp in self.residues}
            dict_parameter = {
                'frames': [i for i in range(self.start, self.stop, self.step)],
                'resids': self.resids,
                'resnames': self.resnames,
                'positions': self.headAtoms.positions,
                'results': self.results.Height,
                'file_path': self.file_path,
                'description': 'Height',
                'parameters': self.parameters,
                'lipids_type': lipids_ratio
            }
            WriteExcelLipids(**dict_parameter).run()
        else:
            return self.results.Height


if __name__ == "__main__":
    import time
    t1 = time.time()
    gro_file = "cases/lnb.gro"
    xtc_file = "cases/md.xtc"
    csv_file = "cases/csv/height_step5_lnb.csv"
    u = mda.Universe(gro_file, xtc_file)
    dict_residue = {'DPPC': (['PO4'], ['C4B', 'C4A']), 'DUPC':(['PO4'], ['C3A', 'C4B']), 'CHOL':(['ROH'], ['R5'])}
    cls2 = Height(u, dict_residue, k=21, file_path=csv_file)
    cls2.run(verbose=True)
    t2 = time.time()
    print(t2 - t1)
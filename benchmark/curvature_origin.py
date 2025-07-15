import warnings
warnings.filterwarnings('ignore')

import MDAnalysis as mda
import numpy as np
from numpy.linalg import lstsq
from scipy.spatial import KDTree

try:
    from .analysis_base import *
except:
    from analysis_base import *

__all__ = ['Curvature']


class Curvature(AnalysisBase):
    """
    A class for calculating the mean curvature of a lipid bilayer.
    
    This class analyzes the local curvature of a lipid bilayer by calculating
    the mean curvature at each point using the positions of lipid head groups.
    The curvature is calculated using the local normal vectors and their derivatives.
    """
    
    def __init__(self, universe, residueGroup: dict, k: int = None, file_path: str = None, method: str = 'mean'):
        """
        Initialize the Curvature analysis class.
        
        Parameters
        ----------
        universe : MDAnalysis.Universe
            The MDAnalysis Universe object containing the molecular dynamics trajectory.
            This should include both the structure file (e.g., .gro) and trajectory file (e.g., .xtc).
            
        residueGroup : dict
            A dictionary specifying the head atoms for each lipid type.
            Format: {'lipid_name': ['head_atom_names']}
            Example: {
                'DPPC': ['PO4'],
                'CHOL': ['ROH']
            }
            
        k : int, optional
            The number of nearest neighbors to use for calculating local curvature.
            A larger k value results in smoother curvature but requires more computation.
            If None, a default value will be used.
            
        file_path : str, optional
            The path where the analysis results will be saved as a CSV file.
            If None, results will not be saved to disk.
            
        method : str, optional
            The method to use for curvature calculation.
            Options: 'mean' (default) or 'gaussian'
            - 'mean': calculates the mean curvature (H)
            - 'gaussian': calculates the Gaussian curvature (K)
            
        Attributes
        ----------
        headAtoms : MDAnalysis.AtomGroup
            The selected head atoms of all lipid molecules.
            
        _n_residues : int
            The total number of lipid molecules in the system.
            
        resids : numpy.ndarray
            The residue IDs of all lipid molecules.
            
        resnames : numpy.ndarray
            The residue names of all lipid molecules.
            
        results.Curvature : numpy.ndarray
            A 2D array storing the curvature values for all lipid molecules across all frames.
            Shape: (n_residues, n_frames)
            
        method : str
            The method used for curvature calculation ('mean' or 'gaussian').
        """
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residueGroup)
        self.k = k
        self.file_path = file_path
        self.method = method

        # Convert head atom names to space-separated strings
        self.headSp = {sp: ' '.join(residueGroup[sp]) for sp in residueGroup}
        print("Head atoms:", self.headSp)

        # Initialize atom selection
        self.headAtoms = self.u.atoms[[]]

        # Select head atoms for all specified lipid types
        for i in range(len(self.residues)):
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (self.residues[i], self.headSp[self.residues[i]]), updating=False)

        # Set basic attributes
        self._n_residues = self.headAtoms.n_residues
        self.resids = self.headAtoms.resids
        self.resnames = self.headAtoms.resnames
        self.results.Curvature = None

        # Record analysis parameters
        self.parameters = str(residueGroup) + 'K:' + str(self.k) + 'Method:' + method

    @property
    def MeanCurvature(self):
        return self.results.MeanCurvature

    @property
    def GussainCurvature(self):
        return self.results.GussainCurvature

    def _prepare(self):
        self.results.MeanCurvature = np.full([self._n_residues, self.n_frames]
                                             , fill_value=np.nan)

        self.results.GussainCurvature = np.full([self._n_residues, self.n_frames]
                                                , fill_value=np.nan)

        self.results.Normal = np.full([self._n_residues, self.n_frames, 3],
                                      fill_value=np.nan)

    def _single_frame(self):
        head_positions = self.headAtoms.positions
        kd_tree = KDTree(head_positions)
        center_head = self.headAtoms.center_of_geometry()
        print(np.mean(np.linalg.norm(head_positions - center_head, axis=-1)))
        center_to_head = center_head - head_positions
        for i, point in enumerate(head_positions):
            _, idxs = kd_tree.query(point, self.k)
            nearest_neighbors = head_positions[idxs]
            normals = fit_plane_and_get_normal(nearest_neighbors)
            if np.dot(center_to_head[i,], normals[2]) < 0:
                normals[2] = -normals[2]
                xxx = cross_product(normals[0], normals[1])
                if np.dot(xxx, normals[2]) < 0:
                    normals[[0, 1]] = normals[[1, 0]]
            self.results.Normal[i, self._frame_index] = normals[2]
            mean = np.mean(nearest_neighbors, axis=0)
            centered_points = nearest_neighbors - mean

            local_points = (np.dot(normals, centered_points.T)).T

            coefficients = fit_quadratic_surface_1(local_points)

            first_style = calculate_coefficients(coefficients)

            if self.method == 'mean':
                k_mean = calculate_mean_curvature(first_style)
                self.results.MeanCurvature[i, self._frame_index] = k_mean * 10

            elif self.method == 'gussain':
                k_gussain = calculate_gaussian_curvature(first_style)
                self.results.GussainCurvature[i, self._frame_index] = k_gussain * 10

    def _conclude(self):
        if self.file_path:
            lipids_ratio = {sp: self.u.select_atoms(f'resname {sp}').n_residues for sp in self.residues}
            dict_parameter = {
                'frames': [i for i in range(self.start, self.stop, self.step)]
                , 'resids': self.resids
                , 'resnames': self.resnames
                , 'positions': self.headAtoms.positions
                , 'results': self.results.MeanCurvature
                , 'file_path': self.file_path
                , 'description': 'Mean Curvature(nm -1)'
                , 'parameters': self.parameters
                , 'lipids_type': lipids_ratio
            }
            WriteExcelLipids(**dict_parameter).run()


def fit_plane_and_get_normal(points):
    mean = np.mean(points, axis=0)
    centered_points = points - mean
    covariance_matrix = np.cov(centered_points.T)
    eigenvalues, eigenvectors = np.linalg.eig(covariance_matrix)
    eigenvalue_index = np.argsort(eigenvalues)
    normal_x = eigenvectors[:, eigenvalue_index[1]]
    normal_y = eigenvectors[:, eigenvalue_index[2]]
    normal_z = eigenvectors[:, eigenvalue_index[0]]
    normal_arr = np.vstack((normal_x, normal_y, normal_z))
    return normal_arr

def find_k_nearest_neighbors_and_normals(k, particles):
    kdtree = KDTree(particles)
    normals = np.zeros((particles.shape[0], 3))
    for i, point in enumerate(particles):
        dists, idxs = kdtree.query(point, k + 1)
        nearest_neighbors = particles[idxs[1:]]
        if nearest_neighbors.shape[0] >= k:
            nearest_neighbors = np.vstack((nearest_neighbors, point))
            normal = fit_plane_and_get_normal(nearest_neighbors)
            normals[i] = normal
        else:
            normals[i] = np.zeros(3)
    return normals


def cross_product(v1, v2):
    x1, y1, z1 = v1
    x2, y2, z2 = v2
    cross_x = y1 * z2 - z1 * y2
    cross_y = z1 * x2 - x1 * z2
    cross_z = x1 * y2 - y1 * x2
    return np.hstack((cross_x, cross_y, cross_z))


def fit_quadratic_surface_1(point_colud):
    A = []
    b = []
    for point in point_colud:
        x, y, z = point
        A.append([x ** 2, y ** 2, x * y, x, y, 1])
        b.append(z)
    A = np.array(A)
    b = np.array(b)
    coefficients, _, _, _ = lstsq(A, b)
    coefficients_points = coefficients
    return coefficients_points


def calculate_coefficients(coefficients):
    E = (coefficients[3] ** 2) + 1
    F = coefficients[3] * coefficients[4]
    G = (coefficients[4] ** 2) + 1
    L = 2 * coefficients[0]
    M = coefficients[2]
    N = coefficients[1] * 2
    EFGLMN = np.hstack((E, F, G, L, M, N))
    return EFGLMN


def calculate_mean_curvature(coefficients):
    EN = coefficients[0] * coefficients[5]
    FM = coefficients[1] * coefficients[4]
    GL = coefficients[2] * coefficients[3]
    EG = coefficients[0] * coefficients[2]
    F2 = coefficients[1] ** 2
    k_mean = (EN - 2 * FM + GL) / (2 * (EG - F2))
    return k_mean


def calculate_gaussian_curvature(coefficients):
    LN = coefficients[3] * coefficients[-1]
    M2 = coefficients[-2] ** 2
    EG = coefficients[0] * coefficients[2]
    F2 = coefficients[1] ** 2
    k_gaussian = (LN - M2) / (EG - F2)
    return k_gaussian


if __name__ == "__main__":
    gro_file = "cases/lnb.gro"
    xtc_file = "cases/md.xtc"
    csv_file = "cases/csv/area_step5_lnb.csv"
    u = mda.Universe(gro_file, xtc_file)
    cls = Curvature(u, {'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}, 20, file_path=csv_file)
    # 执行计算
    cls.run(0, 100, verbose=True)
    


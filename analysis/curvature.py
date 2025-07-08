import warnings
import os
import sys
warnings.filterwarnings('ignore')
import MDAnalysis as mda
import numpy as np
from numpy.linalg import eig, lstsq
from scipy.spatial import KDTree
from joblib import Parallel, delayed

if __name__ == '__main__':
    current_file_path = os.path.abspath(__file__)
    current_dir = os.path.dirname(current_file_path)
    package_root = os.path.abspath(os.path.join(current_dir, '..'))
    if package_root not in sys.path:
        sys.path.insert(0, package_root)
from analysis.analysis_base import * 

__all__ = ['Curvature']


class Curvature(AnalysisBase):
    """
    A class for calculating the mean and Gaussian curvature of a lipid bilayer.
    """
    def __init__(self, universe, residueGroup: dict, k: int = 20, file_path: str = None,
                 method: str = 'mean', parallel: bool = False, n_jobs: int = -1):
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residueGroup)
        self.k = k
        self.file_path = file_path
        self.method = method.lower()
        self.parallel = parallel
        self.n_jobs = n_jobs

        self.headSp = {sp: ' '.join(residueGroup[sp]) for sp in residueGroup}
        print("Head atoms:", self.headSp)

        self.headAtoms = self.u.atoms[[]]
        for res_name in self.residues:
            selection_str = f'resname {res_name} and name {self.headSp[res_name]}'
            self.headAtoms += self.u.select_atoms(selection_str, updating=False)

        self._n_residues = self.headAtoms.n_residues
        self.resids = self.headAtoms.resids
        self.resnames = self.headAtoms.resnames

        self.results.MeanCurvature = None
        self.results.GaussianCurvature = None
        self.results.Normal = None

        self.parameters = f"{residueGroup}, K:{self.k}, Method:{self.method}"

    @property
    def MeanCurvature(self):
        return self.results.MeanCurvature

    @property
    def GaussianCurvature(self):
        return self.results.GaussianCurvature

    def _prepare(self):
        self.results.MeanCurvature = np.full([self._n_residues, self.n_frames], fill_value=np.nan)
        self.results.GaussianCurvature = np.full([self._n_residues, self.n_frames], fill_value=np.nan)
        self.results.Normal = np.full([self._n_residues, self.n_frames, 3], fill_value=np.nan)

    @staticmethod
    def _calculate_curvature_for_frame(head_positions, k):
        n_residues = head_positions.shape[0]
        if n_residues < k:
            nan_array = np.full(n_residues, np.nan)
            nan_normals = np.full((n_residues, 3), np.nan)
            return (nan_array, nan_array, nan_normals)

        kd_tree = KDTree(head_positions)
        _, idxs = kd_tree.query(head_positions, k)
        nearest_neighbors = head_positions[idxs]

        means = np.mean(nearest_neighbors, axis=1, keepdims=True)
        centered_points = nearest_neighbors - means
        cov_matrices = np.einsum('nki,nkj->nij', centered_points, centered_points)
        eigenvalues, eigenvectors = eig(cov_matrices)

        eigenvalue_indices = np.argsort(eigenvalues, axis=1)
        n_range = np.arange(n_residues)

        normals_x = eigenvectors[n_range, :, eigenvalue_indices[:, 1]]
        normals_y = eigenvectors[n_range, :, eigenvalue_indices[:, 2]]
        normals_z = eigenvectors[n_range, :, eigenvalue_indices[:, 0]]

        center_head = np.mean(head_positions, axis=0)
        center_to_head = center_head - head_positions
        dots = np.einsum('ij,ij->i', normals_z, center_to_head)

        flip_mask = dots < 0
        normals_z[flip_mask] *= -1

        cross_products = np.cross(normals_x, normals_y, axis=1)
        handedness_dots = np.einsum('ij,ij->i', cross_products, normals_z)
        flip_handedness_mask = handedness_dots < 0
        normals_x[flip_handedness_mask], normals_y[flip_handedness_mask] = \
            normals_y[flip_handedness_mask], normals_x[flip_handedness_mask]

        normals_batch = np.stack((normals_x, normals_y, normals_z), axis=1)
        local_points = np.einsum('nki,nji->nkj', centered_points, normals_batch)

        k_mean = np.full(n_residues, np.nan)
        k_gaussian = np.full(n_residues, np.nan)

        for i in range(n_residues):
            x_i, y_i, z_i = local_points[i].T
            A_i = np.column_stack([x_i**2, y_i**2, x_i*y_i, x_i, y_i, np.ones_like(x_i)])
            
            try:
                coeffs, _, _, _ = lstsq(A_i, z_i, rcond=None)
            except np.linalg.LinAlgError:
                continue

            c0, c1, c2, c3, c4, _ = coeffs
            E, F, G = c3**2 + 1, c3 * c4, c4**2 + 1
            L, M, N = 2 * c0, c2, 2 * c1

            denominator_mean = 2 * (E * G - F**2)
            denominator_gauss = E * G - F**2

            if np.abs(denominator_mean) > 1e-9:
                k_mean[i] = (E * N - 2 * F * M + G * L) / denominator_mean

            if np.abs(denominator_gauss) > 1e-9:
                k_gaussian[i] = (L * N - M**2) / denominator_gauss
        
        return (k_mean * 10, k_gaussian * 10, normals_z)

    def _single_frame(self):
        head_positions = self.headAtoms.positions
        k_mean, k_gaussian, normals = self._calculate_curvature_for_frame(head_positions, self.k)

        self.results.MeanCurvature[:, self._frame_index] = k_mean
        self.results.GaussianCurvature[:, self._frame_index] = k_gaussian
        self.results.Normal[:, self._frame_index] = normals

    def _get_inputs_generator(self):
        for _ in self.u.trajectory[self.start:self.stop:self.step]:
            yield (self.headAtoms.positions.copy(), self.k)

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
                delayed(self._calculate_curvature_for_frame)(*inputs) for inputs in inputs_generator
            )

            if results_list:
                mean_curv_list, gauss_curv_list, normal_list = zip(*results_list)
                self.results.MeanCurvature = np.array(mean_curv_list).T
                self.results.GaussianCurvature = np.array(gauss_curv_list).T
                self.results.Normal = np.transpose(np.array(normal_list), (1, 0, 2))
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose)
        
        self._conclude()

    def _conclude(self):
        if self.file_path:
            if self.method == 'mean':
                results_to_save = self.results.MeanCurvature
                description = 'Mean Curvature (nm-1)'
            elif self.method == 'gaussian':
                results_to_save = self.results.GaussianCurvature
                description = 'Gaussian Curvature (nm-2)'
            else:
                print(f"Warning: Unknown method '{self.method}'. Defaulting to save Mean Curvature.")
                results_to_save = self.results.MeanCurvature
                description = 'Mean Curvature (nm-1)'

            lipids_ratio = {sp: self.u.select_atoms(f'resname {sp}').n_residues for sp in self.residues}
            dict_parameter = {
                'frames': list(range(self.start, self.stop, self.step)),
                'resids': self.resids,
                'resnames': self.resnames,
                'positions': self.headAtoms.positions,
                'results': results_to_save,
                'file_path': self.file_path,
                'description': description,
                'parameters': self.parameters,
                'lipids_type': lipids_ratio
            }
            WriteExcelLipids(**dict_parameter).run()
            print(f"Analysis complete. Results for '{self.method}' curvature will be saved to {self.file_path}")


if __name__ == "__main__":

    gro_file = "cases/lnb.gro"
    xtc_file = "cases/md.xtc"
    csv_file_mean = "cases/csv/curvature_mean_results.csv"
    csv_file_gauss = "cases/csv/curvature_gauss_results.csv"

    u = mda.Universe(gro_file, xtc_file)
    residue_group = {'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}

    # 1. 串行
    print("--- Running Serial Mean Curvature Analysis ---")
    analysis_serial = Curvature(u, residue_group, k=20, file_path=csv_file_mean, method='mean', parallel=False)
    analysis_serial.run(verbose=True)

    # 2. 并行
    print("--- Running Parallel Gaussian Curvature Analysis ---")
    analysis_parallel = Curvature(u, residue_group, k=20, file_path=csv_file_gauss, method='mean', parallel=True, n_jobs=10)
    analysis_parallel.run(verbose=True)
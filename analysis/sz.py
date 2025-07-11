import warnings
import os
import sys
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

__all__ = ['SZ']


def get_normals(k, points):
    n_points = points.shape[0]
    if n_points < k:
        return np.full_like(points, np.nan)

    kdtree = KDTree(points)
    _, idxs = kdtree.query(points, k=k)
    neighborhoods = points[idxs]

    means = np.mean(neighborhoods, axis=1, keepdims=True)
    centered = neighborhoods - means
    cov_matrices = np.einsum('nki,nkj->nij', centered, centered)

    eigenvalues, eigenvectors = eigh(cov_matrices)
    normals = eigenvectors[:, :, 0]

    geo_center = np.mean(points, axis=0)
    center_vectors = geo_center - points
    dots = np.einsum('ij,ij->i', normals, center_vectors)
    
    flip_mask = dots < 0
    normals[flip_mask] *= -1

    return normals

def ensure_consistent_normals(normals, vectors):
    dots = np.einsum('ij,ij->i', normals, vectors)
    flip_mask = dots < 0
    normals[flip_mask] *= -1
    return normals


class SZ(AnalysisBase):
    """
    A class for calculating the Sz order parameter for lipid acyl chains.
    """
    
    def __init__(self, universe, residues_group: dict, chain: str = 'both', k: int = 15, file_path: str = None,
                 parallel: bool = False, n_jobs: int = -1):
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residues_group)
        self.chain = chain.lower()
        self.k = k
        self.file_path = file_path
        self.parallel = parallel
        self.n_jobs = n_jobs

        self.headAtoms = self.u.atoms[[]]
        self.tailAtoms = self.u.atoms[[]]

        head_sp_map = {sp: ' '.join(residues_group[sp]) for sp in self.residues}
        
        if self.chain == 'sn1':
            tail_selection_str = 'name ??A'
        elif self.chain == 'sn2':
            tail_selection_str = 'name ??B'
        else: # both
            tail_selection_str = 'name ??A ??B'
        
        for sp in self.residues:
            self.headAtoms += self.u.select_atoms(f'resname {sp} and name {head_sp_map[sp]}')
            self.tailAtoms += self.u.select_atoms(f'resname {sp} and ({tail_selection_str})')

        self.n_residues = self.headAtoms.n_residues
        self.resids = self.headAtoms.resids
        self.resnames = self.headAtoms.resnames
        
        self.headMasks = {sp: self.resnames == sp for sp in self.residues}
        self.numSp = {sp: np.sum(mask) for sp, mask in self.headMasks.items()}
        
        self.tail_residue_map = {res.resid: res for res in self.tailAtoms.residues}
        
        # Pre-calculate the static residue ID map to avoid passing the whole universe object
        self.sp_resids_map = {
            sp: np.unique(self.u.select_atoms(f'resname {sp}').residues.resids)
            for sp in self.residues
        }

        self.parameters = f"{residues_group}, Chain:{self.chain}, K:{self.k}"
        self.results.SZ = None

    def _prepare(self):
        self.results.SZ = np.full([self.n_residues, self.n_frames], fill_value=np.nan)

    def _get_tail_pos_map(self):
        tail_pos_map = {}
        for rid, res in self.tail_residue_map.items():
            residue_chains = {}
            sn1_atoms = res.atoms.select_atoms('name ??A')
            sn2_atoms = res.atoms.select_atoms('name ??B')
            if len(sn1_atoms) > 0:
                residue_chains['??A'] = sn1_atoms.positions
            if len(sn2_atoms) > 0:
                residue_chains['??B'] = sn2_atoms.positions
            tail_pos_map[rid] = residue_chains
        return tail_pos_map

    @staticmethod
    def _calculate_sz_for_frame(k, head_positions, tail_residue_pos_map, 
                                residues_list, head_masks, num_sp_map, chain_type, sp_resids_map):
        
        normals = get_normals(k, head_positions)
        sz_results = np.full(head_positions.shape[0], np.nan)

        def compute(tail_positions_sp, head_positions_sp, normals_sp):
            if tail_positions_sp.ndim != 3 or tail_positions_sp.shape[0] == 0:
                return np.full(head_positions_sp.shape[0], np.nan)
            
            tail_end_pos = tail_positions_sp[:, -1, :]
            head_to_tail = head_positions_sp - tail_end_pos
            consistent_normals = ensure_consistent_normals(normals_sp.copy(), head_to_tail)
            
            vectors = tail_positions_sp[:, :-1, :] - tail_positions_sp[:, 1:, :]
            vectors_norm = np.linalg.norm(vectors, axis=2)

            normals_bcast = consistent_normals[:, np.newaxis, :]
            dot_products = np.sum(vectors * normals_bcast, axis=2)
            
            cos_theta = np.divide(dot_products, vectors_norm, 
                                  out=np.zeros_like(dot_products, dtype=float), 
                                  where=vectors_norm!=0)
            
            sz_per_bond = (3 * cos_theta**2 - 1) / 2
            
            return np.mean(sz_per_bond, axis=1)

        for sp in residues_list:
            mask = head_masks[sp]
            num_res_sp = num_sp_map[sp]
            if num_res_sp == 0:
                continue

            sp_normals = normals[mask]
            sp_head_pos = head_positions[mask]
            sp_resids = sp_resids_map[sp]

            sz_sp = np.full(num_res_sp, np.nan)
            if chain_type == 'both':
                sz_sn1_list, sz_sn2_list = [], []
                
                sn1_pos_list = [tail_residue_pos_map.get(rid, {}).get('??A') for rid in sp_resids]
                sn1_pos_list = [pos for pos in sn1_pos_list if pos is not None]

                sn2_pos_list = [tail_residue_pos_map.get(rid, {}).get('??B') for rid in sp_resids]
                sn2_pos_list = [pos for pos in sn2_pos_list if pos is not None]

                if sn1_pos_list:
                    sz_sn1_list = compute(np.array(sn1_pos_list), sp_head_pos, sp_normals)
                if sn2_pos_list:
                    sz_sn2_list = compute(np.array(sn2_pos_list), sp_head_pos, sp_normals)
                
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    sz_sp = np.nanmean([sz_sn1_list, sz_sn2_list], axis=0)
            else:
                chain_name = '??A' if chain_type == 'sn1' else '??B'
                tail_pos_list = [tail_residue_pos_map.get(rid, {}).get(chain_name) for rid in sp_resids]
                tail_pos_list = [pos for pos in tail_pos_list if pos is not None]

                if tail_pos_list:
                    sz_sp = compute(np.array(tail_pos_list), sp_head_pos, sp_normals)
            
            sz_results[mask] = sz_sp
            
        return sz_results

    def _single_frame(self):
        tail_pos_map = self._get_tail_pos_map()
        sz_values = SZ._calculate_sz_for_frame(
            self.k, self.headAtoms.positions, tail_pos_map,
            self.residues, self.headMasks, self.numSp, self.chain, self.sp_resids_map
        )
        self.results.SZ[:, self._frame_index] = sz_values

    def _get_inputs_generator(self):
        for ts in self.u.trajectory[self.start:self.stop:self.step]:
            tail_pos_map = self._get_tail_pos_map()
            yield (self.k, self.headAtoms.positions.copy(), tail_pos_map,
                   self.residues, self.headMasks, self.numSp, self.chain, self.sp_resids_map)

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
                delayed(SZ._calculate_sz_for_frame)(*inputs) for inputs in inputs_generator
            )
            if results_list:
                self.results.SZ = np.array(results_list).T
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose)

        self._conclude()

    def _conclude(self):
        if self.file_path:
            lipids_ratio = {sp: self.u.select_atoms(f'resname {sp}').n_residues for sp in self.residues}
            dict_parameter = {
                'frames': list(range(self.start, self.stop, self.step)),
                'resids': self.resids,
                'resnames': self.resnames,
                'positions': self.headAtoms.positions,
                'results': self.results.SZ,
                'file_path': self.file_path,
                'description': f'Sz Order Parameter (chain: {self.chain})',
                'parameters': self.parameters,
                'lipids_type': lipids_ratio
            }
            # Assuming WriteExcelLipids is a defined class and available
            WriteExcelLipids(**dict_parameter).run()
            print(f"Analysis complete. Results will be saved to {self.file_path}")


if __name__ == "__main__":
    import time
    gro_file = "cases/lnb.gro"
    xtc_file = "cases/md.xtc"
    csv_file_serial = "cases/csv/sz_serial.csv"
    csv_file_parallel = "cases/csv/sz_parallel.csv"

    u = mda.Universe(gro_file, xtc_file)
    residues_group = {'DPPC': ['PO4'], 'DUPC': ['PO4']}
    #1. 串行
    t1 = time.time()
    analysis_serial = SZ(u, residues_group, chain='sn1', k=15, file_path=csv_file_serial, parallel=False)
    analysis_serial.run(verbose=True)
    t2 = time.time()
    print(f"Serial analysis completed in {t2 - t1:.2f} seconds.")
    
    #2. 并行
    t3 = time.time()
    analysis_parallel = SZ(u, residues_group, chain='sn1', k=15, file_path=csv_file_parallel, parallel=True, n_jobs=2)
    analysis_parallel.run(verbose=True)
    t4 = time.time()
    print(f"Parallel analysis completed in {t4 - t3:.2f} seconds.")
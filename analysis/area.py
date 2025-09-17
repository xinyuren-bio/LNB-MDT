import warnings
warnings.filterwarnings('ignore')
import sys
import os
import argparse
import ast     

import numpy as np
from scipy.spatial import Voronoi
from typing import List, Tuple
import math
from joblib import Parallel, delayed

import MDAnalysis as mda
from scipy.spatial import KDTree
from scipy.linalg import eigh

if __name__ == '__main__':
    current_file_path = os.path.abspath(__file__)
    current_dir = os.path.dirname(current_file_path)
    package_root = os.path.abspath(os.path.join(current_dir, '..'))
    if package_root not in sys.path:
        sys.path.insert(0, package_root)
from analysis.analysis_base import *

__all__ = ['Area']


class Point:
    def __init__(self, x: float, y: float):
        self.x = x
        self.y = y

class Polygon:
    def __init__(self):
        self.points: List[Point] = []
    
    def append(self, point: Point):
        self.points.append(point)
    
    def empty(self):
        self.points.clear()
    
    @property
    def size(self) -> int:
        return len(self.points)

def same_side(p1: Point, p2: Point, line_p1: Point, line_p2: Point) -> bool:
    def cross_product(p1: Point, p2: Point, p3: Point) -> float:
        return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x)
    
    cp1 = cross_product(line_p1, line_p2, p1)
    cp2 = cross_product(line_p1, line_p2, p2)
    return (cp1 * cp2) >= 0

def are_parallel(p1: Point, p2: Point, p3: Point, p4: Point) -> bool:
    def get_direction(p1: Point, p2: Point) -> Tuple[float, float]:
        return (p2.x - p1.x, p2.y - p1.y)
    
    dir1 = get_direction(p1, p2)
    dir2 = get_direction(p3, p4)
    
    return abs(dir1[0] * dir2[1] - dir1[1] * dir2[0]) < 1e-10

def get_intersection(p1: Point, p2: Point, p3: Point, p4: Point) -> Point:
    def det(a: float, b: float, c: float, d: float) -> float:
        return a * d - b * c
    
    x1, y1 = p1.x, p1.y
    x2, y2 = p2.x, p2.y
    x3, y3 = p3.x, p3.y
    x4, y4 = p4.x, p4.y
    
    denominator = det(x1 - x2, y1 - y2, x3 - x4, y3 - y4)
    if abs(denominator) < 1e-10:
        return Point((x1 + x2) / 2, (y1 + y2) / 2)
    
    t = det(x1 - x3, y1 - y3, x3 - x4, y3 - y4) / denominator
    x = x1 + t * (x2 - x1)
    y = y1 + t * (y2 - y1)
    
    return Point(x, y)

def fast_clip_zoi(zoi: Polygon, ref_pt: Point, clipping_pt: Point, buffer: Polygon = None) -> Polygon:
    if buffer is None:
        buffer = Polygon()
    
    middle_pt = Point(
        0.5 * (clipping_pt.x + ref_pt.x),
        0.5 * (clipping_pt.y + ref_pt.y)
    )
    
    delta_x = ref_pt.x - clipping_pt.x
    delta_y = ref_pt.y - clipping_pt.y
    
    line_dir_x = delta_y
    line_dir_y = -delta_x
    
    other_pt = Point(
        middle_pt.x + line_dir_x,
        middle_pt.y + line_dir_y
    )
    
    buffer.empty()

    for second_vid in range(zoi.size):
        first_vid = zoi.size - 1 if second_vid == 0 else second_vid - 1
        
        second_same_side = same_side(ref_pt, zoi.points[second_vid], middle_pt, other_pt)
        
        if same_side(zoi.points[first_vid], zoi.points[second_vid], middle_pt, other_pt):
            if not second_same_side:
                continue
        else:
            if not are_parallel(zoi.points[first_vid], zoi.points[second_vid], middle_pt, other_pt):
                inter_pt = get_intersection(zoi.points[first_vid], zoi.points[second_vid], middle_pt, other_pt)
                buffer.append(inter_pt)
        
        if second_same_side:
            buffer.append(zoi.points[second_vid])
    
    return buffer

def compute_voronoi_area(i, points):
    """计算Voronoi面积"""
    try:
        vor = Voronoi(points)
        region_index = vor.point_region[0] 
        region_vertices = vor.regions[region_index]
        
        if -1 in region_vertices: 
            # print(f"点 {i} 的Voronoi区域是无限的，尝试裁剪...")

            zoi = Polygon()
            for vertex in vor.vertices:
                zoi.append(Point(vertex[0], vertex[1]))

            center = Point(points[0, 0], points[0, 1])
            max_dist = 0
            for point in points[1:]:
                dist = math.sqrt((point[0] - center.x)**2 + (point[1] - center.y)**2)
                max_dist = max(max_dist, dist)

            clip_dist = max_dist * 2
            clip_points = [
                Point(center.x + clip_dist, center.y),
                Point(center.x - clip_dist, center.y),
                Point(center.x, center.y + clip_dist),
                Point(center.x, center.y - clip_dist)
            ]
            
            clipped_polygon = zoi
            for clip_pt in clip_points:
                clipped_polygon = fast_clip_zoi(clipped_polygon, center, clip_pt)
            
            if clipped_polygon.size >= 3:
                vertices = np.array([[p.x, p.y] for p in clipped_polygon.points])
                area = polygon_area(vertices)
                return area
            else:
                print(f"点 {i} 裁剪后的多边形点数不足")
                return np.nan
        else:
            finite_vertices = vor.vertices[region_vertices]
            area = polygon_area(finite_vertices)
            
            if area <= 0:
                print(f"点 {i} 的Voronoi面积为负或零: {area}")
                return np.nan
            
            return area
            
    except Exception as e:
        print(f"点 {i} 的Voronoi计算出错: {str(e)}")
        return np.nan

def polygon_area(vertices):
    x = vertices[:, 0]
    y = vertices[:, 1]
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

class Area(AnalysisBase):
    def __init__(self, universe, residueGroup: dict, k: int = None, filePath: str = None, 
                 max_normal_angle_deg: float = 140, parallel: bool = False, n_jobs: int = -1):
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residueGroup)
        self.k = k if k is not None else 15
        self.file_path = filePath
        self.max_normal_angle_deg = max_normal_angle_deg
        self.parallel = parallel
        self.n_jobs = n_jobs

        self.headSp = {sp: ' '.join(residueGroup[sp]) for sp in residueGroup}
        print("Head atoms:", self.headSp)

        self.headAtoms = self.u.atoms[[]]
        for i in range(len(self.residues)):
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (self.residues[i], self.headSp[self.residues[i]]), updating=False)

        self._n_residues = self.headAtoms.n_residues
        self.resids = self.headAtoms.resids
        self.resnames = self.headAtoms.resnames
        self.results.Area = None

        self.parameters = str(residueGroup) + f'K:{self.k}'
        if self.max_normal_angle_deg is not None:
            self.parameters += f', MaxAngle:{self.max_normal_angle_deg}'

    @property
    def Area(self):
        return self.results.Area

    def _prepare(self):
        self.results.Area = np.full([self._n_residues, self.n_frames],
                                    fill_value=np.nan)

    @staticmethod
    def _calculate_area_for_frame(headPos, k_val, max_normal_angle_deg, resnames):
        n_points = len(headPos)
        
        if n_points == 0:
            return np.full(len(resnames), np.nan)
        kd_tree = KDTree(headPos)
        actual_k_for_query = min(k_val, n_points)
        if actual_k_for_query == 0:
             return np.full(len(resnames), np.nan)

        _, initial_all_idxs = kd_tree.query(headPos, k=actual_k_for_query)
        initial_expanded_neighbors = headPos[initial_all_idxs]

        initial_means = np.mean(initial_expanded_neighbors, axis=1, keepdims=True)
        initial_centered_points = initial_expanded_neighbors - initial_means
        initial_cov_matrices = np.einsum('nki,nkj->nij', initial_centered_points, initial_centered_points) / actual_k_for_query
        _, initial_eigenvectors = np.linalg.eigh(initial_cov_matrices)
        all_particle_normals = initial_eigenvectors[:, :, 0]

        # ================== 根据法向量夹角筛选邻居 ==================
        filtered_neighbor_positions_list = []
        max_filtered_neighbors_count = 0

        for i in range(n_points):
            central_atom_normal = all_particle_normals[i]
            true_neighbor_indices = initial_all_idxs[i, 1:]
            
            kept_true_neighbor_indices = []
            if max_normal_angle_deg is not None and len(true_neighbor_indices) > 0 :
                true_neighbor_normals = all_particle_normals[true_neighbor_indices]
                dot_products = np.einsum('j,kj->k', central_atom_normal, true_neighbor_normals)
                dot_products_clipped = np.clip(dot_products, -1.0, 1.0)
                angles_rad = np.arccos(dot_products_clipped)
                angles_deg = np.degrees(angles_rad)
                
                pass_filter_mask = angles_deg <= max_normal_angle_deg
                kept_true_neighbor_indices = true_neighbor_indices[pass_filter_mask]
            elif len(true_neighbor_indices) > 0:
                kept_true_neighbor_indices = true_neighbor_indices
            
            unique_final_indices = [i]
            for idx_val in kept_true_neighbor_indices:
                if idx_val != i and idx_val not in unique_final_indices:
                     unique_final_indices.append(idx_val)
            
            current_positions = headPos[unique_final_indices]
            filtered_neighbor_positions_list.append(current_positions)
            if len(current_positions) > max_filtered_neighbors_count:
                max_filtered_neighbors_count = len(current_positions)

        if max_filtered_neighbors_count == 0 and n_points > 0:
            expanded_neighbors_final = np.full((n_points, 1, 3), np.nan)
            for i in range(n_points):
                expanded_neighbors_final[i,0,:] = headPos[i]
            max_filtered_neighbors_count = 1
        elif n_points == 0:
            expanded_neighbors_final = np.empty((0,0,3))
        else:
            expanded_neighbors_final = np.full((n_points, max_filtered_neighbors_count, 3), np.nan)
            for i in range(n_points):
                data = filtered_neighbor_positions_list[i]
                if len(data) > 0:
                    expanded_neighbors_final[i, :len(data)] = data
        means = np.nanmean(expanded_neighbors_final, axis=1, keepdims=True)
        centered_points = expanded_neighbors_final - means
        mask = ~np.isnan(centered_points).any(axis=2)
        counts_for_pca = mask.sum(axis=1).astype(float)
        # print(counts_for_pca)   
        cov_matrices = np.zeros((n_points, 3, 3))
        for i in range(n_points):
            if counts_for_pca[i] >= 1:
                active_centered_points = centered_points[i, mask[i], :]
                if active_centered_points.shape[0] > 0:
                    cov_matrices[i] = np.einsum('ki,kj->ij', active_centered_points, active_centered_points) / active_centered_points.shape[0]

        _, eigenvectors_for_projection = np.linalg.eigh(cov_matrices)

        x_axes = eigenvectors_for_projection[:, :, 2]
        y_axes = eigenvectors_for_projection[:, :, 1]

        reference_points = headPos[:, np.newaxis, :]
        relative_points = expanded_neighbors_final - reference_points
        local_x = np.einsum('nkj,nj->nk', relative_points, x_axes)
        local_y = np.einsum('nkj,nj->nk', relative_points, y_axes)
        local_coords_2d = np.stack((local_x, local_y), axis=-1)

        # ================== 批量计算Voronoi面积 ==================
        def compute_voronoi_area_single(i):
            points = local_coords_2d[i]
            if np.isnan(points[0]).any():
                print(f"点 {i} 的中心点在投影后有NaN")
                return np.nan
            try:
                vor = Voronoi(points)
                region_index = vor.point_region[0]
                region_vertices = vor.regions[region_index]
                if -1 in region_vertices:
                    print(f"点 {i} 的Voronoi区域不是有限的")
                    return np.nan
                finite_vertices = vor.vertices[region_vertices]
                area = polygon_area(finite_vertices)
                if area <= 0:
                    print(f"点 {i} 的Voronoi面积为负或零: {area}")
                    return np.nan
                return area
            except Exception as e:
                print(f"点 {i} 的Voronoi计算出错: {str(e)}")
                return np.nan

        area_arr = np.array([compute_voronoi_area_single(i) for i in range(n_points)])
        
        nan_count = np.sum(np.isnan(area_arr))
        # print(f"\n面积计算统计:")
        # print(f"总点数: {n_points}")
        # print(f"NaN面积数: {nan_count}")
        # if n_points > 0:
            # print(f"NaN比例: {nan_count/n_points*100:.2f}%")

        filled_area_arr = np.copy(area_arr)
        filled_count = 0
        if n_points > 0:
            global_mean_area_this_frame = np.nanmean(area_arr)
            if np.isnan(global_mean_area_this_frame):
                global_mean_area_this_frame = 0.0

            unique_resnames_in_frame = np.unique(resnames)

            for res_type in unique_resnames_in_frame:
                type_mask = (resnames == res_type)
                mean_area_for_this_type = np.nanmean(area_arr[type_mask])
                
                fill_value_for_type = mean_area_for_this_type if not np.isnan(mean_area_for_this_type) else global_mean_area_this_frame
                
                nan_mask_for_this_type = type_mask & np.isnan(filled_area_arr)
                filled_area_arr[nan_mask_for_this_type] = fill_value_for_type
                filled_count += np.sum(nan_mask_for_this_type)
        
        if np.any(np.isnan(filled_area_arr)):
             final_fallback_mean = np.nanmean(filled_area_arr)
             if np.isnan(final_fallback_mean): final_fallback_mean = 0.0
             nan_mask = np.isnan(filled_area_arr)
             filled_area_arr[nan_mask] = final_fallback_mean
             filled_count += np.sum(nan_mask)
        
        return filled_area_arr * 0.01

    def _single_frame(self):
        """
        Processes a single frame in serial mode.
        """
        headPos = self.headAtoms.positions
        area_arr = self._calculate_area_for_frame(
            headPos, self.k, self.max_normal_angle_deg, self.resnames
        )
        self.results.Area[:, self._frame_index] = area_arr
    
    def _get_inputs_generator(self):
        """
        Generator that yields the necessary inputs for each frame's calculation.
        """
        for ts in self.u.trajectory[self.start:self.stop:self.step]:
            yield (self.headAtoms.positions.copy(), self.k, self.max_normal_angle_deg, self.resnames)

    def run(self, start=None, stop=None, step=None, verbose=None, callBack=None):
        """
        Execute the analysis.
        Chooses between serial and parallel execution based on self.parallel.
        """
        if self.parallel:
            self.start = start if start is not None else 0
            self.stop = stop if stop is not None and stop < self._trajectory.n_frames else self._trajectory.n_frames
            self.step = step if step is not None else 1
            self.n_frames = len(range(self.start, self.stop, self.step))
            
            self._prepare()

            print(f"Running in parallel on {self.n_jobs} jobs...")
            verbose_level = 10 if verbose else 0
            
            inputs_generator = self._get_inputs_generator()
            
            results_list = Parallel(n_jobs=self.n_jobs, verbose=verbose_level)(
                delayed(Area._calculate_area_for_frame)(*inputs) for inputs in inputs_generator
            )
            
            if results_list:
                results_array = np.array(results_list)
                self.results.Area = results_array.T
            
            self._conclude()
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose, callBack=callBack)

    def _conclude(self):
        if self.file_path:
            lipids_ratio = {sp: self.u.select_atoms(f'resname {sp}').n_residues for sp in self.residues}
            dict_parameter = {
                'frames': [i for i in range(self.start, self.stop, self.step)]
                , 'resids': self.resids
                , 'resnames': self.resnames
                , 'positions': self.headAtoms.positions
                , 'results': self.results.Area
                , 'file_path': self.file_path
                , 'description': 'Area(nm^2)'
                , 'parameters': self.parameters
                , 'lipids_type': lipids_ratio
            }
            WriteExcelLipids(**dict_parameter).run()
            print(f"Analysis complete. Results saved to {self.file_path}")


# --- Command-line Argument Parsing ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform Area analysis on molecular dynamics trajectories."
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
        default="cases/csv/area_results.csv",
        help="Path to the output CSV file for area results."
    )
    parser.add_argument(
        "--residues", "-r",
        type=str,
        default="{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}",
        help="A dictionary string defining residue groups for analysis. E.g., \"{'DPPC': ['PO4'], 'CHOL': ['ROH']}\""
    )
    parser.add_argument(
        "--k-value", "-k",
        type=int,
        default=20,
        help="K value for Voronoi tessellation."
    )
    parser.add_argument(
        "--max-normal-angle",
        type=float,
        default=140,
        help="Maximum normal angle in degrees."
    )
    parser.add_argument(
        "--parallel", "-p",
        action="store_true",
        help="Enable parallel processing for area calculation."
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
        print("Please ensure it's a valid dictionary string, e.g., \"{'DPPC': ['PO4'], 'CHOL': ['ROH']}\"")
        sys.exit(1)

    print("\n--- Initializing MDAnalysis Universe ---")
    try:
        u = mda.Universe(args.gro_file, args.xtc_file)
    except Exception as e:
        print(f"Error loading MDAnalysis Universe: {e}")
        print("Please check if GRO/XTC files exist and are valid.")
        sys.exit(1)

    print("\n--- Running Area Analysis ---")
    area_analysis = Area(
        u,
        residues_group_parsed,
        k=args.k_value,
        max_normal_angle_deg=args.max_normal_angle,
        file_path=args.output_csv,
        parallel=args.parallel,
        n_jobs=args.n_jobs
    )
    area_analysis.run(
        start=args.start_frame,
        stop=args.stop_frame,
        step=args.step_frame,
        verbose=args.verbose
    )

    print("\n--- Analysis Finished ---")
    print(f"Results saved to: {args.output_csv}")

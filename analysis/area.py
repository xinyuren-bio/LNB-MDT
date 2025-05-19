import warnings
warnings.filterwarnings('ignore')
import sys
import os
import numpy as np
from scipy.spatial import Voronoi
from typing import List, Tuple
import math

# 获取当前脚本文件的绝对路径
script_path = os.path.abspath(__file__)
# 获取脚本所在的目录
script_dir = os.path.dirname(script_path)
# 计算项目根目录 (假设脚本在 analysis/test_function/ 下，根目录是上两级)
project_root = os.path.abspath(os.path.join(script_dir, '..', '..'))

# 如果项目根目录不在 sys.path 中，则添加到列表的开头
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import MDAnalysis as mda
from scipy.spatial import KDTree
from scipy.linalg import eigh

from analysis.analysis_base import *

# suppress some MDAnalysis warnings when writing PDB/GRO files
__all__ = ['Area']

# 定义点和多边形结构
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
    """检查两个点是否在直线的同一侧"""
    def cross_product(p1: Point, p2: Point, p3: Point) -> float:
        return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x)
    
    cp1 = cross_product(line_p1, line_p2, p1)
    cp2 = cross_product(line_p1, line_p2, p2)
    return (cp1 * cp2) >= 0

def are_parallel(p1: Point, p2: Point, p3: Point, p4: Point) -> bool:
    """检查两条线是否平行"""
    def get_direction(p1: Point, p2: Point) -> Tuple[float, float]:
        return (p2.x - p1.x, p2.y - p1.y)
    
    dir1 = get_direction(p1, p2)
    dir2 = get_direction(p3, p4)
    
    # 检查方向向量是否平行
    return abs(dir1[0] * dir2[1] - dir1[1] * dir2[0]) < 1e-10

def get_intersection(p1: Point, p2: Point, p3: Point, p4: Point) -> Point:
    """计算两条线的交点"""
    def det(a: float, b: float, c: float, d: float) -> float:
        return a * d - b * c
    
    x1, y1 = p1.x, p1.y
    x2, y2 = p2.x, p2.y
    x3, y3 = p3.x, p3.y
    x4, y4 = p4.x, p4.y
    
    denominator = det(x1 - x2, y1 - y2, x3 - x4, y3 - y4)
    if abs(denominator) < 1e-10:
        return Point((x1 + x2) / 2, (y1 + y2) / 2)  # 如果平行，返回中点
    
    t = det(x1 - x3, y1 - y3, x3 - x4, y3 - y4) / denominator
    x = x1 + t * (x2 - x1)
    y = y1 + t * (y2 - y1)
    
    return Point(x, y)

def fast_clip_zoi(zoi: Polygon, ref_pt: Point, clipping_pt: Point, buffer: Polygon = None) -> Polygon:
    """快速裁剪多边形区域"""
    if buffer is None:
        buffer = Polygon()
    
    # 计算裁剪线
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
    
    # 检查线是否与多边形相交
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
        region_index = vor.point_region[0]  # 第一个点是中心点
        region_vertices = vor.regions[region_index]
        
        if -1 in region_vertices:  # 如果是无限区域
            print(f"点 {i} 的Voronoi区域是无限的，尝试裁剪...")
            
            # 创建多边形对象
            zoi = Polygon()
            for vertex in vor.vertices:
                zoi.append(Point(vertex[0], vertex[1]))
            
            # 计算裁剪点
            center = Point(points[0, 0], points[0, 1])
            max_dist = 0
            for point in points[1:]:
                dist = math.sqrt((point[0] - center.x)**2 + (point[1] - center.y)**2)
                max_dist = max(max_dist, dist)
            
            # 在远处添加裁剪点
            clip_dist = max_dist * 2
            clip_points = [
                Point(center.x + clip_dist, center.y),
                Point(center.x - clip_dist, center.y),
                Point(center.x, center.y + clip_dist),
                Point(center.x, center.y - clip_dist)
            ]
            
            # 对每个裁剪点进行裁剪
            clipped_polygon = zoi
            for clip_pt in clip_points:
                clipped_polygon = fast_clip_zoi(clipped_polygon, center, clip_pt)
            
            # 计算裁剪后的面积
            if clipped_polygon.size >= 3:
                vertices = np.array([[p.x, p.y] for p in clipped_polygon.points])
                area = polygon_area(vertices)
                return area
            else:
                print(f"点 {i} 裁剪后的多边形点数不足")
                return np.nan
        else:
            # 处理有限区域
            finite_vertices = vor.vertices[region_vertices]
            area = polygon_area(finite_vertices)
            
            if area <= 0:
                print(f"点 {i} 的Voronoi面积为负或零: {area}")
                return np.nan
            
            return area
            
    except Exception as e:
        print(f"点 {i} 的Voronoi计算出错: {str(e)}")
        return np.nan

class Area(AnalysisBase):
    """
    A class for calculating the area per lipid in a bilayer system.
    
    This class analyzes the surface area occupied by each lipid molecule
    in a molecular dynamics trajectory. It uses Voronoi tessellation
    to determine the area occupied by each lipid's head group.
    """
    
    def __init__(self, universe, residueGroup: dict, k: int = None, file_path: str = None, max_normal_angle_deg: float = 140):
        """
        Initialize the Area analysis class.
        
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
            The number of nearest neighbors to use for calculating local normal vectors
            and for defining the neighborhood for area calculation.
            If None, a default value will be used (e.g., 15).
            
        file_path : str, optional
            The path where the analysis results will be saved as a CSV file.
            If None, results will not be saved to disk.

        max_normal_angle_deg : float, optional
            The maximum angle (in degrees) between a central atom's normal vector
            and its neighbor's normal vector for the neighbor to be included in
            the projection and Voronoi calculation. If None, no angle filtering is applied.
            A typical value might be 60 or 90 degrees.
            
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
            
        results.Area : numpy.ndarray
            A 2D array storing the area values for all lipid molecules across all frames.
            Shape: (n_residues, n_frames)
        """
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residueGroup)
        self.k = k if k is not None else 15 # Default k if not provided
        self.file_path = file_path
        self.max_normal_angle_deg = max_normal_angle_deg

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
        self.results.Area = None

        # Record analysis parameters
        self.parameters = str(residueGroup) + f'K:{self.k}'
        if self.max_normal_angle_deg is not None:
            self.parameters += f', MaxAngle:{self.max_normal_angle_deg}'

    @property
    def Area(self):
        return self.results.Area

    def _prepare(self):
        self.results.Area = np.full([self._n_residues, self.n_frames],
                                    fill_value=np.nan)

    def _single_frame(self):
        headPos = self.headAtoms.positions
        n_points = len(headPos)
        
        if n_points == 0:
            self.results.Area[:, self._frame_index] = np.nan
            return

        k_val = self.k

        # ================== 初始邻近点查询和法向量计算 ==================
        kd_tree = KDTree(headPos)
        # Query k+1 because the first neighbor is the point itself. We need k true neighbors for a stable normal.
        # If k_val is small (e.g. <3 for PCA), PCA might be unstable. Ensure k_val is reasonable for PCA.
        # Let's use self.k for normal calculation consistently.
        actual_k_for_query = min(k_val, n_points) # Cannot query more neighbors than points available
        if actual_k_for_query == 0 : # if k=0 or n_points=0 (already handled) or n_points=1 and k=1
             self.results.Area[:, self._frame_index] = np.nan # or some default for single/no points
             return


        _, initial_all_idxs = kd_tree.query(headPos, k=actual_k_for_query)
        
        initial_expanded_neighbors = headPos[initial_all_idxs]

        # 计算每个粒子基于初始k邻居的法向量
        initial_means = np.mean(initial_expanded_neighbors, axis=1, keepdims=True)
        initial_centered_points = initial_expanded_neighbors - initial_means
        initial_cov_matrices = np.einsum('nki,nkj->nij', initial_centered_points, initial_centered_points) / actual_k_for_query
        _, initial_eigenvectors = np.linalg.eigh(initial_cov_matrices)
        all_particle_normals = initial_eigenvectors[:, :, 0]  # (n_points, 3) 每个粒子的法向量

        # ================== 根据法向量夹角筛选邻居 ==================
        filtered_neighbor_positions_list = []
        max_filtered_neighbors_count = 0

        for i in range(n_points):
            central_atom_idx = initial_all_idxs[i, 0] # Should be 'i' if n_points > 0
            # Ensure central_atom_idx is indeed i, KDTree.query guarantees this for k>=1
            if central_atom_idx != i and n_points > 1 : # Safety check, though unlikely with KDTree
                 # This case implies headPos[i] was not its own closest neighbor, problematic.
                 # For simplicity, we proceed assuming initial_all_idxs[i,0] == i
                 pass

            central_atom_normal = all_particle_normals[i]
            
            # True neighbors are initial_all_idxs[i, 1:]
            true_neighbor_indices = initial_all_idxs[i, 1:]
            
            kept_true_neighbor_indices = []
            if self.max_normal_angle_deg is not None and len(true_neighbor_indices) > 0 :
                true_neighbor_normals = all_particle_normals[true_neighbor_indices]
                dot_products = np.einsum('j,kj->k', central_atom_normal, true_neighbor_normals)
                # Clip dot product for arccos stability
                dot_products_clipped = np.clip(dot_products, -1.0, 1.0)
                angles_rad = np.arccos(dot_products_clipped)
                angles_deg = np.degrees(angles_rad)
                
                pass_filter_mask = angles_deg <= self.max_normal_angle_deg
                kept_true_neighbor_indices = true_neighbor_indices[pass_filter_mask]
            elif len(true_neighbor_indices) > 0: # No angle filtering, keep all true neighbors
                kept_true_neighbor_indices = true_neighbor_indices
            
            # Final list for projection: central atom + filtered true neighbors
            # The central atom's index is 'i' (or initial_all_idxs[i,0])
            final_indices_for_this_atom = np.concatenate(([i], kept_true_neighbor_indices)).astype(int)
            
            # Remove potential duplicates if 'i' was somehow in kept_true_neighbor_indices (shouldn't be by design)
            # Using unique while preserving order of first element 'i'
            unique_final_indices = [i]
            for idx_val in kept_true_neighbor_indices:
                if idx_val != i: # True neighbors should not be the central atom itself
                    if idx_val not in unique_final_indices: # Avoid duplicates among neighbors
                         unique_final_indices.append(idx_val)
            
            current_positions = headPos[unique_final_indices]
            filtered_neighbor_positions_list.append(current_positions)
            if len(current_positions) > max_filtered_neighbors_count:
                max_filtered_neighbors_count = len(current_positions)

        if max_filtered_neighbors_count == 0 and n_points > 0: # Should not happen if central atom i is always included.
            # This means even the central atom wasn't included, or n_points = 0
            # If it's 1 (only central atom), Voronoi will handle it (yield NaN).
            # If truly 0, create a dummy with at least the central atom
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
                if len(data) > 0: # Ensure data is not empty before assignment
                    expanded_neighbors_final[i, :len(data)] = data
                # If len(data) is 0, it remains all NaNs for this atom, which is fine.

        # ================== 批量计算局部坐标系 (基于筛选后的邻居) ==================
        # Check if any row in expanded_neighbors_final is all NaNs (e.g. if a point had no valid neighbors including self)
        # This can happen if max_filtered_neighbors_count was 0 and then forced to 1 with only headPos[i]
        # or if filtered_neighbor_positions_list[i] was empty.
        
        # Ensure there are enough points for PCA in each neighborhood
        # If a row in expanded_neighbors_final has < 3 actual points (non-NaN rows in the k-dimension), PCA is problematic.
        # The nanmean and subsequent steps should handle this, but resulting x_axes, y_axes might be non-ideal.
        
        means = np.nanmean(expanded_neighbors_final, axis=1, keepdims=True)
        centered_points = expanded_neighbors_final - means

        mask = ~np.isnan(centered_points).any(axis=2) # Mask for valid points (k-dimension)
        
        # Handle cases where all neighbors for a point (including itself) were NaN, making `means` NaN for that point.
        # This would make `centered_points` all NaN for that point.
        means_are_nan = np.isnan(means).all(axis=(1,2)) # Check if mean calculation resulted in NaN for any point i
        if np.any(means_are_nan):
            # For points where mean is NaN, centered_points will be all NaN. Mask will be all False. Counts will be 0.
            # This path needs to result in NaN area for such points.
            pass

        # valid_points_for_pca = np.where(mask[:, :, np.newaxis], centered_points, 0)
        
        counts_for_pca = mask.sum(axis=1).astype(float) # Number of valid neighbors for each central point
        print(counts_for_pca)   
        cov_matrices = np.zeros((n_points, 3, 3))
        for i in range(n_points):
            if counts_for_pca[i] >= 1: # Ensure there's at least one point
                # Get the actual centered points for this atom that are not NaN
                active_centered_points = centered_points[i, mask[i], :] # Shape (counts_for_pca[i], 3)
                
                # Sum of outer products of centered vectors, divided by the number of vectors
                if active_centered_points.shape[0] > 0: # Should be true if counts_for_pca[i] >= 1
                    cov_matrices[i] = np.einsum('ki,kj->ij', active_centered_points, active_centered_points) / active_centered_points.shape[0]
                else: # Fallback, though counts_for_pca[i] >= 1 should prevent this path
                    cov_matrices[i] = np.zeros((3,3))
            else: # No valid points at all for point i (e.g. only NaNs in expanded_neighbors_final for this point)
                cov_matrices[i] = np.zeros((3,3))


        # Handle cases where cov_matrices might be singular or ill-conditioned for some points
        # e.g. if counts_for_pca[i] < 3.
        # _, eigenvectors_for_projection = np.linalg.eigh(cov_matrices) # This might fail if matrices are bad
        eigenvectors_for_projection = np.zeros_like(cov_matrices)
        valid_pca_computation = counts_for_pca >=1 # or counts_for_pca >=3 for meaningful PCA
        
        # Attempt eigh only for valid covariance matrices
        # For invalid ones, eigenvectors can be set to identity or kept as zeros, leading to defined x/y axes.
        # If cov_matrices[i] is all zeros, eigh returns identity matrix as eigenvectors.
        _, eigenvectors_for_projection = np.linalg.eigh(cov_matrices)


        x_axes = eigenvectors_for_projection[:, :, 2]
        y_axes = eigenvectors_for_projection[:, :, 1]
        # normals_of_plane = eigenvectors_for_projection[:, :, 0] # Normals of the projection planes

        # 将点相对于被分析粒子进行平移
        reference_points = headPos[:, np.newaxis, :]
        relative_points = expanded_neighbors_final - reference_points # expanded_neighbors_final[i,0,:] is headPos[i]

        # 转换到局部二维坐标系
        local_x = np.einsum('nkj,nj->nk', relative_points, x_axes)
        local_y = np.einsum('nkj,nj->nk', relative_points, y_axes)
        local_coords_2d = np.stack((local_x, local_y), axis=-1)

        # ================== 批量计算Voronoi面积 ==================
        def compute_voronoi_area(i):
            points = local_coords_2d[i]  # (k, 2)
            # 只检查中心点（第一个点）是否有NaN
            if np.isnan(points[0]).any():
                print(f"点 {i} 的中心点在投影后有NaN")
                return np.nan
                
            try:
                vor = Voronoi(points)
                region_index = vor.point_region[0]  # 第一个点是中心点
                region_vertices = vor.regions[region_index]
                if -1 in region_vertices:  # 检查是否为有限区域
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

        # 逐点计算Voronoi面积
        area_arr = np.array([compute_voronoi_area(i) for i in range(n_points)])
        
        # 打印NaN统计信息
        nan_count = np.sum(np.isnan(area_arr))
        print(f"\n面积计算统计:")
        print(f"总点数: {n_points}")
        print(f"NaN面积数: {nan_count}")
        print(f"NaN比例: {nan_count/n_points*100:.2f}%")

        # 根据不同resname替换NaN值为各自类型的平均值
        filled_area_arr = np.copy(area_arr)
        filled_count = 0  # 记录使用平均值填充的数量
        if n_points > 0: # Ensure there are points to process
            # 计算全局帧内平均面积作为备用
            global_mean_area_this_frame = np.nanmean(area_arr)
            if np.isnan(global_mean_area_this_frame): # 如果所有面积都是NaN
                global_mean_area_this_frame = 0.0 # 或其他合理的默认值

            unique_resnames_in_frame = np.unique(self.resnames) # self.resnames is (n_residues,)

            for res_type in unique_resnames_in_frame:
                type_mask = (self.resnames == res_type)
                mean_area_for_this_type = np.nanmean(area_arr[type_mask])
                
                fill_value_for_type = mean_area_for_this_type
                if np.isnan(fill_value_for_type):
                    fill_value_for_type = global_mean_area_this_frame # Fallback
                
                # 定位当前类型中需要填充的NaN
                nan_mask_for_this_type = type_mask & np.isnan(filled_area_arr) # Check NaNs in the array we are filling
                filled_area_arr[nan_mask_for_this_type] = fill_value_for_type
                filled_count += np.sum(nan_mask_for_this_type)  # 累加填充的数量
        
        # 如果在类型特定填充后仍有NaN（理论上不应该，除非global_mean_area_this_frame也是NaN且未设默认值）
        if np.any(np.isnan(filled_area_arr)):
             final_fallback_mean = np.nanmean(filled_area_arr) # Re-calculate or use global_mean_area_this_frame
             if np.isnan(final_fallback_mean):
                 final_fallback_mean = 0.0 # Ultimate fallback
             nan_mask = np.isnan(filled_area_arr)
             filled_area_arr[nan_mask] = final_fallback_mean
             filled_count += np.sum(nan_mask)  # 累加最终填充的数量

        # 打印填充统计信息
        total_residues = len(area_arr)
        filled_ratio = filled_count / total_residues * 100
        print(f'填充统计信息:')
        print(f'总残基数: {total_residues}')
        print(f'使用平均值填充的残基数: {filled_count}')
        print(f'填充比例: {filled_ratio:.2f}%')
        
        # 打印替换后的标准差
        print('std AFTER type-specific filling:', np.std(filled_area_arr * 0.01))
        
        # 存储结果
        self.results.Area[:, self._frame_index] = filled_area_arr * 0.01

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


def polygon_area(vertices):
    # 计算多边形面积（顺时针或逆时针）
    x = vertices[:, 0]
    y = vertices[:, 1]
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


if __name__ == "__main__":
    import time
    import os
    import matplotlib.pyplot as plt
    print(f"Current working directory: {os.getcwd()}")

    gro_file = "cases/ach.gro"
    # csv_file = "cases/csv/area_test.csv" # 绘图时暂时不保存CSV，以免生成过多文件

    u = mda.Universe(gro_file)
    residue_group = {'DPPC':['PO4'], 'DUPC':['PO4'], 'CHOL':['ROH']}
    k_value = 21

    angles_to_test = [20, 60.0, 140, 180.0]
    residue_group_names = list(residue_group.keys()) # e.g., ['DPPC', 'DUPC', 'CHOL']

    for angle in angles_to_test:
        print(f"Running analysis for max_normal_angle_deg: {angle}")
        cls_instance = Area(u, residue_group, k=k_value, file_path=None, max_normal_angle_deg=angle)
        
        t_start = time.time()
        cls_instance.run() 
        t_end = time.time()
        print(f"Analysis for angle {angle} took {t_end - t_start:.2f} seconds.")

        if cls_instance.Area is not None and cls_instance.Area.size > 0:
            # 为当前角度创建一个新的figure，包含1行3列的子图
            fig, axes = plt.subplots(1, len(residue_group_names), figsize=(18, 6), sharey=True)
            # 如果只有一个残基类型，axes不会是数组，而是单个Axes对象，进行调整
            if len(residue_group_names) == 1:
                axes = [axes]
            
            fig.suptitle(f'Area Distribution by Residue Type (Angle Cutoff: {angle if angle is not None else "None"}°, k={k_value})', fontsize=16)

            all_areas_for_angle = cls_instance.Area # Shape (n_residues, n_frames)
            all_resnames_for_angle = cls_instance.resnames # Shape (n_residues,)

            # 首先计算所有数据的范围，用于统一横坐标
            all_valid_areas = []
            for res_type in residue_group_names:
                type_mask = (all_resnames_for_angle == res_type)
                if np.any(type_mask):
                    areas_for_type = all_areas_for_angle[type_mask, :].flatten()
                    areas_cleaned = areas_for_type[~np.isnan(areas_for_type)]
                    if areas_cleaned.size > 0:
                        all_valid_areas.extend(areas_cleaned)
            
            if all_valid_areas:
                global_min = min(all_valid_areas)
                global_max = max(all_valid_areas)
                # 添加一些边距
                margin = (global_max - global_min) * 0.05
                x_min = global_min - margin
                x_max = global_max + margin
            else:
                x_min, x_max = 0, 1  # 默认范围

            for i, res_type in enumerate(residue_group_names):
                ax = axes[i]
                # 筛选当前残基类型的索引
                type_mask = (all_resnames_for_angle == res_type)
                
                if np.any(type_mask): # 确保该类型存在于分析中
                    areas_for_type = all_areas_for_angle[type_mask, :].flatten()
                    areas_cleaned = areas_for_type[~np.isnan(areas_for_type)]
                    
                    if areas_cleaned.size > 0:
                        # 计算统计信息
                        mean_area = np.mean(areas_cleaned)
                        std_area = np.std(areas_cleaned)
                        
                        # 绘制直方图
                        ax.hist(areas_cleaned, bins=30, alpha=0.75, density=True, label=res_type)
                        
                        # 添加统计信息文本
                        stats_text = f'Mean: {mean_area:.3f}\nStd: {std_area:.3f}'
                        ax.text(0.95, 0.95, stats_text,
                               transform=ax.transAxes,
                               verticalalignment='top',
                               horizontalalignment='right',
                               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                        
                        ax.set_title(f'{res_type}')
                        ax.set_xlabel('Area (nm^2)')
                        if i == 0: # 只给最左边的子图设置Y轴标签
                            ax.set_ylabel('Density')
                        ax.grid(True, linestyle='--', alpha=0.7)
                        
                        # 设置统一的x轴范围
                        ax.set_xlim(0, 2.25)
                    else:
                        ax.set_title(f'{res_type} (No valid data)')
                        print(f"No valid area data for {res_type} at angle {angle} after cleaning.")
                else:
                    ax.set_title(f'{res_type} (Not in selection)')
                    print(f"{res_type} not found in simulation for angle {angle}.")

            plt.tight_layout(rect=[0, 0, 1, 0.96]) # 调整布局以适应suptitle
            plt.show()
        else:
            print(f"No area data generated for angle {angle}.")

    # The main plotting loop is now complete.
    # Original code for a final specific run can be un-commented if needed.

    # 如果您还想运行一次特定设置并保存CSV和打印标准差，可以在这里单独进行：
    # print("\nRunning final analysis with angle 180 for CSV and std dev:")
    # csv_file_final = "cases/csv/area_test_angle180.csv"
    # cls_final = Area(u, residue_group, k=k_value, file_path=csv_file_final, max_normal_angle_deg=180.0)
    # cls_final.run()
    # if cls_final.Area is not None:
    #     print(f'Standard deviation for 180° (from cls_final.Area): {np.std(cls_final.Area)}')
    # print(f'Results saved to {csv_file_final}')

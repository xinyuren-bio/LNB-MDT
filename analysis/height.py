import warnings
import os
import sys
import argparse
import ast
warnings.filterwarnings('ignore')

import MDAnalysis as mda
import numpy as np
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

__all__ = ['Height']

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


class Height(AnalysisBase):
    """
    A class for calculating the height of lipid molecules in LNB system.
    """

    def __init__(self, universe, residuesGroup: dict, k: int = 20, filePath: str = None, 
                 parallel: bool = False, n_jobs: int = -1):
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residuesGroup)
        self.k = k
        self.file_path = filePath
        self.parallel = parallel
        self.n_jobs = n_jobs
        
        # Height分析支持的图表类型
        self.supported_figure_types = ['Line Chart', 'Bar Chart']
        
        # 存储绘图数据，用于后续绘图
        self.plot_data = None

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

        self.parameters = str(residuesGroup) + 'K:' + str(self.k)

    def _prepare(self):
        self.results.Height = np.full([self._n_residues, self.n_frames],
                                      fill_value=np.nan)

    @staticmethod
    def _calculate_height_for_frame(head_pos, tail_pos, k, res_arrange):
        if head_pos.shape[0] == 0:
            return np.array([])
        normals = get_normals(k, head_pos)
        head_to_tail = head_pos - tail_pos
        distance = (np.abs(np.einsum('ij,ij->i', normals, head_to_tail)))[res_arrange]
        return distance * 0.1

    def _single_frame(self):
        head_pos = self.headAtoms.positions
        tail_pos = self.tailAtoms.center_of_geometry(compound='residues')[self.resArrange]
        heights = Height._calculate_height_for_frame(head_pos, tail_pos, self.k, self.resArrange)
        self.results.Height[:, self._frame_index] = heights

    @property
    def Height(self):
        return self.results.Height

    def _get_inputs_generator(self):
        for _ in self.u.trajectory[self.start:self.stop:self.step]:
            head_pos = self.headAtoms.positions.copy()
            tail_pos = self.tailAtoms.center_of_geometry(compound='residues')
            yield (head_pos, tail_pos, self.k, self.resArrange)

    def run(self, start=None, stop=None, step=None, verbose=None, callBack=None):
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
                delayed(Height._calculate_height_for_frame)(*inputs) for inputs in inputs_generator
            )
            if results_list:
                self.results.Height = np.array(results_list).T
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose, callBack=callBack)

        self._conclude()

    def _conclude(self):
        if self.file_path:
            lipids_ratio = {sp: self.u.select_atoms(f'resname {sp}').n_residues for sp in self.residues}
            dict_parameter = {
                'frames': list(range(self.start, self.stop, self.step)),
                'resids': self.resids,
                'resnames': self.resnames,
                'positions': self.headAtoms.positions,
                'results': self.results.Height,
                'file_path': self.file_path,
                'description': 'Height (nm)',
                'parameters': self.parameters,
                'lipids_type': lipids_ratio,
                'trajectory': self._trajectory
            }
            # Assuming WriteExcelLipids is defined and available
            WriteExcelLipids(**dict_parameter).run()
            print(f"Analysis complete. Results saved to {self.file_path}")
            
            # 准备绘图数据
            self._prepare_plot_data()
        else:
            return self.results.Height
    
    def _prepare_plot_data(self):
        """准备绘图数据，格式化为DataFrame"""
        import pandas as pd
        
        # 创建绘图数据
        frames = list(range(self.start, self.stop, self.step))
        data = []
        
        for i, resname in enumerate(self.resnames):
            row = {
                'Resid': self.resids[i],
                'Resname': resname,
                'Coordinates': f"{self.headAtoms.positions[i][0]:.3f},{self.headAtoms.positions[i][1]:.3f},{self.headAtoms.positions[i][2]:.3f}"
            }
            # 添加每个时间帧的高度数据
            for j, frame in enumerate(frames):
                row[str(frame)] = self.results.Height[i, j]
            data.append(row)
        
        self.plot_data = pd.DataFrame(data)
        print(f"Plot data prepared: {self.plot_data.shape}")
    
    def plot_line(self, figure_settings=None):
        """绘制Height数据的折线图"""
        if self.plot_data is None:
            self._prepare_plot_data()
        
        if figure_settings is None:
            figure_settings = {
                'x_title': 'Time (ns)',
                'y_title': 'Height (nm)',
                'axis_text': 12,
                'marker_shape': 'o',
                'marker_size': 0,
                'line_color': 'darkblue'
            }
        
        # 直接使用matplotlib绘制
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        
        # 获取时间列（第4列及以后）
        time_columns = self.plot_data.columns[3:]
        x_axis = time_columns.astype(float) * self._trajectory.dt / 1000  # 转换为ns
        
        # 按残基类型分组绘制
        residue_groups = self.plot_data.groupby('Resname')
        colors = ['darkblue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
        
        for i, (resname, group) in enumerate(residue_groups):
            color = colors[i % len(colors)]
            # 计算该残基类型的平均值
            mean_values = group.iloc[:, 3:].mean(axis=0).values
            plt.plot(x_axis, mean_values,
                    marker=figure_settings.get('marker_shape', 'o'),
                    markersize=figure_settings.get('marker_size', 0),
                    color=color,
                    label=f'{resname}')
        
        plt.xlabel(figure_settings.get('x_title', 'Time (ns)'), fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Height (nm)'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Height Analysis Results', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    
    def plot_bar(self, figure_settings=None):
        """绘制Height数据的条形图"""
        if self.plot_data is None:
            self._prepare_plot_data()
        
        if figure_settings is None:
            figure_settings = {
                'x_title': 'Time (ns)',
                'y_title': 'Height (nm)',
                'axis_text': 12,
                'bar_color': 'steelblue'
            }
        
        # 直接使用matplotlib绘制
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        
        # 获取时间列（第4列及以后）
        time_columns = self.plot_data.columns[3:]
        x_axis = time_columns.astype(float) * self._trajectory.dt / 1000  # 转换为ns
        
        # 按残基类型分组计算平均值
        residue_groups = self.plot_data.groupby('Resname')
        residue_names = []
        mean_heights = []
        colors = ['steelblue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
        
        for i, (resname, group) in enumerate(residue_groups):
            residue_names.append(resname)
            # 计算该残基类型在所有时间点的平均高度
            mean_height = group.iloc[:, 3:].mean().mean()
            mean_heights.append(mean_height)
        
        bars = plt.bar(residue_names, mean_heights,
                      color=colors[:len(residue_names)],
                      alpha=0.7)
        
        plt.xlabel('Residue Type', fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Height (nm)'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Height Analysis Results (Average)', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    
    def plot_scatter(self, figure_settings=None):
        """绘制Height数据的散点图"""
        if self.plot_data is None:
            self._prepare_plot_data()
        
        if figure_settings is None:
            figure_settings = {
                'x_title': 'Time (ns)',
                'y_title': 'Height (nm)',
                'axis_text': 12,
                'scatter_color': 'navy',
                'scatter_size': 50
            }
        
        # 直接使用matplotlib绘制
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        
        # 获取时间列（第4列及以后）
        time_columns = self.plot_data.columns[3:]
        x_axis = time_columns.astype(float) * self._trajectory.dt / 1000  # 转换为ns
        
        # 按残基类型分组绘制散点
        residue_groups = self.plot_data.groupby('Resname')
        colors = ['navy', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
        
        for i, (resname, group) in enumerate(residue_groups):
            color = colors[i % len(colors)]
            # 计算该残基类型的平均值
            mean_values = group.iloc[:, 3:].mean(axis=0).values
            plt.scatter(x_axis, mean_values,
                       color=color,
                       s=figure_settings.get('scatter_size', 50),
                       alpha=0.7,
                       label=f'{resname}')
        
        plt.xlabel(figure_settings.get('x_title', 'Time (ns)'), fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Height (nm)'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Height Analysis Results', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()


# --- Command-line Argument Parsing ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform Height analysis on molecular dynamics trajectories."
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
        default="cases/csv/height_results.csv",
        help="Path to the output CSV file for height results."
    )
    parser.add_argument(
        "--residues", "-r",
        type=str,
        default="{'DPPC': (['PO4'], ['C4B', 'C4A']), 'DUPC':(['PO4'], ['C3A', 'C4B']), 'CHOL':(['ROH'], ['R5'])}",
        help="A dictionary string defining residue groups for analysis. E.g., \"{'DPPC': (['PO4'], ['C4B', 'C4A']), 'CHOL':(['ROH'], ['R5'])}\""
    )
    parser.add_argument(
        "--k-value", "-k",
        type=int,
        default=20,
        help="K value for height calculation."
    )
    parser.add_argument(
        "--parallel", "-p",
        action="store_true",
        help="Enable parallel processing for height calculation."
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
        print("Please ensure it's a valid dictionary string, e.g., \"{'DPPC': (['PO4'], ['C4B', 'C4A']), 'CHOL':(['ROH'], ['R5'])}\"")
        sys.exit(1)

    print("\n--- Initializing MDAnalysis Universe ---")
    try:
        u = mda.Universe(args.gro_file, args.xtc_file)
    except Exception as e:
        print(f"Error loading MDAnalysis Universe: {e}")
        print("Please check if GRO/XTC files exist and are valid.")
        sys.exit(1)

    print("\n--- Running Height Analysis ---")
    height_analysis = Height(
        u,
        residues_group_parsed,
        k=args.k_value,
        filePath=args.output_csv,
        parallel=args.parallel,
        n_jobs=args.n_jobs
    )
    height_analysis.run(
        start=args.start_frame,
        stop=args.stop_frame,
        step=args.step_frame,
        verbose=args.verbose
    )

    print("\n--- Analysis Finished ---")
    print(f"Results saved to: {args.output_csv}")
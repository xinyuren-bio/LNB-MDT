import warnings
import os
import sys
import argparse
import ast
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

    def __init__(self, universe, residues_group: dict, cutoff: float = 8.0, filePath: str = None,
                 parallel: bool = False, n_jobs: int = -1):
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residues_group)
        self.cutoff = cutoff
        self.file_path = filePath
        self.parallel = parallel
        self.n_jobs = n_jobs
        
        # Cluster分析支持的图表类型
        self.supported_figure_types = ['Line Chart', 'Bar Chart']
        
        # 存储绘图数据，用于后续绘图
        self.plot_data = None

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

        # 过滤掉超出范围的索引
        valid_mask = (residue_indices_pairs[:, 0] < n_residues) & (residue_indices_pairs[:, 1] < n_residues)
        residue_indices_pairs = residue_indices_pairs[valid_mask]
        
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

    def run(self, start=None, stop=None, step=None, verbose=None, callBack=None):
        # 处理stop参数：-1应该转换为None，让MDAnalysis的check_slice_indices正确处理
        if stop == -1:
            stop = None  # 转换为None，让MDAnalysis处理为完整的轨迹长度
        
        # 对于并行模式，需要先计算实际的stop值用于n_frames计算
        if stop is None:
            actual_stop = self._trajectory.n_frames
        elif stop < 0:
            actual_stop = self._trajectory.n_frames + stop + 1
        else:
            actual_stop = stop if stop <= self._trajectory.n_frames else self._trajectory.n_frames
        
        self.start = start if start is not None else 0
        self.stop = actual_stop
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
            super().run(start=start, stop=stop, step=step, verbose=verbose, callBack=callBack)

        self._conclude()

    def _conclude(self):
        if self.file_path:
            dict_parameter = {
                'frames': list(range(self.start, self.stop, self.step)),
                'results': self.results.LargestClusterSize,
                'file_path': self.file_path,
                'description': 'Largest Cluster Size',
                'parameters': self.parameters,
                'trajectory': self._trajectory
            }
            WriteExcelBubble(**dict_parameter).run()
            print(f"Analysis complete. Results will be saved to {self.file_path}")
            
            # 准备绘图数据
            self._prepare_plot_data()
        else:
            return self.results.LargestClusterSize
    
    def _prepare_plot_data(self):
        """准备绘图数据，格式化为DataFrame（BubbleFigure格式）"""
        import pandas as pd
        
        # 创建绘图数据（BubbleFigure格式）
        frames = list(range(self.start, self.stop, self.step))
        data = []
        
        # Cluster是整体数据，不是按残基分组
        for j, frame in enumerate(frames):
            row = {
                'Time(ns)': frame * self._trajectory.dt / 1000,  # 转换为ns
                'Value': self.results.LargestClusterSize[j]
            }
            data.append(row)
        
        self.plot_data = pd.DataFrame(data)
        print(f"Cluster plot data prepared: {self.plot_data.shape}")
    
    def plot_line(self, figure_settings=None):
        """绘制Cluster数据的折线图"""
        if self.plot_data is None:
            self._prepare_plot_data()
        
        if figure_settings is None:
            figure_settings = {
                'x_title': 'Time (ns)',
                'y_title': 'Largest Cluster Size',
                'axis_text': 12,
                'marker_shape': 'o',
                'marker_size': 0,
                'line_color': 'purple'
            }
        
        # 直接使用matplotlib绘制
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        plt.plot(self.plot_data['Time(ns)'], self.plot_data['Value'],
                marker=figure_settings.get('marker_shape', 'o'),
                markersize=figure_settings.get('marker_size', 0),
                color=figure_settings.get('line_color', 'purple'),
                label='Cluster Size')
        
        plt.xlabel(figure_settings.get('x_title', 'Time (ns)'), fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Largest Cluster Size'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Cluster Analysis Results', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    
    def plot_bar(self, figure_settings=None):
        """绘制Cluster数据的条形图"""
        if self.plot_data is None:
            self._prepare_plot_data()
        
        if figure_settings is None:
            figure_settings = {
                'x_title': 'Time (ns)',
                'y_title': 'Largest Cluster Size',
                'axis_text': 12,
                'bar_color': 'plum'
            }
        
        # 直接使用matplotlib绘制
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        plt.bar(self.plot_data['Time(ns)'], self.plot_data['Value'],
               color=figure_settings.get('bar_color', 'plum'),
               alpha=0.7,
               label='Cluster Size')
        
        plt.xlabel(figure_settings.get('x_title', 'Time (ns)'), fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Largest Cluster Size'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Cluster Analysis Results', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()


# --- Command-line Argument Parsing ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform Cluster analysis on molecular dynamics trajectories."
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
        default="cases/csv/cluster_results.csv",
        help="Path to the output CSV file for cluster results."
    )
    parser.add_argument(
        "--residues", "-r",
        type=str,
        default="{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}",
        help="A dictionary string defining residue groups for analysis. E.g., \"{'DPPC': ['PO4'], 'CHOL': ['ROH']}\""
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=8.0,
        help="Cutoff distance for clustering (in Angstroms)."
    )
    parser.add_argument(
        "--parallel", "-p",
        action="store_true",
        help="Enable parallel processing for cluster calculation."
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

    print("\n--- Running Cluster Analysis ---")
    cluster_analysis = Cluster(
        u,
        residues_group_parsed,
        cutoff=args.cutoff,
        filePath=args.output_csv,
        parallel=args.parallel,
        n_jobs=args.n_jobs
    )
    cluster_analysis.run(
        start=args.start_frame,
        stop=args.stop_frame,
        step=args.step_frame,
        verbose=args.verbose
    )

    print("\n--- Analysis Finished ---")
    print(f"Results saved to: {args.output_csv}")
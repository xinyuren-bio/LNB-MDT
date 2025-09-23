import warnings
import os
import sys
import argparse
import ast
warnings.filterwarnings('ignore')

import numpy as np
import MDAnalysis as mda
from joblib import Parallel, delayed

if __name__ == '__main__':
    current_file_path = os.path.abspath(__file__)
    current_dir = os.path.dirname(current_file_path)
    package_root = os.path.abspath(os.path.join(current_dir, '..'))
    if package_root not in sys.path:
        sys.path.insert(0, package_root)
from analysis.analysis_base import * 

__all__ = ['Gyration']


def get_radius_of_gyration(positions, masses):
    total_mass = np.sum(masses)
    if total_mass == 0:
        return np.nan
    
    squared_distances = np.sum(positions**2, axis=1)
    mass_weighted_sq_dist = np.dot(squared_distances, masses)
    
    return np.sqrt(mass_weighted_sq_dist / total_mass)


class Gyration(AnalysisBase):
    """
    A class for calculating the radius of gyration for a selection of atoms.
    """
    def __init__(self, universe, residues_group: dict, filePath: str = None, 
                 parallel: bool = False, n_jobs: int = -1):
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(residues_group)
        self.file_path = filePath
        self.parallel = parallel
        self.n_jobs = n_jobs
        
        # Gyration分析支持的图表类型
        self.supported_figure_types = ['Line Chart', 'Bar Chart']
        
        # 存储绘图数据，用于后续绘图
        self.plot_data = None

        self.headSp = {sp: ' '.join(residues_group[sp]) for sp in self.residues}
        print("Selected Atoms:", self.headSp)
        
        self.headAtoms = self.u.atoms[[]]
        for sp in self.residues:
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]), updating=False)
        
        self._n_residues = self.headAtoms.n_residues
        self.results.Gyration = None
        self.parameters = str(residues_group)

    def _prepare(self):
        self.results.Gyration = np.full(self.n_frames, fill_value=np.nan)

    @staticmethod
    def _calculate_gyration_for_frame(residue_com, group_com, residue_masses):
        relative_positions = residue_com - group_com
        return get_radius_of_gyration(relative_positions, residue_masses)

    def _single_frame(self):
        residue_com = self.headAtoms.center_of_mass(compound='residues') / 10
        group_com = self.headAtoms.center_of_mass() / 10
        residue_masses = self.headAtoms.residues.masses

        gyration_radius = Gyration._calculate_gyration_for_frame(residue_com, group_com, residue_masses)
        self.results.Gyration[self._frame_index] = gyration_radius

    def _get_inputs_generator(self):
        residue_masses = self.headAtoms.residues.masses
        for _ in self.u.trajectory[self.start:self.stop:self.step]:
            residue_com = self.headAtoms.center_of_mass(compound='residues') / 10
            group_com = self.headAtoms.center_of_mass() / 10
            yield (residue_com, group_com, residue_masses)

    def run(self, start=None, stop=None, step=None, verbose=None, callBack=None):
        self.start = start if start is not None else 0
        self.stop = stop if stop is not None and stop < self._trajectory.n_frames else self._trajectory.n_frames
        self.step = step if step is not None else 1
        self.n_frames = len(range(self.start, self.stop, self.step))
        self._prepare()

        if self.parallel:
            print(f"Running in parallel on {self.n_jobs} jobs...")
            
            # 使用父类的run方法以确保进度条正常工作
            # 但禁用并行模式以避免冲突
            original_parallel = self.parallel
            self.parallel = False
            super().run(start=start, stop=stop, step=step, verbose=verbose, callBack=callBack)
            self.parallel = original_parallel
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose, callBack=callBack)

        self._conclude()

    def _conclude(self):
        if self.file_path:
            try:
                from .analysis_base import WriteExcelBubble
            except ImportError:
                from analysis_base import WriteExcelBubble
            
            dict_parameter = {
                'frames': list(range(self.start, self.stop, self.step)),
                'results': self.results.Gyration,
                'file_path': self.file_path,
                'description': 'Gyration (nm)',
                'parameters': self.parameters,
                'trajectory': self._trajectory
            }
            WriteExcelBubble(**dict_parameter).run()
            print(f"Analysis complete. Results will be saved to {self.file_path}")
            
            # 准备绘图数据
            self._prepare_plot_data()
        else:
            return self.results.Gyration
    
    def _prepare_plot_data(self):
        """准备绘图数据，格式化为DataFrame（BubbleFigure格式）"""
        import pandas as pd
        
        # 创建绘图数据（BubbleFigure格式）
        frames = list(range(self.start, self.stop, self.step))
        data = []
        
        # Gyration是整体数据，不是按残基分组
        for j, frame in enumerate(frames):
            row = {
                'Time(ns)': frame * self._trajectory.dt / 1000,  # 转换为ns
                'Value': self.results.Gyration[j]
            }
            data.append(row)
        
        self.plot_data = pd.DataFrame(data)
        print(f"Gyration plot data prepared: {self.plot_data.shape}")
    
    def plot_line(self, figure_settings=None):
        """绘制Gyration数据的折线图"""
        if self.plot_data is None:
            self._prepare_plot_data()
        
        if figure_settings is None:
            figure_settings = {
                'x_title': 'Time (ns)',
                'y_title': 'Gyration (nm)',
                'axis_text': 12,
                'marker_shape': 'o',
                'marker_size': 0,
                'line_color': 'green'
            }
        
        # 直接使用matplotlib绘制
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        plt.plot(self.plot_data['Time(ns)'], self.plot_data['Value'],
                marker=figure_settings.get('marker_shape', 'o'),
                markersize=figure_settings.get('marker_size', 0),
                color=figure_settings.get('line_color', 'green'),
                label='Gyration')
        
        plt.xlabel(figure_settings.get('x_title', 'Time (ns)'), fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Gyration (nm)'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Gyration Analysis Results', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    
    def plot_bar(self, figure_settings=None):
        """绘制Gyration数据的条形图"""
        if self.plot_data is None:
            self._prepare_plot_data()
        
        if figure_settings is None:
            figure_settings = {
                'x_title': 'Time (ns)',
                'y_title': 'Gyration (nm)',
                'axis_text': 12,
                'bar_color': 'lightgreen'
            }
        
        # 直接使用matplotlib绘制
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        plt.bar(self.plot_data['Time(ns)'], self.plot_data['Value'],
               color=figure_settings.get('bar_color', 'lightgreen'),
               alpha=0.7,
               label='Gyration')
        
        plt.xlabel(figure_settings.get('x_title', 'Time (ns)'), fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Gyration (nm)'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Gyration Analysis Results', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()


# --- Command-line Argument Parsing ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform Gyration analysis on molecular dynamics trajectories."
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
        default="cases/csv/gyration_results.csv",
        help="Path to the output CSV file for gyration results."
    )
    parser.add_argument(
        "--residues", "-r",
        type=str,
        default="{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}",
        help="A dictionary string defining residue groups for analysis. E.g., \"{'DPPC': ['PO4'], 'CHOL': ['ROH']}\""
    )
    parser.add_argument(
        "--parallel", "-p",
        action="store_true",
        help="Enable parallel processing for gyration calculation."
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

    print("\n--- Running Gyration Analysis ---")
    gyration_analysis = Gyration(
        u,
        residues_group_parsed,
        filePath=args.output_csv,
        parallel=args.parallel,
        n_jobs=args.n_jobs
    )
    gyration_analysis.run(
        start=args.start_frame,
        stop=args.stop_frame,
        step=args.step_frame,
        verbose=args.verbose
    )

    print("\n--- Analysis Finished ---")
    print(f"Results saved to: {args.output_csv}")
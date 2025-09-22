import warnings
import os
import sys
import argparse 
import ast      
warnings.filterwarnings('ignore')

import pandas as pd
from scipy.linalg import eigh
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

__all__ = ['Anisotropy']


class Anisotropy(AnalysisBase):
    def __init__(self, universe: mda.Universe, residueGroup: dict, filePath: str = None, parallel: bool = False, n_jobs: int = -1):
        super().__init__(universe.trajectory)
        self.u = universe
        self.file_path = filePath
        self.parallel = parallel
        self.n_jobs = n_jobs
        
        # Anisotropy分析支持的图表类型
        self.supported_figure_types = ['Line Chart', 'Bar Chart']
        
        # 存储绘图数据，用于后续绘图
        self.plot_data = None
        
        self.residues = list(residueGroup)
        self.headSp = {sp: ' '.join(residueGroup[sp]) for sp in residueGroup}
        print("Head atoms:", self.headSp)

        self.headAtoms = self.u.atoms[[]]

        # 为了保持原子顺序，循环选择指定脂质类型的头部原子并添加到 AtomGroup
        for i in range(len(self.residues)):
            self.headAtoms += self.u.select_atoms(
                f'resname {self.residues[i]} and name {self.headSp[self.residues[i]]}', 
                updating=False
            )
        
        if self.headAtoms.n_atoms == 0:
            raise ValueError("Atom selection is empty. Please check your `residueGroup` dictionary and atomic names.")

        self._n_residues = self.headAtoms.n_residues
        self.resids = self.headAtoms.resids
        self.resnames = self.headAtoms.resnames
        self.results.Anisotropy = None

        self.parameters = str(residueGroup)

    @property
    def Anisotropy(self):
        return self.results.Anisotropy

    def _prepare(self):
        self.results.Anisotropy = np.full(self.n_frames, fill_value=np.nan)

    @staticmethod
    def _calculate_anisotropy(positions: np.ndarray) -> float: 
        covariance_matrix = np.cov(positions.T)
        eigenvalues, _ = eigh(covariance_matrix)
        anisotropy = np.sqrt(1 - eigenvalues[0]**2 / eigenvalues[2]**2)
        return anisotropy

    def _single_frame(self):
        positions = self.headAtoms.positions
        self.results.Anisotropy[self._frame_index] = self._calculate_anisotropy(positions)

    def _get_positions_generator(self):
        for ts in self.u.trajectory[self.start:self.stop:self.step]:
            yield self.headAtoms.positions

    def run(self, start=None, stop=None, step=None, verbose=None, callBack=None):
        if self.parallel:
            self.start = start if start is not None else 0
            self.stop = stop if stop is not None and stop < self._trajectory.n_frames else self._trajectory.n_frames
            self.step = step if step is not None else 1
            self.n_frames = len(range(self.start, self.stop, self.step))
            self._prepare()
            print(f"Running in parallel on {self.n_jobs} jobs...")
            positions_generator = self._get_positions_generator()
            verbose_level = 10 if verbose else 0
            results_list = Parallel(n_jobs=self.n_jobs, verbose=verbose_level)(
                delayed(self._calculate_anisotropy)(pos) for pos in positions_generator
            )
            self.results.Anisotropy = np.array(results_list)
            self._conclude()
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose, callBack=callBack)
            
    def _conclude(self):
        if self.file_path:
            dict_parameter = {
                'frames': list(range(self.start, self.stop, self.step)),
                'results': self.results.Anisotropy,
                'file_path': self.file_path,
                'description': 'Anisotropy',
                'parameters': self.parameters,
                'trajectory': self._trajectory
            }
            WriteExcelBubble(**dict_parameter).run()
            print(f"Results saved to {self.file_path}")
            
            # 准备绘图数据
            self._prepare_plot_data()
        else:
            return self.results.Anisotropy
    
    def _prepare_plot_data(self):
        """准备绘图数据，格式化为DataFrame（BubbleFigure格式）"""
        import pandas as pd
        
        # 创建绘图数据（BubbleFigure格式）
        frames = list(range(self.start, self.stop, self.step))
        data = []
        
        # Anisotropy是整体数据，不是按残基分组
        for j, frame in enumerate(frames):
            row = {
                'Time(ns)': frame * self._trajectory.dt / 1000,  # 转换为ns
                'Value': self.results.Anisotropy[j]
            }
            data.append(row)
        
        self.plot_data = pd.DataFrame(data)
        print(f"Anisotropy plot data prepared: {self.plot_data.shape}")
    
    def plot_line(self, figure_settings=None):
        """绘制Anisotropy数据的折线图"""
        if self.plot_data is None:
            self._prepare_plot_data()
        
        if figure_settings is None:
            figure_settings = {
                'x_title': 'Time (ns)',
                'y_title': 'Anisotropy',
                'axis_text': 12,
                'marker_shape': 'o',
                'marker_size': 0,
                'line_color': 'blue'
            }
        
        # 直接使用matplotlib绘制
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        plt.plot(self.plot_data['Time(ns)'], self.plot_data['Value'],
                marker=figure_settings.get('marker_shape', 'o'),
                markersize=figure_settings.get('marker_size', 0),
                color=figure_settings.get('line_color', 'blue'),
                label='Anisotropy')
        
        plt.xlabel(figure_settings.get('x_title', 'Time (ns)'), fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Anisotropy'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Anisotropy Analysis Results', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    
    def plot_bar(self, figure_settings=None):
        """绘制Anisotropy数据的条形图"""
        if self.plot_data is None:
            self._prepare_plot_data()
        
        if figure_settings is None:
            figure_settings = {
                'x_title': 'Time (ns)',
                'y_title': 'Anisotropy',
                'axis_text': 12,
                'bar_color': 'lightblue'
            }
        
        # 直接使用matplotlib绘制
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        plt.bar(self.plot_data['Time(ns)'], self.plot_data['Value'],
               color=figure_settings.get('bar_color', 'lightblue'),
               alpha=0.7,
               label='Anisotropy')
        
        plt.xlabel(figure_settings.get('x_title', 'Time (ns)'), fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Anisotropy'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Anisotropy Analysis Results', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()


# --- Command-line Argument Parsing ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform Anisotropy analysis on molecular dynamics trajectories."
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
        default="cases/csv/anisotropy_results.csv",
        help="Path to the output CSV file for anisotropy results."
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
        help="Enable parallel processing for anisotropy calculation."
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

    print("\n--- Running Anisotropy Analysis ---")
    anisotropy_analysis = Anisotropy(
        u,
        residues_group_parsed,
        filePath=args.output_csv,
        parallel=args.parallel,
        n_jobs=args.n_jobs
    )
    anisotropy_analysis.run(
        start=args.start_frame,
        stop=args.stop_frame,
        step=args.step_frame,
        verbose=args.verbose
    )

    print("\n--- Analysis Finished ---")
    print(f"Results saved to: {args.output_csv}")


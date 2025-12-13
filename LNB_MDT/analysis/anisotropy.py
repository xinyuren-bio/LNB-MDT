import warnings
import os
import sys
import argparse 
import ast      
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import MDAnalysis as mda
from joblib import Parallel, delayed

if __name__ == '__main__':
    current_file_path = os.path.abspath(__file__)
    current_dir = os.path.dirname(current_file_path)
    package_root = os.path.abspath(os.path.join(current_dir, '..'))
    if package_root not in sys.path:
        sys.path.insert(0, package_root)
from .analysis_base import AnalysisBase, WriteExcelBubble

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
        
        # 处理原子选择：支持单个原子、多个原子的列表，或 'all' 关键字
        self.headSp = {}
        self.use_all_atoms = {}  # 标记哪些残基类型使用 'all'
        
        for sp in residueGroup:
            value = residueGroup[sp]
            self.use_all_atoms[sp] = False
            
            # 检查是否是字符串 'all'
            if isinstance(value, str) and value.lower() == 'all':
                self.use_all_atoms[sp] = True
                self.headSp[sp] = None
                print(f"  {sp}: selecting ALL atoms")
                continue
            
            # 检查是否是列表/元组
            if isinstance(value, (list, tuple)):
                # 如果是嵌套结构
                if len(value) > 0 and isinstance(value[0], (list, tuple)):
                    # 检查嵌套列表的第一个元素是否包含 'all'
                    if 'all' in [item.lower() for item in value[0] if isinstance(item, str)]:
                        self.use_all_atoms[sp] = True
                        self.headSp[sp] = None
                        print(f"  {sp}: selecting ALL atoms")
                        continue
                    else:
                        # 使用第一个嵌套列表（头部原子）
                        self.headSp[sp] = ' '.join([str(item) for item in value[0]])
                # 如果是普通列表
                else:
                    # 检查列表中是否包含 'all'
                    if 'all' in [item.lower() if isinstance(item, str) else str(item).lower() for item in value]:
                        self.use_all_atoms[sp] = True
                        self.headSp[sp] = None
                        print(f"  {sp}: selecting ALL atoms")
                        continue
                    else:
                        # 使用整个列表
                        self.headSp[sp] = ' '.join([str(item) for item in value])
            else:
                # 单个原子（向后兼容）
                self.headSp[sp] = str(value)
        
        # 只显示非 'all' 的原子选择
        non_all_atoms = {k: v for k, v in self.headSp.items() if v is not None}
        if non_all_atoms:
            print("Head atoms:", non_all_atoms)
        all_atoms_species = {k: v for k, v in self.use_all_atoms.items() if v}
        if all_atoms_species:
            print("Using 'all' atoms for:", list(all_atoms_species.keys()))

        self.headAtoms = self.u.atoms[[]]

        # 为了保持原子顺序，循环选择指定脂质类型的头部原子并添加到 AtomGroup
        for i in range(len(self.residues)):
            sp = self.residues[i]
            if self.use_all_atoms.get(sp, False):
                # 选择该残基类型的所有原子
                selected_atoms = self.u.select_atoms(
                    f'resname {sp}', 
                    updating=False
                )
            else:
                # 选择指定的原子
                selected_atoms = self.u.select_atoms(
                    f'resname {sp} and name {self.headSp[sp]}', 
                    updating=False
                )
            self.headAtoms += selected_atoms
        
        if self.headAtoms.n_atoms == 0:
            raise ValueError("Atom selection is empty. Please check your `residueGroup` dictionary and atomic names.")

        # 检查每个残基是否选择了多个原子
        self._n_residues = self.headAtoms.n_residues
        self.resids = self.headAtoms.resids
        self.resnames = self.headAtoms.resnames
        
        # 检查是否有残基选择了多个原子
        residue_atom_counts = {}
        for resid in np.unique(self.headAtoms.resids):
            residue_atoms = self.headAtoms[self.headAtoms.resids == resid]
            residue_atom_counts[resid] = residue_atoms.n_atoms
        
        max_atoms_per_residue = max(residue_atom_counts.values()) if residue_atom_counts else 1
        
        # 标记是否使用几何中心（避免每帧重复检查）
        self._use_center_of_geometry = max_atoms_per_residue > 1
        
        if self._use_center_of_geometry:
            print(f"Note: Some residues have multiple atoms selected (max: {max_atoms_per_residue}).")
            print("Will compute center of geometry for each residue.")
        else:
            print("Note: Each residue has one atom selected.")
        
        self.results.Anisotropy = None

        self.parameters = str(residueGroup)

    @property
    def Anisotropy(self):
        return self.results.Anisotropy

    def _prepare(self):
        self.results.Anisotropy = np.full(self.n_frames, fill_value=np.nan)

    @staticmethod
    def _calculate_asphericity(positions: np.ndarray) -> float:
        """
        计算κ²值 (相对形状各向异性)
        
        算法步骤：
        1. 计算几何中心
        2. 中心化坐标
        3. 构建回转半径张量S
        4. 计算特征值
        5. 计算不变量I1和I2
        6. 计算κ² = 1 - 3*I2/(I1²)
        
        Args:
            positions: Nx3的坐标数组
            
        Returns:
            κ²值
        """
        # 1. 计算几何中心（平均值）
        cog = np.mean(positions, axis=0)
        
        # 2. 中心化坐标：从每个原子的坐标中减去几何中心
        centered_positions = positions - cog
        
        # 3. 构建回转半径张量S (3x3)
        # S[i,j] = mean(centered_positions[:,i] * centered_positions[:,j])
        S = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                S[i, j] = np.mean(centered_positions[:, i] * centered_positions[:, j])
        
        # 4. 计算对称张量S的特征值
        eigenvalues = np.linalg.eigvalsh(S)
        
        # 5. 计算不变量
        # I1 = λ1 + λ2 + λ3
        I1 = np.sum(eigenvalues)
        
        # I2 = λ1*λ2 + λ1*λ3 + λ2*λ3
        I2 = (eigenvalues[0] * eigenvalues[1] + 
              eigenvalues[0] * eigenvalues[2] + 
              eigenvalues[1] * eigenvalues[2])
        
        # 6. 计算κ²
        if I1 != 0:
            kappa2 = 1 - 3 * I2 / (I1 ** 2)
        else:
            kappa2 = np.nan
        
        return kappa2

    def _single_frame(self):
        # 如果每个残基有多个原子，计算几何中心；否则直接使用原子位置
        if self._use_center_of_geometry:
            # 多个原子：计算每个残基的几何中心
            positions = self.headAtoms.center_of_geometry(compound='residues')
        else:
            # 单个原子：直接使用原子位置
            positions = self.headAtoms.positions
        
        self.results.Anisotropy[self._frame_index] = self._calculate_asphericity(positions)

    def _get_positions_generator(self):
        for ts in self.u.trajectory[self.start:self.stop:self.step]:
            if self._use_center_of_geometry:
                # 多个原子：计算每个残基的几何中心
                yield self.headAtoms.center_of_geometry(compound='residues')
            else:
                # 单个原子：直接使用原子位置
                yield self.headAtoms.positions

    def run(self, start=None, stop=None, step=None, verbose=None, callBack=None, plot_type='none', plot_output=None):
        # 处理stop参数：-1应该转换为None
        if stop == -1:
            stop = None
        
        if self.parallel:
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
            print(f"Running in parallel on {self.n_jobs} jobs...")
            
            # 使用父类的run方法以确保进度条正常工作
            # 但禁用并行模式以避免冲突
            original_parallel = self.parallel
            self.parallel = False
            # 临时保存参数，在_conclude中会使用
            self._plot_type = plot_type
            self._plot_output = plot_output
            super().run(start=start, stop=stop, step=step, verbose=verbose, callBack=callBack)
            self.parallel = original_parallel
        else:
            print("Running in serial mode...")
            # 临时保存参数，在_conclude中会使用
            self._plot_type = plot_type
            self._plot_output = plot_output
            # 串行模式下，super().run()会在最后自动调用_conclude()，会使用保存的参数
            super().run(start=start, stop=stop, step=step, verbose=verbose, callBack=callBack)
            
    def _conclude(self, plot_type=None, plot_output=None):
        # 如果plot_type为None，尝试从实例属性获取（由run方法设置）
        if plot_type is None:
            plot_type = getattr(self, '_plot_type', 'none')
        if plot_output is None:
            plot_output = getattr(self, '_plot_output', None)
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
        
        # 准备绘图数据（无论是否保存文件）
        self._prepare_plot_data()
        
        # 自动绘图（无论是否保存文件）
        if plot_type != 'none':
            self.plot_anisotropy(plot_type, plot_output)
        
        if not self.file_path:
            return self.results.Anisotropy
    
    def plot_anisotropy(self, plot_type='both', plot_output=None):
        """
        绘制非球形度分析结果
        
        Args:
            plot_type: 绘图类型，'line', 'bar', 'both', 或 'none'
            plot_output: 输出文件路径（可选）
        """
        import matplotlib.pyplot as plt
        
        if self.plot_data is None:
            self._prepare_plot_data()
        
        if plot_type == 'none':
            return
        
        # 用于跟踪是否需要在最后显示图形
        should_show = plot_output is None
        
        if plot_type in ['line', 'both']:
            self.plot_line()
            
            if plot_output:
                # 如果是 'both'，添加后缀；如果是 'line'，直接使用或添加扩展名
                if plot_type == 'both':
                    output_file = plot_output.replace('.png', '_line.png') if plot_output.endswith('.png') else f"{plot_output}_line.png"
                else:
                    output_file = plot_output if plot_output.endswith('.png') else f"{plot_output}.png"
                
                # 获取当前图形并保存
                current_fig = plt.gcf()
                current_fig.savefig(output_file, dpi=300, bbox_inches='tight')
                print(f"Line chart saved to: {output_file}")
                plt.close()
        
        if plot_type in ['bar', 'both']:
            self.plot_bar()
            
            if plot_output:
                # 如果是 'both'，添加后缀；如果是 'bar'，直接使用或添加扩展名
                if plot_type == 'both':
                    output_file = plot_output.replace('.png', '_bar.png') if plot_output.endswith('.png') else f"{plot_output}_bar.png"
                else:
                    output_file = plot_output if plot_output.endswith('.png') else f"{plot_output}.png"
                
                # 获取当前图形并保存
                current_fig = plt.gcf()
                current_fig.savefig(output_file, dpi=300, bbox_inches='tight')
                print(f"Bar chart saved to: {output_file}")
                plt.close()
        
        # 如果没有保存文件，显示图形（只显示一次）
        if should_show:
            plt.show(block=False)  # 非阻塞显示，不阻塞GUI线程
    
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
                label='Anisotropy', linewidth=2)
        
        plt.xlabel(figure_settings.get('x_title', 'Time (ns)'), fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Anisotropy'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Anisotropy Analysis Results', fontsize=14, fontweight='bold')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show(block=False)  # 非阻塞显示，不阻塞GUI线程
    
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
               label='Anisotropy',
               edgecolor='black',
               linewidth=0.5)
        
        plt.xlabel(figure_settings.get('x_title', 'Time (ns)'), fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Anisotropy'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Anisotropy Analysis Results', fontsize=14, fontweight='bold')
        plt.legend()
        plt.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        plt.show(block=False)  # 非阻塞显示，不阻塞GUI线程


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
        default=None,
        help="Path to the XTC file (trajectory file). Optional. If not provided, only GRO file will be analyzed (single frame)."
    )
    parser.add_argument(
        "--output-csv", "-o",
        type=str,
        default="cases/csv/asphericity_results.csv",
        help="Path to the output CSV file for asphericity results."
    )
    parser.add_argument(
        "--residues", "-r",
        type=str,
        default="{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}",
        help="A dictionary string defining residue groups for analysis. Supports multiple atoms per residue by providing a list. E.g., \"{'DPPC': ['PO4', 'C1', 'C2'], 'CHOL': ['ROH']}\" or \"{'DPPC': ['PO4'], 'CHOL': ['ROH']}\". When multiple atoms are provided, their center of geometry will be calculated."
    )
    parser.add_argument(
        "--parallel", "-p",
        action="store_true",
        help="Enable parallel processing for asphericity calculation."
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
        default=True,  # Default to True for CLI mode
        help="Enable verbose output during analysis (default: True in CLI mode)."
    )
    parser.add_argument(
        "--no-verbose",
        action="store_false",
        dest="verbose",
        help="Disable verbose output (disable progress bar)."
    )
    parser.add_argument(
        "--plot",
        type=str,
        choices=['line', 'bar', 'both', 'none'],
        default='none',
        help="Type of plot to generate. 'line' for line chart, 'bar' for bar chart, 'both' for both, 'none' to skip plotting."
    )
    parser.add_argument(
        "--plot-output",
        type=str,
        default=None,
        help="Output file path for the plot. If not specified and plot is enabled, plot will only be displayed (not saved). For 'both', use filename without extension (e.g., 'anisotropy_plot' will create 'anisotropy_plot_line.png' and 'anisotropy_plot_bar.png')."
    )

    parser.add_argument(
        "--test", "-test",
        action="store_true",
        help="Run test mode: automatically use test files from LNB_MDT/cases_lnb/ directory."
    )

    return parser.parse_args()


def main():
    """Main function for command-line execution"""
    args = parse_args()

    # Handle test mode
    if args.test:
        # Get the path to LNB_MDT package
        import LNB_MDT
        package_dir = os.path.dirname(os.path.abspath(LNB_MDT.__file__))
        test_dir = os.path.join(package_dir, 'cases_lnb')
        
        # Set test file paths
        args.gro_file = os.path.join(test_dir, 'lnb.gro')
        args.xtc_file = os.path.join(test_dir, 'lnb.xtc')
        args.output_csv = os.path.join(test_dir, 'anisotropy.csv')
        
        # Set default residues for test
        args.residues = "{'DPPC': ['PO4'], 'DAPC': ['PO4'], 'CHOL': ['ROH']}"
        
        # Verify test files exist
        if not os.path.exists(args.gro_file):
            print(f"Error: Test GRO file not found: {args.gro_file}")
            sys.exit(1)
        if not os.path.exists(args.xtc_file):
            print(f"Error: Test XTC file not found: {args.xtc_file}")
            sys.exit(1)
        
        print("\n=== Running in TEST mode ===")
        print(f"GRO file: {args.gro_file}")
        print(f"XTC file: {args.xtc_file}")
        print(f"Output CSV: {args.output_csv}")
        print(f"Residues: {args.residues}")
        print("=" * 40)

    # Parse residues_group from string

    # Parse residues_group from string
    try:
        residues_group_parsed = ast.literal_eval(args.residues)
        if not isinstance(residues_group_parsed, dict):
            raise ValueError("Residues argument must be a dictionary string.")
    except (ValueError, SyntaxError) as e:
        print(f"Error: Could not parse residues argument: {e}")
        print("Please ensure it's a valid dictionary string, e.g., \"{'DPPC': ['PO4'], 'CHOL': ['ROH']}\"")
        sys.exit(1)

    # Check if files exist before trying to load
    if not os.path.exists(args.gro_file):
        print(f"\nError: GRO file not found: {args.gro_file}")
        print(f"Current working directory: {os.getcwd()}")
        print(f"Absolute path attempted: {os.path.abspath(args.gro_file)}")
        print("\nPlease:")
        print("  1. Check if the file path is correct")
        print("  2. Use --gro-file or -g to specify the correct path")
        print("  3. Or use --test to run with test files")
        sys.exit(1)
    
    if args.xtc_file and not os.path.exists(args.xtc_file):
        print(f"\nError: XTC file not found: {args.xtc_file}")
        print(f"Current working directory: {os.getcwd()}")
        print(f"Absolute path attempted: {os.path.abspath(args.xtc_file)}")
        print("\nPlease:")
        print("  1. Check if the file path is correct")
        print("  2. Use --xtc-file or -x to specify the correct path")
        print("  3. Or use --test to run with test files")
        sys.exit(1)
    
    print("\n--- Initializing MDAnalysis Universe ---")
    try:
        # 如果提供了xtc文件，则同时加载gro和xtc文件
        if args.xtc_file and os.path.exists(args.xtc_file):
            u = mda.Universe(args.gro_file, args.xtc_file)
            print(f"Loaded both GRO and XTC files: {args.gro_file}, {args.xtc_file}")
        else:
            # 只使用gro文件
            u = mda.Universe(args.gro_file)
            print(f"Loaded only GRO file: {args.gro_file}")
            print(f"Note: Analyzing single frame (frame 0) from GRO file")
    except Exception as e:
        print(f"\nError loading MDAnalysis Universe: {e}")
        print(f"GRO file: {args.gro_file}")
        if args.xtc_file:
            print(f"XTC file: {args.xtc_file}")
        print("\nPlease check if:")
        print("  1. Files exist and are readable")
        print("  2. Files are valid GROMACS format")
        print("  3. Use --test to run with test files")
        sys.exit(1)

    print("\n--- Running Anisotropy Analysis ---")
    anisotropy_analysis = Anisotropy(
        u,
        residues_group_parsed,
        filePath=args.output_csv,
        parallel=args.parallel,
        n_jobs=args.n_jobs
    )
    # Set verbose for progress bar display (default True in CLI mode)
    anisotropy_analysis._verbose = args.verbose
    anisotropy_analysis.run(
        start=args.start_frame,
        stop=args.stop_frame,
        step=args.step_frame,
        verbose=args.verbose,
        plot_type=args.plot,
        plot_output=args.plot_output
    )

    print("\n--- Analysis Finished ---")
    if args.output_csv:
        print(f"Results saved to: {args.output_csv}")


if __name__ == "__main__":
    main()


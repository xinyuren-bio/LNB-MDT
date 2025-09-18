import os
import sys
import argparse
import ast
import time
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
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
from analysis.parameter_utils import parse_residues_simple, parse_gas_group_simple 

__all__ = ['Density', 'DensityMultiRadius', 'DensityVisualizer']

class Density(AnalysisBase):
    """单半径密度分析类"""
    def __init__(self, universe, ResiudeGroup: dict, GasGroup: dict, MW: float = 14, radius: float = 50, filePath: str = None, parallel: bool = False, n_jobs: int = -1):
        self.radius = radius
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(ResiudeGroup)
        self.gas_group = list(GasGroup)
        self.MW = MW
        self.file_path = filePath
        self.parallel = parallel
        self.n_jobs = n_jobs
        self.headSp = {sp: ' '.join(ResiudeGroup[sp]) for sp in ResiudeGroup}
        self.gas_sp = {sp: ' '.join(GasGroup[sp]) for sp in GasGroup}

        self.headAtoms = self.u.atoms[[]]
        self.gas_atoms = self.u.atoms[[]]
        for sp in self.residues:
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]), updating=False)
        for sp in self.gas_group:
            self.gas_atoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.gas_sp[sp]), updating=False)

        self._n_residues = self.headAtoms.n_residues
        self.results.DensityWithtimes = None

        self.resids = self.headAtoms.resids
        self.resnames = self.headAtoms.resnames
        self.parameters = f"{ResiudeGroup}, {GasGroup}, Radius:{radius}, MW:{MW}"

    def _single_frame(self):
        """计算单帧的密度"""
        # 获取当前帧的原子位置
        head_positions = self.headAtoms.positions
        
        density_values = np.zeros(self._n_residues)
        
        for i, head_pos in enumerate(head_positions):
            # 检查位置是否有效
            if head_pos is None or np.any(np.isnan(head_pos)):
                density_values[i] = 0
                continue
                
            # 使用MDAnalysis的point选择语句
            x, y, z = head_pos
            selection_string = f"point {x:.3f} {y:.3f} {z:.3f} {self.radius:.3f}"
            
            try:
                # 选择半径内的气体原子
                selected_gas = self.gas_atoms.select_atoms(selection_string)
                n_atoms = selected_gas.n_atoms
                
                if n_atoms > 0:
                    # 计算体积 (球体体积)
                    volume = (4/3) * np.pi * (self.radius ** 3) * 1e-30  # 转换为m³
                    
                    # 计算质量 (假设气体分子量为MW)
                    mass = n_atoms * self.MW * 1.66054e-27  # 转换为kg
                    
                    # 计算密度
                    if volume > 0:
                        density_values[i] = mass / volume
                    else:
                        density_values[i] = 0
                else:
                    density_values[i] = 0
            except Exception as e:
                print(f"Warning: Selection failed for residue {i}: {e}")
                density_values[i] = 0
        
        # 存储结果
        if self.results.DensityWithtimes is None:
            self.results.DensityWithtimes = np.zeros((self._n_residues, self.n_frames))
        
        self.results.DensityWithtimes[:, self._frame_index] = density_values

    def run(self, start=None, stop=None, step=None, frames=None,
            verbose=None, *, progressbar_kwargs={}, callBack=None):
        """运行密度分析"""
        verbose = getattr(self, '_verbose', False) if verbose is None else verbose
        
        self._setup_frames(self._trajectory, start=start, stop=stop, step=step, frames=frames)
        self._prepare()
        
        if self.parallel and self.n_jobs != 1:
            print("Running in parallel mode...")
            # 并行处理
            results_list = Parallel(n_jobs=self.n_jobs)(
                delayed(self._process_frame_parallel)(i, ts) 
                for i, ts in enumerate(self._sliced_trajectory)
            )
            if results_list:
                self.results.DensityWithtimes = np.array(results_list).T
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose, callBack=callBack)

        self._conclude()

    def _process_frame_parallel(self, frame_idx, ts):
        """并行处理单帧"""
        # 跳转到指定帧
        self._trajectory[ts.frame]
        
        # 获取当前帧的原子位置
        head_positions = self.headAtoms.positions
        
        density_values = np.zeros(self._n_residues)
        
        for i, head_pos in enumerate(head_positions):
            # 检查位置是否有效
            if head_pos is None or np.any(np.isnan(head_pos)):
                density_values[i] = 0
                continue
                
            # 使用MDAnalysis的point选择语句
            x, y, z = head_pos
            selection_string = f"point {x:.3f} {y:.3f} {z:.3f} {self.radius:.3f}"
            
            try:
                # 选择半径内的气体原子
                selected_gas = self.gas_atoms.select_atoms(selection_string)
                n_atoms = selected_gas.n_atoms
                
                if n_atoms > 0:
                    # 计算体积 (球体体积)
                    volume = (4/3) * np.pi * (self.radius ** 3) * 1e-30
                    
                    # 计算质量 (假设气体分子量为MW)
                    mass = n_atoms * self.MW * 1.66054e-27
                    
                    # 计算密度
                    if volume > 0:
                        density_values[i] = mass / volume
                    else:
                        density_values[i] = 0
                else:
                    density_values[i] = 0
            except Exception as e:
                print(f"Warning: Selection failed for residue {i}: {e}")
                density_values[i] = 0
        
        return density_values

    def _conclude(self):
        if self.file_path:
            # 计算平均密度
            mean_density = np.mean(self.results.DensityWithtimes, axis=0)
            dict_parameter = {
                'frames': list(range(self.start, self.stop, self.step)),
                'results': mean_density,
                'file_path': self.file_path,
                'description': 'Density With Times',
                'parameters': self.parameters,
                'trajectory': self._trajectory
            }
            # Assuming WriteExcelBubble is defined and available
            WriteExcelBubble(**dict_parameter).run()


class DensityMultiRadius(AnalysisBase):
    """多半径密度分析类"""
    def __init__(self, universe, ResiudeGroup: dict, GasGroup: dict, 
                 MW: float = 14, max_radius: float = 50, number_segments: int = 5, 
                 filePath: str = None, parallel: bool = False, n_jobs: int = -1):
        # 生成半径分段
        self.max_radius = max_radius
        self.number_segments = number_segments
        self.radii = self._generate_radius_segments(max_radius, number_segments)
        super().__init__(universe.trajectory)
        self.u = universe
        self.residues = list(ResiudeGroup)
        self.gas_group = list(GasGroup)
        self.MW = MW
        self.file_path = filePath
        self.parallel = parallel
        self.n_jobs = n_jobs
        self.headSp = {sp: ' '.join(ResiudeGroup[sp]) for sp in ResiudeGroup}
        self.gas_sp = {sp: ' '.join(GasGroup[sp]) for sp in GasGroup}

        self.headAtoms = self.u.atoms[[]]
        self.gas_atoms = self.u.atoms[[]]
        for sp in self.residues:
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]), updating=False)
        for sp in self.gas_group:
            self.gas_atoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.gas_sp[sp]), updating=False)

        self._n_residues = self.headAtoms.n_residues
        self.results.DensityMultiRadius = None

        self.resids = self.headAtoms.resids
        self.resnames = self.headAtoms.resnames
        self.parameters = f"{ResiudeGroup}, {GasGroup}, MaxRadius:{max_radius}, Segments:{number_segments}, MW:{MW}"

    def _generate_radius_segments(self, max_radius, number_segments):
        """生成半径分段"""
        return np.linspace(0, max_radius, number_segments + 1)

    def _single_frame(self):
        """计算单帧的多半径密度"""
        # 计算脂质头部的整体质心
        cog = self.headAtoms.center_of_mass()
        
        # 检查质心是否有效
        if cog is None or np.any(np.isnan(cog)):
            density_values = np.zeros(len(self.radii) - 1)
        else:
            density_values = np.zeros(len(self.radii) - 1)
            
            # 使用MDAnalysis的sphlayer选择语句
            x, y, z = cog
            
            for j in range(len(self.radii) - 1):
                inner_radius = self.radii[j]
                outer_radius = self.radii[j + 1]
                
                # 使用逻辑运算符选择环形区域内的原子
                # 语法: (point x y z outer_radius) and not (point x y z inner_radius)
                selection_string = f"(point {x:.3f} {y:.3f} {z:.3f} {outer_radius:.3f}) and not (point {x:.3f} {y:.3f} {z:.3f} {inner_radius:.3f})"
                
                try:
                    # 选择环形区域内的气体原子
                    selected_gas = self.gas_atoms.select_atoms(selection_string)
                    n_atoms = selected_gas.n_atoms
                    
                    if n_atoms > 0:
                        # 计算环形体积
                        inner_vol = (4/3) * np.pi * (inner_radius ** 3) * 1e-30
                        outer_vol = (4/3) * np.pi * (outer_radius ** 3) * 1e-30
                        ring_volume = outer_vol - inner_vol
                        
                        # 计算质量
                        mass = n_atoms * self.MW * 1.66054e-27
                        
                        # 计算密度
                        if ring_volume > 0:
                            density_values[j] = mass / ring_volume
                        else:
                            density_values[j] = 0
                    else:
                        density_values[j] = 0
                except Exception as e:
                    print(f"Warning: Selection failed for radius {j}: {e}")
                    density_values[j] = 0
        
        # 存储结果
        if self.results.DensityMultiRadius is None:
            self.results.DensityMultiRadius = np.zeros((len(self.radii) - 1, self.n_frames))
        
        self.results.DensityMultiRadius[:, self._frame_index] = density_values

    def run(self, start=None, stop=None, step=None, frames=None,
            verbose=None, *, progressbar_kwargs={}, callBack=None):
        """运行多半径密度分析"""
        verbose = getattr(self, '_verbose', False) if verbose is None else verbose
        
        self._setup_frames(self._trajectory, start=start, stop=stop, step=step, frames=frames)
        self._prepare()
        
        if self.parallel and self.n_jobs != 1:
            print("Running in parallel mode...")
            results_list = Parallel(n_jobs=self.n_jobs)(
                delayed(self._process_frame_parallel)(i, ts) 
                for i, ts in enumerate(self._sliced_trajectory)
            )
            if results_list:
                self.results.DensityMultiRadius = np.array(results_list).transpose(1, 0)
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose, callBack=callBack)

        self._conclude()

    def _process_frame_parallel(self, frame_idx, ts):
        """并行处理单帧"""
        self._trajectory[ts.frame]
        
        # 计算脂质头部的整体质心
        cog = self.headAtoms.center_of_mass()
        
        # 检查质心是否有效
        if cog is None or np.any(np.isnan(cog)):
            density_values = np.zeros(len(self.radii) - 1)
        else:
            density_values = np.zeros(len(self.radii) - 1)
            
            # 使用MDAnalysis的sphlayer选择语句
            x, y, z = cog
            
            for j in range(len(self.radii) - 1):
                inner_radius = self.radii[j]
                outer_radius = self.radii[j + 1]
                
                # 使用逻辑运算符选择环形区域内的原子
                # 语法: (point x y z outer_radius) and not (point x y z inner_radius)
                selection_string = f"(point {x:.3f} {y:.3f} {z:.3f} {outer_radius:.3f}) and not (point {x:.3f} {y:.3f} {z:.3f} {inner_radius:.3f})"
                
                try:
                    # 选择环形区域内的气体原子
                    selected_gas = self.gas_atoms.select_atoms(selection_string)
                    n_atoms = selected_gas.n_atoms
                    
                    if n_atoms > 0:
                        # 计算环形体积
                        inner_vol = (4/3) * np.pi * (inner_radius ** 3) * 1e-30
                        outer_vol = (4/3) * np.pi * (outer_radius ** 3) * 1e-30
                        ring_volume = outer_vol - inner_vol
                        
                        # 计算质量
                        mass = n_atoms * self.MW * 1.66054e-27
                        
                        # 计算密度
                        if ring_volume > 0:
                            density_values[j] = mass / ring_volume
                        else:
                            density_values[j] = 0
                    else:
                        density_values[j] = 0
                except Exception as e:
                    print(f"Warning: Selection failed for radius {j}: {e}")
                    density_values[j] = 0
        
        return density_values

    def _conclude(self):
        if self.file_path:
            # 保存多半径密度数据
            self._save_multi_radius_data()

    def _save_multi_radius_data(self):
        """保存多半径密度数据到CSV，使用WriteExcelLipids的格式"""
        # 转换帧数为时间单位
        frames = list(range(self.start, self.stop, self.step))
        time_values, time_unit = WriteExcel._convert_frames_to_time(frames, self._trajectory)
        
        # 确保time_values是数值列表
        if len(time_values) > 0:
            if isinstance(time_values[0], str):
                time_values = [float(t) for t in time_values]
            elif hasattr(time_values[0], 'time'):
                # 如果是Timestep对象，提取time属性
                time_values = [float(t.time) for t in time_values]
        
        # 创建类似WriteExcelLipids的数据框
        # 每个半径层作为一行，包含所有时间点的密度值
        data_rows = []
        
        for radius_idx in range(len(self.radii) - 1):
            inner_radius = self.radii[radius_idx]
            outer_radius = self.radii[radius_idx + 1]
            
            # 构建行数据：半径信息 + 各时间点的密度值
            row_data = {
                'Radius_Layer': f"{inner_radius:.1f}-{outer_radius:.1f}Å",
                'Inner_Radius': inner_radius,
                'Outer_Radius': outer_radius,
                'Volume': (4/3) * np.pi * (outer_radius**3 - inner_radius**3) * 1e-30,  # m³
                **{f"{time_val:.2f}": self.results.DensityMultiRadius[radius_idx, frame_idx]
                   for frame_idx, time_val in enumerate(time_values) if frame_idx < self.results.DensityMultiRadius.shape[1]}
            }
            data_rows.append(row_data)
        
        df_lipid = pd.DataFrame(data_rows)
        
        # 创建时间序列数据框（类似WriteExcelLipids的FRAME部分）
        # 确保数组长度匹配
        n_frames_actual = min(len(time_values), self.results.DensityMultiRadius.shape[1])
        df_frames = pd.DataFrame({
            f'Time ({time_unit})': time_values[:n_frames_actual],
            'Values': np.mean(self.results.DensityMultiRadius[:, :n_frames_actual], axis=0)  # 所有半径层的平均密度
        })
        
        # 保存到CSV，使用WriteExcelLipids的格式
        with open(self.file_path, 'w') as f:
            f.write(f"# Created by LNB-MDT v1.0\n")
            f.write(f"# Multi-Radius Density Analysis\n")
            f.write(f"# TYPE:Lipids\n")
            f.write(f"# Parameters:{self.parameters}\n")
            f.write(f"# Time Unit: {time_unit}\n")
            f.write(f"# Radius Segments: {len(self.radii)-1}\n")
            f.write(f"# Total Frames: {self.n_frames}\n")
            f.write(f"\n")
            
            # 写入半径层数据
            df_lipid.to_csv(f, index=False)
            f.write(f"\n")
            
            # 写入时间序列数据
            f.write(f"# FRAME\n")
            df_frames.to_csv(f, index=False)
        
        print(f"Multi-radius density data saved to: {self.file_path}")


class DensityVisualizer:
    """密度数据可视化类"""
    def __init__(self, csv_file_path):
        """
        初始化可视化器
        
        Args:
            csv_file_path: CSV文件路径
        """
        self.csv_file_path = csv_file_path
        self.df = pd.read_csv(csv_file_path)
        
    def plot_line_chart(self, save_path=None, figsize=(12, 8)):
        """
        绘制密度折线图
        
        Args:
            save_path: 保存图片的路径
            figsize: 图片尺寸
        """
        plt.figure(figsize=figsize)
        
        # 为每个半径绘制一条线
        radii = sorted(self.df['radius'].unique())
        colors = plt.cm.viridis(np.linspace(0, 1, len(radii)))
        
        for i, radius in enumerate(radii):
            radius_data = self.df[self.df['radius'] == radius]
            if i == 0:
                label = f'0-{radius:.1f} Å'
            else:
                prev_radius = radii[i-1]
                label = f'{prev_radius:.1f}-{radius:.1f} Å'
            plt.plot(radius_data['frame'], radius_data['density'], 
                    label=label, color=colors[i], linewidth=2)
        
        plt.xlabel('时间帧', fontsize=12)
        plt.ylabel('密度 (kg/m³)', fontsize=12)
        plt.title('不同半径下气体密度随时间变化', fontsize=14, fontweight='bold')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"折线图已保存到: {save_path}")
        
        plt.show()
    
    def plot_heatmap(self, save_path=None, figsize=(12, 8)):
        """
        绘制密度热力图
        
        Args:
            save_path: 保存图片的路径
            figsize: 图片尺寸
        """
        # 创建透视表
        pivot_table = self.df.pivot(index='radius', columns='frame', values='density')
        
        plt.figure(figsize=figsize)
        sns.heatmap(pivot_table, cmap='viridis', cbar_kws={'label': '密度 (kg/m³)'})
        
        plt.xlabel('时间帧', fontsize=12)
        plt.ylabel('半径范围 (Å)', fontsize=12)
        plt.title('气体密度热力图 (环形区域)', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"热力图已保存到: {save_path}")
        
        plt.show()
    
    def plot_3d_surface(self, save_path=None, figsize=(12, 8)):
        """
        绘制3D表面图
        
        Args:
            save_path: 保存图片的路径
            figsize: 图片尺寸
        """
        # 创建透视表
        pivot_table = self.df.pivot(index='radius', columns='frame', values='density')
        
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
        
        X, Y = np.meshgrid(pivot_table.columns, pivot_table.index)
        Z = pivot_table.values
        
        surf = ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8)
        
        ax.set_xlabel('时间帧', fontsize=12)
        ax.set_ylabel('半径 (Å)', fontsize=12)
        ax.set_zlabel('密度 (kg/m³)', fontsize=12)
        ax.set_title('气体密度3D表面图', fontsize=14, fontweight='bold')
        
        fig.colorbar(surf, shrink=0.5, aspect=5)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"3D表面图已保存到: {save_path}")
        
        plt.show()


# --- Command-line Argument Parsing ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform Density analysis on molecular dynamics trajectories."
    )

    parser.add_argument(
        "--method", "-M",
        type=str,
        choices=['frame', 'radius'],
        default='frame',
        help="Analysis method: 'frame' for single radius, 'radius' for multi-radius analysis."
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
        default="cases/csv/density_results.csv",
        help="Path to the output CSV file for density results."
    )
    parser.add_argument(
        "--residues", "-r",
        type=str,
        default="DPPC:PO4,CHOL:ROH",
        help="Residue groups for analysis. Format: 'DPPC:PO4,CHOL:ROH'"
    )
    parser.add_argument(
        "--gas-group", "-gas",
        type=str,
        default="N2:N2",
        help="Gas group for analysis. Format: 'N2:N2'"
    )
    parser.add_argument(
        "--radius", "-rad",
        type=float,
        default=50.0,
        help="Radius for density calculation (in Angstroms). Used for single radius method."
    )
    parser.add_argument(
        "--max-radius", "-max-rad",
        type=float,
        default=50.0,
        help="Maximum radius for density calculation (in Angstroms). Used for multi-radius method."
    )
    parser.add_argument(
        "--number-segments", "-segments",
        type=int,
        default=5,
        help="Number of radius segments. Used for multi-radius method."
    )
    parser.add_argument(
        "--MW", "-mw",
        type=float,
        default=14.0,
        help="Molecular weight of gas molecules."
    )
    parser.add_argument(
        "--parallel", "-p",
        action="store_true",
        help="Enable parallel processing."
    )
    parser.add_argument(
        "--n-jobs", "-j",
        type=int,
        default=-1,
        help="Number of parallel jobs. -1 means use all available cores."
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
        default=None,
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
        help="Enable verbose output."
    )

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    start_time = time.time()
    
    try:
        residues_group_parsed = parse_residues_simple(args.residues)
        gas_group_parsed = parse_gas_group_simple(args.gas_group)
        
        if not isinstance(residues_group_parsed, dict):
            raise ValueError("Residues argument must be a valid format.")
        if not isinstance(gas_group_parsed, dict):
            raise ValueError("Gas group argument must be a valid format.")
            
    except Exception as e:
        print(f"Error: Could not parse arguments: {e}")
        print("Supported formats:")
        print("  Simple: 'DPPC:PO4,CHOL:ROH' or 'N2:N2'")
        print("  Dictionary: \"{'DPPC': ['PO4']}\"")
        print("  Config file: '@config.json'")
        sys.exit(1)

    print("Selected Residues:")
    for residue, atoms in residues_group_parsed.items():
        print(f"  {residue}: {atoms}")
    print("Selected Gas Group:")
    for residue, atoms in gas_group_parsed.items():
        print(f"  {residue}: {atoms}")

    try:
        u = mda.Universe(args.gro_file, args.xtc_file)
    except Exception as e:
        print(f"Error loading MDAnalysis Universe: {e}")
        print("Please check if GRO/XTC files exist and are valid.")
        sys.exit(1)

    if args.method == 'radius':
        print(f"\nRunning Multi-Radius Density Analysis on {u.trajectory.n_frames} frames...")
        density_analysis = DensityMultiRadius(
            u,
            ResiudeGroup=residues_group_parsed,
            GasGroup=gas_group_parsed,
            filePath=args.output_csv,
            parallel=args.parallel,
            n_jobs=args.n_jobs,
            max_radius=args.max_radius,
            number_segments=args.number_segments,
            MW=args.MW
        )
    else:
        print(f"\nRunning Single-Radius Density Analysis on {u.trajectory.n_frames} frames...")
        density_analysis = Density(
            u,
            ResiudeGroup=residues_group_parsed,
            GasGroup=gas_group_parsed,
            filePath=args.output_csv,
            parallel=args.parallel,
            n_jobs=args.n_jobs,
            radius=args.radius,
            MW=args.MW
        )
    
    density_analysis.run(
        start=args.start_frame,
        stop=args.stop_frame,
        step=args.step_frame,
        verbose=args.verbose
    )

    end_time = time.time()
    elapsed_time = end_time - start_time
    
    print(f"\nAnalysis completed in {elapsed_time:.2f} seconds")
    print(f"Results saved to: {args.output_csv}")

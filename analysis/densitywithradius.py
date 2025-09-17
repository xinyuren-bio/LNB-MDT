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
from joblib import Parallel, delayed

if __name__ == '__main__':
    current_file_path = os.path.abspath(__file__)
    current_dir = os.path.dirname(current_file_path)
    package_root = os.path.abspath(os.path.join(current_dir, '..'))
    if package_root not in sys.path:
        sys.path.insert(0, package_root)
from analysis.analysis_base import * 

__all__ = ['DensityMultiRadius', 'DensityVisualizer']

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
        self.results.DensityWithTimesAndRadii = None
        self.parameters = str(GasGroup) + f'MaxRadius:{self.max_radius}, Segments:{self.number_segments}, MW:{self.MW}'

    def _generate_radius_segments(self, max_radius, number_segments):
        """生成半径分段"""
        segment_size = max_radius / number_segments
        radii = []
        for i in range(1, number_segments + 1):
            radii.append(i * segment_size)
        return radii

    def _prepare(self):
        """初始化结果数组"""
        self.results.DensityWithTimesAndRadii = np.full(
            [self.n_frames, len(self.radii)], fill_value=np.nan
        )

    def _single_frame(self):
        """计算单帧的多半径密度"""
        for i, radius in enumerate(self.radii):
            self.results.DensityWithTimesAndRadii[self._frame_index, i] = \
                self._calculate_density_with_radius(radius)

    def _calculate_density_with_radius(self, radius):
        """计算指定半径环形区域内的气体密度 (kg/m³)"""
        bubble_center = self.headAtoms.center_of_mass()
        
        # 计算当前半径和上一个半径
        segment_index = self.radii.index(radius)
        if segment_index == 0:
            inner_radius = 0
        else:
            inner_radius = self.radii[segment_index - 1]
        
        outer_radius = radius
        
        # 选择环形区域内的气体分子
        # 先选择外半径内的所有分子
        outer_atoms = self.gas_atoms.select_atoms(
            f'point {bubble_center[0]} {bubble_center[1]} {bubble_center[2]} {outer_radius}'
        )
        
        # 再选择内半径内的分子（需要排除）
        if inner_radius > 0:
            inner_atoms = self.gas_atoms.select_atoms(
                f'point {bubble_center[0]} {bubble_center[1]} {bubble_center[2]} {inner_radius}'
            )
            # 计算环形区域内的分子数量
            number_gas = outer_atoms.n_residues - inner_atoms.n_residues
        else:
            number_gas = outer_atoms.n_residues
        
        # 计算环形体积
        outer_volume = 4 / 3 * np.pi * (outer_radius ** 3)
        inner_volume = 4 / 3 * np.pi * (inner_radius ** 3)
        ring_volume_in_A3 = outer_volume - inner_volume
        
        if ring_volume_in_A3 <= 0:
            return 0.0
        
        # 计算密度
        mass_in_grams = number_gas * self.MW / 6.022e23
        density_g_A3 = mass_in_grams / ring_volume_in_A3
        density_kg_m3 = density_g_A3 * 1e27        
        return density_kg_m3

    def run(self, start=None, stop=None, step=None, verbose=None, callBack=None):
        """运行多半径密度分析"""
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
                delayed(self._calculate_frame_densities)() for inputs in inputs_generator
            )
            if results_list:
                self.results.DensityWithTimesAndRadii = np.array(results_list)
        else:
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose, callBack=callBack)

        self._conclude()

    def _calculate_frame_densities(self):
        """计算单帧所有半径的密度"""
        densities = []
        for radius in self.radii:
            densities.append(self._calculate_density_with_radius(radius))
        return densities

    def _conclude(self):
        """保存结果"""
        if self.file_path:
            # 创建DataFrame保存数据
            data = []
            frames = list(range(self.start, self.stop, self.step))
            
            for frame_idx, frame in enumerate(frames):
                for radius_idx, radius in enumerate(self.radii):
                    data.append({
                        'frame': frame,
                        'radius': radius,
                        'density': self.results.DensityWithTimesAndRadii[frame_idx, radius_idx]
                    })
            
            df = pd.DataFrame(data)
            df.to_csv(self.file_path, index=False)
            print(f"Results saved to: {self.file_path}")


class DensityVisualizer:
    """密度数据可视化类"""
    
    def __init__(self, data_file=None, dataframe=None):
        """
        初始化可视化器
        
        Args:
            data_file: CSV文件路径
            dataframe: pandas DataFrame对象
        """
        if dataframe is not None:
            self.df = dataframe
        elif data_file:
            self.df = pd.read_csv(data_file)
        else:
            raise ValueError("必须提供data_file或dataframe参数")
    
    def plot_line_chart(self, save_path=None, figsize=(12, 8)):
        """
        绘制多半径密度随时间变化的折线图
        
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
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
        
        # 准备数据
        frames = sorted(self.df['frame'].unique())
        radii = sorted(self.df['radius'].unique())
        
        X, Y = np.meshgrid(frames, radii)
        Z = np.zeros_like(X)
        
        for i, radius in enumerate(radii):
            for j, frame in enumerate(frames):
                density = self.df[(self.df['radius'] == radius) & 
                                (self.df['frame'] == frame)]['density'].values
                if len(density) > 0:
                    Z[i, j] = density[0]
        
        # 绘制表面
        surf = ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8)
        
        ax.set_xlabel('时间帧', fontsize=12)
        ax.set_ylabel('半径范围 (Å)', fontsize=12)
        ax.set_zlabel('密度 (kg/m³)', fontsize=12)
        ax.set_title('气体密度3D表面图 (环形区域)', fontsize=14, fontweight='bold')
        
        # 添加颜色条
        fig.colorbar(surf, shrink=0.5, aspect=5)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"3D图已保存到: {save_path}")
        
        plt.show()
    
    def plot_all(self, save_dir=None):
        """
        绘制所有类型的图表
        
        Args:
            save_dir: 保存图片的目录
        """
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
        
        # 折线图
        line_path = os.path.join(save_dir, 'density_line_chart.png') if save_dir else None
        self.plot_line_chart(save_path=line_path)
        
        # 热力图
        heatmap_path = os.path.join(save_dir, 'density_heatmap.png') if save_dir else None
        self.plot_heatmap(save_path=heatmap_path)
        
        # 3D图
        surface_path = os.path.join(save_dir, 'density_3d_surface.png') if save_dir else None
        self.plot_3d_surface(save_path=surface_path)


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="Perform Multi-Radius Density Analysis on molecular dynamics trajectories."
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
        default="cases/csv/density_multi_radius_results.csv",
        help="Path to the output CSV file for multi-radius density results."
    )
    parser.add_argument(
        "--residues", "-r",
        type=str,
        default="{'DPPC': ['PO4'], 'DUPC': ['PO4'], 'CHOL': ['ROH']}",
        help="A dictionary string defining residue groups for analysis."
    )
    parser.add_argument(
        "--gas-group", "-a",
        type=str,
        default="{'N2': ['N2']}",
        help="A dictionary string defining gas groups for analysis."
    )
    parser.add_argument(
        "--MW", "-m",
        type=float,
        default=14,
        help="Molecular weight for density calculation."
    )
    parser.add_argument(
        "--max-radius", "-R",
        type=float,
        default=50,
        help="Maximum radius for analysis (Å)."
    )
    parser.add_argument(
        "--number-segments", "-n",
        type=int,
        default=5,
        help="Number of radius segments to divide the analysis into."
    )
    parser.add_argument(
        "--plot-type", "-P",
        type=str,
        choices=['line', 'heatmap', '3d', 'all'],
        default='all',
        help="Type of plot to generate."
    )
    parser.add_argument(
        "--plot-dir", "-d",
        type=str,
        default="cases/plots",
        help="Directory to save plots."
    )
    parser.add_argument(
        "--parallel", "-p",
        action="store_true",
        help="Enable parallel processing."
    )
    parser.add_argument(
        "--n-jobs", "-j",
        type=int,
        default=2,
        help="Number of jobs to run in parallel."
    )
    parser.add_argument(
        "--start-frame", "-s",
        type=int,
        default=0,
        help="Starting frame for analysis."
    )
    parser.add_argument(
        "--stop-frame", "-e",
        type=int,
        help="Stopping frame for analysis."
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

    # 解析参数
    try:
        residues_group_parsed = ast.literal_eval(args.residues)
        gas_group_parsed = ast.literal_eval(args.gas_group)
        
        if not isinstance(residues_group_parsed, dict):
            raise ValueError("Residues argument must be a dictionary string.")
        if not isinstance(gas_group_parsed, dict):
            raise ValueError("Gas group argument must be a dictionary string.")
    except (ValueError, SyntaxError) as e:
        print(f"Error: Could not parse arguments: {e}")
        sys.exit(1)

    # 加载轨迹
    try:
        u = mda.Universe(args.gro_file, args.xtc_file)
    except Exception as e:
        print(f"Error loading MDAnalysis Universe: {e}")
        sys.exit(1)

    print(f"\nRunning Multi-Radius Density Analysis on {u.trajectory.n_frames} frames...")
    print(f"Max radius: {args.max_radius} Å, Number of segments: {args.number_segments}")
    
    # 运行分析
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
    
    density_analysis.run(
        start=args.start_frame,
        stop=args.stop_frame,
        step=args.step_frame,
        verbose=args.verbose
    )

    # 生成可视化
    print("\nGenerating visualizations...")
    visualizer = DensityVisualizer(data_file=args.output_csv)
    
    if args.plot_type == 'all':
        visualizer.plot_all(save_dir=args.plot_dir)
    elif args.plot_type == 'line':
        visualizer.plot_line_chart(save_path=os.path.join(args.plot_dir, 'density_line_chart.png'))
    elif args.plot_type == 'heatmap':
        visualizer.plot_heatmap(save_path=os.path.join(args.plot_dir, 'density_heatmap.png'))
    elif args.plot_type == '3d':
        visualizer.plot_3d_surface(save_path=os.path.join(args.plot_dir, 'density_3d_surface.png'))

    end_time = time.time()
    elapsed_time = end_time - start_time
    
    print(f"\nAnalysis completed in {elapsed_time:.2f} seconds")
    print(f"Results saved to: {args.output_csv}")
    print(f"Plots saved to: {args.plot_dir}")

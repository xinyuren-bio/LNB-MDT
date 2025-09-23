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

class DensityTime(AnalysisBase):
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
        
        # Density分析支持的图表类型
        self.supported_figure_types = ['Line Chart', 'Heatmap']
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
                'description': 'Density With Time',
                'parameters': self.parameters,
                'trajectory': self._trajectory,
                'type_name': 'Density With Time'
            }
            # Assuming WriteExcelBubble is defined and available
            WriteExcelBubble(**dict_parameter).run()
            
            # 准备绘图数据
            self._prepare_plot_data()
        else:
            return self.results.DensityWithtimes
    
    def _prepare_plot_data(self):
        """准备绘图数据，格式化为DataFrame"""
        import pandas as pd
        
        # 创建绘图数据
        frames = list(range(self.start, self.stop, self.step))
        data = []
        
        # Density是整体数据，不是按残基分组
        mean_density = np.mean(self.results.DensityWithtimes, axis=0)
        for j, frame in enumerate(frames):
            row = {
                'Time(ns)': frame * self._trajectory.dt / 1000,  # 转换为ns
                'Value': mean_density[j]
            }
            data.append(row)
        
        self.plot_data = pd.DataFrame(data)
        print(f"Density plot data prepared: {self.plot_data.shape}")
        print(f"Columns: {list(self.plot_data.columns)}")
    
    def plot_line(self, figure_settings=None):
        """绘制Density数据的折线图"""
        if self.plot_data is None:
            self._prepare_plot_data()
        
        if figure_settings is None:
            figure_settings = {
                'x_title': 'Time(ns)',
                'y_title': 'Density',
                'axis_text': 12,
                'marker_shape': 'o',
                'marker_size': 0,
                'bubble_color': 'red'
            }
        
        # 直接使用matplotlib绘制
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        plt.plot(self.plot_data['Time(ns)'], self.plot_data['Value'],
                marker=figure_settings.get('marker_shape', 'o'),
                markersize=figure_settings.get('marker_size', 0),
                color=figure_settings.get('bubble_color', 'red'),
                label='Density')
        
        plt.xlabel(figure_settings.get('x_title', 'Time(ns)'), fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Density'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Density Analysis Results', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    
    def plot_bar(self, figure_settings=None):
        """绘制Density数据的条形图"""
        if self.plot_data is None:
            self._prepare_plot_data()
        
        if figure_settings is None:
            figure_settings = {
                'x_title': 'Time(ns)',
                'y_title': 'Density',
                'axis_text': 12,
                'bar_color': 'skyblue',
                'bar_min': None,
                'bar_max': None
            }
        
        # 直接使用matplotlib绘制
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        plt.bar(self.plot_data['Time(ns)'], self.plot_data['Value'],
               color=figure_settings.get('bar_color', 'skyblue'),
               alpha=0.7,
               label='Density')
        
        plt.xlabel(figure_settings.get('x_title', 'Time(ns)'), fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Density'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Density Analysis Results', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    
    def plot_heatmap(self, figure_settings=None):
        """绘制Density数据的热图（基于时间序列数据）"""
        if self.plot_data is None:
            self._prepare_plot_data()
        
        # 直接使用matplotlib绘制简化的热图
        import matplotlib.pyplot as plt
        import numpy as np
        
        if figure_settings is None:
            figure_settings = {
                'x_title': 'Time(ns)',
                'y_title': 'Density',
                'axis_text': 12,
                'cmap': 'viridis'
            }
        
        # 创建热图数据：将1D数据转换为2D矩阵
        frames = self.plot_data['Time(ns)'].values
        values = self.plot_data['Value'].values
        
        # 创建2D数据矩阵（每行代表一个时间点，每列代表一个密度值）
        # 这里我们创建一个简化的热图，显示密度随时间的变化
        data_matrix = values.reshape(1, -1)  # 1行，多列
        
        plt.figure(figsize=(12, 6))
        
        # 绘制热图
        im = plt.imshow(data_matrix, cmap=figure_settings.get('cmap', 'viridis'), 
                       aspect='auto', interpolation='nearest')
        
        # 设置坐标轴
        plt.xticks(range(len(frames)), frames)
        plt.yticks([0], ['Density'])
        
        # 添加颜色条
        cbar = plt.colorbar(im)
        cbar.set_label('Density Value', fontsize=figure_settings.get('axis_text', 12))
        
        plt.xlabel(figure_settings.get('x_title', 'Time(ns)'), fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Density'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Density Analysis Heatmap', fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.show()
    
    def plot_3d_surface(self, figure_settings=None):
        """绘制Density数据的3D表面图（基于时间序列数据）"""
        if self.plot_data is None:
            self._prepare_plot_data()
        
        # 直接使用matplotlib绘制简化的3D图
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        import numpy as np
        
        if figure_settings is None:
            figure_settings = {
                'x_title': 'Time(ns)',
                'y_title': 'Density',
                'z_title': 'Value',
                'axis_text': 12,
                'cmap': 'viridis'
            }
        
        # 创建3D数据
        frames = self.plot_data['Time(ns)'].values
        values = self.plot_data['Value'].values
        
        # 创建3D表面数据
        X = frames
        Y = np.zeros_like(frames)  # 创建一个平面
        Z = values
        
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # 绘制3D表面
        ax.plot_surface(X, Y, Z, cmap=figure_settings.get('cmap', 'viridis'), 
                       alpha=0.8, linewidth=0, antialiased=True)
        
        # 设置标签
        ax.set_xlabel(figure_settings.get('x_title', 'Time(ns)'), fontsize=figure_settings.get('axis_text', 12))
        ax.set_ylabel(figure_settings.get('y_title', 'Density'), fontsize=figure_settings.get('axis_text', 12))
        ax.set_zlabel(figure_settings.get('z_title', 'Value'), fontsize=figure_settings.get('axis_text', 12))
        ax.set_title('Density Analysis 3D Surface', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        plt.show()


class DensityRadius(AnalysisBase):
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
        
        # DensityMultiRadius分析支持的图表类型
        self.supported_figure_types = ['Line Chart', 'Heatmap']
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
        cog = self.headAtoms.center_of_mass()
        
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
            
            # 准备绘图数据
            self._prepare_plot_data()
        else:
            return self.results.DensityMultiRadius
    
    def _prepare_plot_data(self):
        """准备绘图数据，格式化为DataFrame（按半径层分组）"""
        import pandas as pd
        
        # 创建绘图数据（按半径层分组）
        frames = list(range(self.start, self.stop, self.step))
        data = []
        
        # 按半径层分组
        for i, radius_range in enumerate(self.radii[:-1]):  # 排除最后一个半径
            row = {
                'Radius': f'{radius_range:.1f}-{self.radii[i+1]:.1f}Å',
                'RadiusIndex': i,
                'RadiusRange': f'{radius_range:.1f}-{self.radii[i+1]:.1f}'
            }
            # 添加每个时间帧的该半径层密度数据（不是平均密度）
            for j, frame in enumerate(frames):
                if self.results.DensityMultiRadius is not None:
                    layer_density = self.results.DensityMultiRadius[i, j]  # 第i个半径层，第j帧
                else:
                    layer_density = 0.0
                row[str(frame)] = layer_density
            data.append(row)
        
        self.plot_data = pd.DataFrame(data)
        print(f"DensityMultiRadius plot data prepared: {self.plot_data.shape}")
    
    def plot_line(self, figure_settings=None):
        """绘制DensityMultiRadius数据的折线图：纵坐标是数值，横坐标是时间，每个折线代表每个半径层"""
        if self.plot_data is None:
            self._prepare_plot_data()
        
        if figure_settings is None:
            figure_settings = {
                'x_title': 'Time(ns)',
                'y_title': 'Density Value',
                'axis_text': 12,
                'marker_shape': 'o',
                'marker_size': 0
            }
        
        # 直接使用matplotlib绘制
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        
        # 获取时间列（第4列及以后）
        time_columns = self.plot_data.columns[3:]
        x_axis = time_columns.astype(float) * self._trajectory.dt / 1000  # 转换为ns
        
        # 按半径层分组绘制
        colors = ['darkred', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'cyan', 'magenta']
        
        for i, (_, row) in enumerate(self.plot_data.iterrows()):
            color = colors[i % len(colors)]
            radius_range = row['Radius']
            # 获取该半径层在所有时间点的平均密度值
            density_values = []
            for col in self.plot_data.columns[3:]:  # 从第4列开始是时间列
                try:
                    density_values.append(float(row[col]))
                except (ValueError, TypeError):
                    density_values.append(0.0)  # 如果转换失败，使用0
            
            plt.plot(x_axis, density_values,
                    marker=figure_settings.get('marker_shape', 'o'),
                    markersize=figure_settings.get('marker_size', 0),
                    color=color,
                    label=f'Layer {i+1}: {radius_range}')
        
        plt.xlabel(figure_settings.get('x_title', 'Time(ns)'), fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Density Value'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Multi-Radius Density Analysis - Line Chart', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
    
    
    
    def plot_heatmap(self, figure_settings=None):
        """绘制DensityMultiRadius数据的热图：每个半径层横着占一层，每层根据时间变化颜色"""
        if self.plot_data is None:
            self._prepare_plot_data()
        
        # 直接使用matplotlib绘制热图
        import matplotlib.pyplot as plt
        import numpy as np
        
        if figure_settings is None:
            figure_settings = {
                'x_title': 'Time(ns)',
                'y_title': 'Radius Layer',
                'axis_text': 12,
                'cmap': 'viridis'
            }
        
        # 创建热图数据矩阵
        # 行：半径层，列：时间点
        time_columns = self.plot_data.columns[3:]  # 时间列（第4列及以后）
        x_axis = time_columns.astype(float) * self._trajectory.dt / 1000  # 转换为ns
        
        # 创建数据矩阵：每行是一个半径层，每列是一个时间点
        data_matrix = []
        radius_labels = []
        
        for i, (_, row) in enumerate(self.plot_data.iterrows()):
            radius_labels.append(f'Layer {i+1}: {row["Radius"]}')
            # 获取时间序列数据（跳过前3列：Radius, RadiusIndex, RadiusRange）
            density_values = []
            for col in self.plot_data.columns[3:]:  # 从第4列开始是时间列
                try:
                    density_values.append(float(row[col]))
                except (ValueError, TypeError):
                    density_values.append(0.0)  # 如果转换失败，使用0
            data_matrix.append(density_values)
        
        data_matrix = np.array(data_matrix, dtype=float)
        
        # 为了确保每个半径层占相等的垂直空间，我们需要扩展数据矩阵
        # 每个半径层扩展为多行，这样看起来像横条
        expanded_data_matrix = []
        rows_per_layer = 10  # 每个半径层扩展为10行
        
        for i in range(data_matrix.shape[0]):  # 对每个半径层
            layer_data = data_matrix[i, :]  # 获取该层的所有时间数据
            for _ in range(rows_per_layer):  # 扩展为多行
                expanded_data_matrix.append(layer_data)
        
        expanded_data_matrix = np.array(expanded_data_matrix)
        
        plt.figure(figsize=(12, 8))
        im = plt.imshow(expanded_data_matrix, cmap=figure_settings.get('cmap', 'viridis'), 
                       aspect='auto', interpolation='nearest')
        
        # 设置坐标轴
        plt.xlabel(figure_settings.get('x_title', 'Time(ns)'), fontsize=figure_settings.get('axis_text', 12))
        plt.ylabel(figure_settings.get('y_title', 'Radius Layer'), fontsize=figure_settings.get('axis_text', 12))
        plt.title('Multi-Radius Density Analysis - Heatmap', fontsize=14)
        
        # 设置x轴标签（时间）
        x_ticks = np.linspace(0, len(x_axis)-1, min(10, len(x_axis)))
        plt.xticks(x_ticks, [f'{x_axis[int(i)]:.1f}' for i in x_ticks])
        
        # 设置y轴标签（半径层）- 每个半径层占多行，所以标签位置需要调整
        y_ticks = []
        y_labels = []
        for i in range(len(radius_labels)):
            tick_pos = (i * rows_per_layer) + (rows_per_layer // 2)  # 每个层的中间位置
            y_ticks.append(tick_pos)
            y_labels.append(radius_labels[i])
        
        plt.yticks(y_ticks, y_labels)
        
        # 添加水平分割线，区分不同的半径层
        for i in range(1, len(radius_labels)):
            plt.axhline(y=i * rows_per_layer - 0.5, color='white', linewidth=1)
        
        # 添加颜色条
        cbar = plt.colorbar(im, label='Density Value')
        cbar.ax.tick_params(labelsize=figure_settings.get('axis_text', 12))
        
        plt.tight_layout()
        plt.show()
    

    def _save_multi_radius_data(self):
        """保存多半径密度数据到CSV，前4列保持不变，后面每一列是每一帧的该半径层密度"""
        # 创建半径层数据框
        data_rows = []
        
        for radius_idx in range(len(self.radii) - 1):
            inner_radius = self.radii[radius_idx]
            outer_radius = self.radii[radius_idx + 1]
            
            # 计算环形体积
            inner_vol = (4/3) * np.pi * (inner_radius ** 3) * 1e-30
            outer_vol = (4/3) * np.pi * (outer_radius ** 3) * 1e-30
            ring_volume = outer_vol - inner_vol
            
            # 构建行数据：前4列保持不变
            row_data = {
                'Radius_Layer': radius_idx + 1,
                'Inner_Radius_A': inner_radius,
                'Outer_Radius_A': outer_radius,
                'Radius_Range': f"{inner_radius:.1f}-{outer_radius:.1f}Å"
            }
            
            # 添加每一帧的该半径层密度作为列
            frames = list(range(self.start, self.stop, self.step))
            for frame_idx, frame in enumerate(frames):
                if self.results.DensityMultiRadius is not None:
                    layer_density = self.results.DensityMultiRadius[radius_idx, frame_idx]
                else:
                    layer_density = 0.0
                row_data[str(frame)] = layer_density
            
            data_rows.append(row_data)
        
        df_radius = pd.DataFrame(data_rows)
        
        # 保存到CSV
        with open(self.file_path, 'w') as f:
            f.write(f"# Created by LNB-MDT v1.0\n")
            f.write(f"# Multi-Radius Density Analysis - Radius Layer Density Values\n")
            f.write(f"# Parameters:{self.parameters}\n")
            f.write(f"# TYPE:Density With Radius\n")
            f.write(f"# Total Radius Layers: {len(self.radii)-1}\n")
            f.write(f"# Max Radius: {self.max_radius}Å\n")
            f.write(f"# Number Segments: {self.number_segments}\n")
            f.write(f"\n")
            
            # 写入半径层数据
            df_radius.to_csv(f, index=False)
            f.write(f"\n")

        
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
        density_analysis = DensityRadius(
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
        density_analysis = DensityTime(
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
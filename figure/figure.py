from abc import abstractmethod, ABC

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import AutoLocator
from .density_figure import DensityFigure

__all__ = ['LipidsFigure', 'BubbleFigure', 'DensityFigure', 'read_excel']

TYPE = {
                'Height (nm)': 0
               , 'Sz Order Parameter (chain: sn1 and sn2)': 0  # 添加Sz Order Parameter支持
               , 'Sz Order Parameter (chain: sn1)': 0  # 添加Sz Order Parameter支持
               , 'Sz Order Parameter (chain: sn2)': 0  # 添加Sz Order Parameter支持
               , 'Bubble Sz Order Parameter (chain: sn1 and sn2)': 0  # 添加Sz Order Parameter支持
               , 'Bubble Sz Order Parameter (chain: sn1)': 0  # 添加Sz Order Parameter支持
               , 'Bubble Sz Order Parameter (chain: sn2)': 0  # 添加Sz Order Parameter支持
               , 'Mean Curvature(nm -1)': 0
               , 'Area(nm^2)': 0
               , 'Density With Times': 1  # 添加Density支持
               , 'Density With Time': 1  # 新的Density With Time类型
               , 'Density With Radius': 3  # 新的Density With Radius类型
               , 'DensityTime': 3  # DensityTime类型
               , 'Density Radius': 3  # DensityRadius类型
               , 'Anisotropy': 1
               , 'Gyration(nm)': 1
               , 'Gyration (nm)': 1
               , 'Cluster': 1
               , 'PCA': 1
               , 'RadialDistribution': 2
               , 'merge': 2
               , 'Bubble Height(nm)': 1
               , 'Bubble Height (nm)': 1  # 添加带空格的版本
               , 'Bubble SZ': 1
               , 'Bubble Mean Curvature(nm -1)': 1
               , 'Bubble Area(nm^2)': 1
               , 'Number of Cluster': 1
               , 'Largest Cluster Size': 1
           }

class Figure(ABC):
    """基础图表类"""
    def __init__(self, description, excel_data, figure_settings):
        """
        :param description: # 分析的性质
        :param excel_data: # 读取excel数据结果
        :param figure_settings: # 参数信息
        """
        self.description = description
        self.data_type = TYPE[description]
        self.excel_data = excel_data
        self.figure_settings = figure_settings

    @abstractmethod
    def plot(self):
        pass

    @staticmethod
    def _unsupported_plot():
        print('该结果不支持绘制当前的类型图！')

    def _set_axes(self):
        """设置坐标轴"""
        # 只有当参数存在且有效时才设置坐标轴范围
        if ('x_min' in self.figure_settings and 'x_max' in self.figure_settings and 
            self.figure_settings['x_min'] < self.figure_settings['x_max']):
            plt.xlim(self.figure_settings['x_min'], self.figure_settings['x_max'])
        if ('y_min' in self.figure_settings and 'y_max' in self.figure_settings and 
            self.figure_settings['y_min'] < self.figure_settings['y_max']):
            plt.ylim(self.figure_settings['y_min'], self.figure_settings['y_max'])


class LipidsFigure(Figure):
    """脂质数据图表类 - 包含Line、Bar、Scatter方法"""
    
    def plot(self):
        """默认绘制方法 - 绘制Line图"""
        self.plot_line()
    
    def plot_line(self):
        """绘制脂质数据的折线图"""
        residue_groups = self.excel_data.groupby('Resname')
        # 获取时间轴数据
        time_columns = self.excel_data.columns[3:]
        x_axis = time_columns.astype(float).to_numpy()
        
        for type_name, df_group in residue_groups:
            result = df_group.iloc[:, 3:].mean(axis=0)
            
            # 获取用户设置的颜色，否则使用默认颜色
            user_color = self.figure_settings.get('color', {}).get(type_name)
            if user_color:
                color = user_color
            else:
                # 使用默认颜色
                default_colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
                color = default_colors[list(residue_groups.groups.keys()).index(type_name) % len(default_colors)]
            
            plt.plot(x_axis, result,
                    marker=self.figure_settings.get('marker_shape', 'o'),
                    markersize=self.figure_settings.get('marker_size', 0),
                    color=color,
                    label=type_name)
        
        # 设置轴标题
        x_title = self.figure_settings.get('x_title')
        y_title = self.figure_settings.get('y_title')
        axis_text_size = self.figure_settings.get('axis_text', 12)
        
        if x_title:
            plt.xlabel(x_title, fontsize=axis_text_size)
        if y_title:
            plt.ylabel(y_title, fontsize=axis_text_size)
        
        # 设置y轴范围
        y_min = self.figure_settings.get('y_min')
        y_max = self.figure_settings.get('y_max')
        
        if y_min is not None and y_max is not None and y_min < y_max:
            plt.ylim(y_min, y_max)
        
        # 设置轴刻度字体大小
        if axis_text_size:
            plt.tick_params(axis='both', which='major', labelsize=axis_text_size)
        
        plt.legend()
        self._set_axes()
        plt.gca().xaxis.set_major_locator(AutoLocator())
        plt.show()

    def plot_bar(self):
        """绘制脂质数据的条形图"""
        residue_groups = self.excel_data.groupby('Resname')
        # 获取时间轴数据
        time_columns = self.excel_data.columns[3:]
        x_axis = time_columns.astype(float).to_numpy()
        
        # 为每个脂质类型计算统计数据
        lipid_names = []
        mean_values = []
        error_values = []
        colors = []
        
        # 默认颜色列表，确保不同脂质类型有不同颜色
        default_colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
        
        for i, (type_name, df_group) in enumerate(residue_groups):
            # 计算每个时间列的平均值
            time_means = df_group.iloc[:, 3:].mean(axis=0)  # 每列的平均值
            
            # 计算总体平均值（所有时间列平均值的平均值）
            overall_mean = time_means.mean()
            
            # 计算标准差（用于error bar）
            overall_std = time_means.std()
            
            lipid_names.append(type_name)
            mean_values.append(overall_mean)
            error_values.append(overall_std)
            
            # 优先使用用户设置的颜色，否则使用默认颜色
            user_color = self.figure_settings.get('color', {}).get(type_name)
            if user_color:
                colors.append(user_color)
            else:
                colors.append(default_colors[i % len(default_colors)])
        
        # 创建条形图
        x_positions = range(len(lipid_names))
        error_bar_enabled = self.figure_settings.get('error_deci', False)
        
        bars = plt.bar(x_positions, mean_values, color=colors, 
                      yerr=error_values if error_bar_enabled else None,
                      capsize=5)
        
        # 设置x轴标签
        plt.xticks(x_positions, lipid_names, rotation=45, ha='right')
        
        # 设置轴标题
        x_title = self.figure_settings.get('x_title')
        y_title = self.figure_settings.get('y_title')
        axis_text_size = self.figure_settings.get('axis_text', 12)
        
        if x_title:
            plt.xlabel(x_title, fontsize=axis_text_size)
        if y_title:
            plt.ylabel(y_title, fontsize=axis_text_size)
        
        # 设置y轴范围
        y_min = self.figure_settings.get('y_min')
        y_max = self.figure_settings.get('y_max')
        
        if y_min is not None and y_max is not None and y_min < y_max:
            plt.ylim(y_min, y_max)
        
        # 设置轴刻度字体大小
        if axis_text_size:
            plt.tick_params(axis='both', which='major', labelsize=axis_text_size)
        
        # 添加图例
        handles = [plt.Rectangle((0,0),1,1, color=color) for color in colors]
        plt.legend(handles, lipid_names)
        plt.show()

    def plot_scatter(self):
        """绘制脂质数据的3D散点图"""
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        result = self.excel_data.iloc[:, 3:].mean(axis=1)
        norm, cmap, boundaries = self._normalize_and_color_map(result)

        residue_groups = self.excel_data.groupby('Resname')
        legend_handles = []
        for type_name, df_group in residue_groups:
            positions = df_group['Coordinates'].str.split(',', expand=True).astype(float)
            positions = positions.to_numpy()
            colors = [cmap(norm(val)) for val in df_group.iloc[:, 3:].mean(axis=1)]
            sc = ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
                            c=colors,
                            s=self.figure_settings.get('shape_size', 50),
                            marker=self.figure_settings.get('shape', {}).get(type_name, 'o'),
                            label=type_name)
            legend_handles.append(sc)

        ax.legend(handles=legend_handles, loc='best',
                  fontsize=self.figure_settings.get('grid_size', 12) if self.figure_settings.get('grid_size', 0) > 0 else None)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical', boundaries=boundaries)
        plt.show()

    def plot_map(self):
        """绘制脂质数据的Map图 - 基于您提供的MocartooScatter插值方法"""
        # 解析Coordinates列获取X、Y、Z坐标
        coordinates = self.excel_data['Coordinates'].str.split(',', expand=True).astype(float)
        positions = coordinates.values  # 3D坐标 (x, y, z)
        
        # 根据每个原子对应的类别赋值（如DPPC=1，DUPC=0）
        # 从Resname列获取原子类型，然后根据类型赋值
        resnames = self.excel_data['Resname'].values
        values = []
        for resname in resnames:
            if resname == 'DPPC':
                values.append(1)
            elif resname == 'DUPC':
                values.append(0)
            else:
                # 其他类型可以设置不同的值
                values.append(0.5)  # 默认值
        values = np.array(values)
        
        # 显示原子类型分布
        unique_resnames, counts = np.unique(resnames, return_counts=True)
        
        # 创建Mocartoo风格的球面投影
        try:
            import cartopy.crs as ccrs
            from scipy.interpolate import griddata
            
            # 计算质心
            center_pos = np.mean(positions, axis=0)
            
            # 使用您代码中的旋转逻辑
            v1 = positions[0]  # 选择第一个原子作为参考
            distance = np.linalg.norm(center_pos - v1)
            v2 = np.array([center_pos[0], center_pos[1], center_pos[2] + distance])
            
            # 计算旋转矩阵（完全按照您的代码）
            v1_vec = v1 - center_pos
            v2_vec = v2 - center_pos
            axis = np.cross(v1_vec, v2_vec)
            axis_norm = np.linalg.norm(axis)
            if axis_norm < 1e-9:
                rotation_matrix = np.identity(3)
            else:
                axis = axis / axis_norm
                angle_cos = np.dot(v1_vec, v2_vec) / (np.linalg.norm(v1_vec) * np.linalg.norm(v2_vec))
                angle = np.arccos(np.clip(angle_cos, -1.0, 1.0))
                I = np.identity(3)
                skew = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
                rotation_matrix = I + np.sin(angle) * skew + (1 - np.cos(angle)) * np.dot(skew, skew)
            
            # 旋转坐标
            rotated_vectors = np.dot(positions - center_pos, rotation_matrix.T)
            
            # 转换为球坐标（按照您的代码）
            r = np.linalg.norm(rotated_vectors, axis=1)
            r[r == 0] = 1e-9  # 防止除以零
            lats_rad = np.arcsin(rotated_vectors[:, 2] / r)
            lons_rad = np.arctan2(rotated_vectors[:, 1], rotated_vectors[:, 0])
            lats = np.degrees(lats_rad)
            lons = np.degrees(lons_rad)
            
            # 创建图形和坐标系
            fig = plt.figure(figsize=(12, 6))
            ax = plt.axes(projection=ccrs.Mollweide())
            ax.set_global()
            fig.patch.set_facecolor('white')
            ax.set_facecolor('white')
            ax.gridlines(linestyle='--', alpha=0.5, color='lightgray')
            ax.spines['geo'].set_edgecolor('black')
            ax.spines['geo'].set_linewidth(1.0)
            
            # 定义颜色映射
            color_map = self.figure_settings.get('color_map', 'coolwarm')
            value_min = self.figure_settings.get('value_min', values.min())
            value_max = self.figure_settings.get('value_max', values.max())
            
            # --- 关键区域：插值和边界处理 ---
            
            # 创建插值网格
            # 定义网格密度，数值越大，图像越精细
            grid_lon = np.linspace(-180, 180, 500)
            grid_lat = np.linspace(-90, 90, 250)
            mesh_lon, mesh_lat = np.meshgrid(grid_lon, grid_lat)
            points = np.vstack((lons, lats)).T
            
            # 执行两步插值来填充边界
            # 步骤1: 首先进行线性插值，保证内部区域的平滑性
            grid_values_linear = griddata(points, values, (mesh_lon, mesh_lat), method='linear')
            
            # 步骤2: 找到线性插值失败的区域 (NaNs)，并用最近邻插值填充它们
            # 这是确保边界被完全覆盖的关键步骤
            nan_mask = np.isnan(grid_values_linear)
            grid_values_nearest = griddata(points, values, (mesh_lon[nan_mask], mesh_lat[nan_mask]), method='nearest')
            
            # 将填充后的值赋给原数组中的 NaN 位置，得到最终的网格数据
            grid_values_final = grid_values_linear
            grid_values_final[nan_mask] = grid_values_nearest
            
            # 绘制插值后的填充图
            # 使用 pcolormesh 替换 contourf 来避免 'GeometryCollection' 错误
            mesh = ax.pcolormesh(
                mesh_lon, 
                mesh_lat, 
                grid_values_final, 
                cmap=color_map, 
                transform=ccrs.PlateCarree(),  # 输入数据是标准的经纬度坐标
                shading='auto',
                vmin=value_min,
                vmax=value_max
            )
            
            # 添加颜色条
            color_bar_enabled = self.figure_settings.get('color_bar', True)
            if color_bar_enabled:
                cbar = plt.colorbar(mesh, ax=ax, orientation='horizontal', pad=0.05, shrink=0.7)
                cbar.set_label('Value', fontsize=12)
                # 针对 0/1 的离散值，设置颜色条的刻度
                if np.array_equal(np.unique(values), [0, 1]):
                    cbar.set_ticks([0.25, 0.75])
                    cbar.set_ticklabels(['0', '1'])
            
            # 添加标题
            title = self.figure_settings.get('title', 'Lipid Distribution Map (Mocartoo Projection)')
            ax.set_title(title, fontsize=16, pad=20)
            
            # 添加统计信息
            stats_text = f'Total Atoms: {len(values)}\nDPPC: {np.sum(values == 1)}, DUPC: {np.sum(values == 0)}'
            plt.figtext(0.02, 0.02, stats_text, fontsize=10, 
                       bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9,
                               edgecolor='gray', linewidth=1))
            
            plt.tight_layout()
            plt.show()
            
        except ImportError:
            # 回退到传统的2D投影
            self._plot_map_2d_fallback(positions, values)
    
    def _plot_map_2d_fallback(self, positions, values):
        """2D投影回退方法"""
        # 使用XY平面投影
        x_coords = positions[:, 0]
        y_coords = positions[:, 1]
        
        # 创建网格用于插值
        x_min, x_max = x_coords.min(), x_coords.max()
        y_min, y_max = y_coords.min(), y_coords.max()
        
        xi = np.linspace(x_min, x_max, 100)
        yi = np.linspace(y_min, y_max, 100)
        xi_grid, yi_grid = np.meshgrid(xi, yi)
        
        zi = griddata((x_coords, y_coords), values, (xi_grid, yi_grid), method='cubic')
        
        fig, ax = plt.subplots(figsize=(12, 9))
        
        value_min = self.figure_settings.get('value_min', values.min())
        value_max = self.figure_settings.get('value_max', values.max())
        color_map = self.figure_settings.get('color_map', 'viridis')
        
        im = ax.imshow(zi, extent=[x_min, x_max, y_min, y_max], 
                      origin='lower', cmap=color_map, 
                      vmin=value_min, vmax=value_max)
        
        ax.set_xlabel('X Position (nm)', fontsize=14)
        ax.set_ylabel('Y Position (nm)', fontsize=14)
        ax.set_title('Lipid Distribution Map (2D Projection)', fontsize=16)
        
        if self.figure_settings.get('color_bar', True):
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label('Value', fontsize=12)
        
        plt.tight_layout()
        plt.show()

    def _normalize_and_color_map(self, result):
        """标准化和颜色映射"""
        bar_min = self.figure_settings.get('bar_min')
        bar_max = self.figure_settings.get('bar_max')
        bar_color = self.figure_settings.get('bar_color', 'viridis')
        
        if bar_min is not None and bar_max is not None and bar_min < bar_max:
            norm = Normalize(bar_min, bar_max)
        else:
            norm = Normalize(result.min(), result.max())
        
        cmap = plt.cm.get_cmap(bar_color)
        boundaries = np.linspace(norm.vmin, norm.vmax, 256)
        return norm, cmap, boundaries


class BubbleFigure(Figure):
    """气泡数据图表类 - 包含Line、Bar方法"""
    
    def plot(self):
        """默认绘制方法 - 绘制Line图"""
        self.plot_line()
    
    def plot_line(self):
        """绘制气泡数据的折线图"""
        # 检查数据格式
        if len(self.excel_data.columns) < 2:
            print("错误：数据列数不足，无法绘制图表")
            return
            
        
        # 对于简单的Time-Values格式，使用前两列
        if len(self.excel_data.columns) == 2:
            # 简单格式：Time, Values
            x_axis = self.excel_data.iloc[:, 0].astype(float).to_numpy()
            y_values = self.excel_data.iloc[:, 1].astype(float).to_numpy()
        else:
            # 复杂格式：从第4列开始是时间数据
            time_columns = self.excel_data.columns[3:]
            x_axis = time_columns.astype(float).to_numpy()
            result = self.excel_data.iloc[:, 3:].mean(axis=0)
            y_values = result.to_numpy()
        
        # 使用气泡颜色
        bubble_color = self.figure_settings.get('bubble_color', 'blue')
        
        plt.plot(x_axis, y_values,
                marker=self.figure_settings.get('marker_shape', 'o'),
                markersize=self.figure_settings.get('marker_size', 0),
                color=bubble_color,
                label='Bubble')
        
        
        # 设置轴标题
        x_title = self.figure_settings.get('x_title')
        y_title = self.figure_settings.get('y_title')
        axis_text_size = self.figure_settings.get('axis_text', 12)
        
        if x_title:
            plt.xlabel(x_title, fontsize=axis_text_size)
        if y_title:
            plt.ylabel(y_title, fontsize=axis_text_size)
        
        # 设置y轴范围
        y_min = self.figure_settings.get('y_min')
        y_max = self.figure_settings.get('y_max')
        
        if y_min is not None and y_max is not None and y_min < y_max:
            plt.ylim(y_min, y_max)
        
        # 设置轴刻度字体大小
        if axis_text_size:
            plt.tick_params(axis='both', which='major', labelsize=axis_text_size)
        
        plt.legend()
        self._set_axes()
        plt.gca().xaxis.set_major_locator(AutoLocator())
        
        plt.show()

    def plot_bar(self):
        """绘制气泡数据的条形图"""
        
        # 检查数据格式
        if len(self.excel_data.columns) < 2:
            print("错误：数据列数不足，无法绘制图表")
            return
            
        
        # 对于简单的Time-Values格式，使用前两列
        if len(self.excel_data.columns) == 2:
            # 简单格式：Time, Values
            y_values = self.excel_data.iloc[:, 1].astype(float).to_numpy()
            overall_mean = y_values.mean()
            overall_std = y_values.std()
        else:
            # 复杂格式：从第4列开始是时间数据
            time_columns = self.excel_data.columns[3:]
            result = self.excel_data.iloc[:, 3:].mean(axis=0)
            overall_mean = result.mean()
            overall_std = result.std()
        
        # 使用气泡颜色
        bubble_color = self.figure_settings.get('bubble_color', 'blue')
        
        # 创建条形图
        error_bar_enabled = self.figure_settings.get('error_deci', False)
        
        bars = plt.bar([0], [overall_mean], color=bubble_color,
                      yerr=[overall_std] if error_bar_enabled else None,
                      capsize=5)
        
        
        # 设置x轴标签
        plt.xticks([0], ['Bubble'])
        
        # 设置轴标题
        x_title = self.figure_settings.get('x_title')
        y_title = self.figure_settings.get('y_title')
        axis_text_size = self.figure_settings.get('axis_text', 12)
        
        if x_title:
            plt.xlabel(x_title, fontsize=axis_text_size)
        if y_title:
            plt.ylabel(y_title, fontsize=axis_text_size)
        
        # 设置y轴范围
        y_min = self.figure_settings.get('y_min')
        y_max = self.figure_settings.get('y_max')
        
        if y_min is not None and y_max is not None and y_min < y_max:
            plt.ylim(y_min, y_max)
        
        # 设置轴刻度字体大小
        if axis_text_size:
            plt.tick_params(axis='both', which='major', labelsize=axis_text_size)
        
        # 添加图例
        handles = [plt.Rectangle((0,0),1,1, color=bubble_color)]
        plt.legend(handles, ['Bubble'])
        
        plt.show()


# 工厂类已删除 - 现在直接使用LipidsFigure和BubbleFigure类


def read_excel(file_path):
    """读取Excel文件并返回描述、数据和时间单位"""
    try:
        comments = []
        time_unit = "ns"  # 默认时间单位
        
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    comment = line.strip()[2:]
                    comments.append(comment)
                    
                    if comment.startswith('TIME_UNIT:'):
                        time_unit = comment.split('TIME_UNIT:')[1].strip()
                else:
                    break
        
        df = pd.read_csv(file_path, skiprows=len(comments), header=0)
        valid_comments = comments[1]
        return valid_comments, df, time_unit
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return None, None, "frame"
"""
Density绘图模块 - 独立设计，与其他代码解耦
支持DensityTime和DensityRadius两种类型的Line和Heatmap绘图
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.ticker import AutoLocator
import seaborn as sns
from scipy.interpolate import griddata

class DensityFigure:
    """Density数据图表类 - 支持Line和Heatmap两种绘图方式"""
    
    def __init__(self, description, excel_data, figure_settings):
        """
        :param description: 分析的性质 ('DensityTime' 或 'DensityRadius')
        :param excel_data: 读取的CSV数据
        :param figure_settings: 参数信息
        """
        self.description = description
        self.excel_data = excel_data
        self.figure_settings = figure_settings
        
        # 根据描述确定数据类型
        if ('DensityTime' in description or 'Density With Times' in description or 
            'Density With Time' in description or 'Density Analysis' in description or
            'Time Series' in description):
            self.data_type = 'time'
        elif ('DensityRadius' in description or 'Multi-Radius Density Analysis' in description or 
              'Density With Radius' in description or 'Multi-Radius Density' in description):
            self.data_type = 'radius'
        else:
            # 如果无法从描述判断，尝试从数据列判断
            if 'Time(ns)' in excel_data.columns or 'Time' in excel_data.columns:
                self.data_type = 'time'
            elif 'Radius' in excel_data.columns or 'radius' in excel_data.columns:
                self.data_type = 'radius'
            else:
                self.data_type = 'unknown'
    
    def plot_line(self):
        """绘制Density数据的折线图"""
        if self.data_type == 'time':
            self._plot_time_line()
        elif self.data_type == 'radius':
            self._plot_radius_line()
        else:
            print(f"未知的Density类型: {self.description}")
    
    def plot_bar(self):
        """绘制Density数据的柱状图（仅支持time类型）"""
        if self.data_type == 'time':
            self._plot_time_bar()
        elif self.data_type == 'radius':
            print("Density Radius类型不支持Bar图表，请使用Line或Heatmap")
        else:
            print(f"未知的Density类型: {self.description}")
    
    def plot_heatmap(self):
        """绘制Density数据的热力图（仅支持radius类型）"""
        if self.data_type == 'time':
            print("Density Time类型不支持Heatmap图表，请使用Line或Bar")
            raise ValueError("Heatmap is not supported for Density Time (1D time series). Use Line or Bar instead.")
        elif self.data_type == 'radius':
            self._plot_radius_heatmap()
        else:
            print(f"未知的Density类型: {self.description}")
    
    def _plot_time_line(self):
        """绘制DensityTime的折线图 - 类似Bubble风格"""
        # 检查数据格式
        if len(self.excel_data.columns) < 2:
            print("错误：数据列数不足，无法绘制图表")
            return
        
        # 检查是否有Time和Density列
        if 'Time(ns)' in self.excel_data.columns and 'Density' in self.excel_data.columns:
            # 标准格式：Frame, Time(ns), Density
            x_axis = self.excel_data['Time(ns)'].astype(float).to_numpy()
            y_values = self.excel_data['Density'].astype(float).to_numpy()
        elif 'Time' in self.excel_data.columns and 'Density' in self.excel_data.columns:
            # 变体格式：Frame, Time, Density
            x_axis = self.excel_data['Time'].astype(float).to_numpy()
            y_values = self.excel_data['Density'].astype(float).to_numpy()
        elif len(self.excel_data.columns) == 2:
            # 简单格式：只有两列
            x_axis = self.excel_data.iloc[:, 0].astype(float).to_numpy()
            y_values = self.excel_data.iloc[:, 1].astype(float).to_numpy()
        elif len(self.excel_data.columns) == 3:
            # 三列格式：假设是 Frame, Time, Density
            x_axis = self.excel_data.iloc[:, 1].astype(float).to_numpy()  # 第2列是时间
            y_values = self.excel_data.iloc[:, 2].astype(float).to_numpy()  # 第3列是密度
        else:
            # 复杂格式：从第4列开始是时间数据（多行数据，每行代表一个时间点）
            time_columns = self.excel_data.columns[3:]
            if len(time_columns) > 0:
                x_axis = time_columns.astype(float).to_numpy()
                result = self.excel_data.iloc[:, 3:].mean(axis=0)
                y_values = result.to_numpy()
            else:
                print("错误：无法识别数据格式")
                return
        
        # 使用DensityTime颜色
        density_color = self.figure_settings.get('density_color', '#2E8B57')  # 海绿色
        
        plt.plot(x_axis, y_values,
                marker=self.figure_settings.get('marker_shape', 'o'),
                linestyle=self.figure_settings.get('line_style', '-'),
                linewidth=self.figure_settings.get('line_width', 2),
                markersize=self.figure_settings.get('marker_size', 6),
                color=density_color,
                label=f'{self.description} Density')
        
        plt.xlabel(self.figure_settings.get('x_label', 'Time (ns)'))
        plt.ylabel(self.figure_settings.get('y_label', 'Density'))
        plt.title(self.figure_settings.get('title', f'{self.description} Density Over Time'))
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 设置坐标轴范围
        self._set_axes()
        plt.tight_layout()
        plt.show()
    
    def _plot_time_bar(self):
        """绘制DensityTime的柱状图"""
        # 检查数据格式
        if len(self.excel_data.columns) < 2:
            print("错误：数据列数不足，无法绘制图表")
            return
        
        # 检查是否有Time和Density列
        if 'Time(ns)' in self.excel_data.columns and 'Density' in self.excel_data.columns:
            # 标准格式：Frame, Time(ns), Density
            x_axis = self.excel_data['Time(ns)'].astype(float).to_numpy()
            y_values = self.excel_data['Density'].astype(float).to_numpy()
        elif 'Time' in self.excel_data.columns and 'Density' in self.excel_data.columns:
            # 变体格式：Frame, Time, Density
            x_axis = self.excel_data['Time'].astype(float).to_numpy()
            y_values = self.excel_data['Density'].astype(float).to_numpy()
        elif len(self.excel_data.columns) == 2:
            # 简单格式：只有两列
            x_axis = self.excel_data.iloc[:, 0].astype(float).to_numpy()
            y_values = self.excel_data.iloc[:, 1].astype(float).to_numpy()
        elif len(self.excel_data.columns) == 3:
            # 三列格式：假设是 Frame, Time, Density
            x_axis = self.excel_data.iloc[:, 1].astype(float).to_numpy()  # 第2列是时间
            y_values = self.excel_data.iloc[:, 2].astype(float).to_numpy()  # 第3列是密度
        else:
            # 复杂格式：从第4列开始是时间数据（多行数据，每行代表一个时间点）
            time_columns = self.excel_data.columns[3:]
            if len(time_columns) > 0:
                x_axis = time_columns.astype(float).to_numpy()
                result = self.excel_data.iloc[:, 3:].mean(axis=0)
                y_values = result.to_numpy()
            else:
                print("错误：无法识别数据格式")
                return
        
        # 使用DensityTime颜色
        density_color = self.figure_settings.get('density_color', '#2E8B57')  # 海绿色
        
        # 计算柱状图宽度（基于时间间隔）
        if len(x_axis) > 1:
            bar_width = (x_axis.max() - x_axis.min()) / len(x_axis) * 0.8
        else:
            bar_width = None
        
        plt.bar(x_axis, y_values,
                width=bar_width if bar_width else self.figure_settings.get('bar_width', None),
                color=density_color,
                alpha=self.figure_settings.get('bar_alpha', 0.7),
                label=f'{self.description} Density')
        
        plt.xlabel(self.figure_settings.get('x_label', 'Time (ns)'))
        plt.ylabel(self.figure_settings.get('y_label', 'Density'))
        plt.title(self.figure_settings.get('title', f'{self.description} Density Over Time'))
        plt.legend()
        plt.grid(True, alpha=0.3, axis='y')
        
        # 设置坐标轴范围
        self._set_axes()
        plt.tight_layout()
        plt.show()
    
    def _plot_radius_line(self):
        """绘制DensityRadius的折线图：纵坐标是数值，横坐标是时间，每个折线代表每个半径层"""
        # DensityRadius数据格式：Radius_Layer, Inner_Radius_A, Outer_Radius_A, Radius_Range, 0, 1, 2, 3, ...
        
        # 检查数据格式
        if 'Radius_Layer' in self.excel_data.columns and len(self.excel_data.columns) > 4:
            # 使用DensityRadius的CSV格式：前4列是基本信息，后面是时间帧数据
            plt.figure(figsize=(10, 6))
            
            # 获取时间列（第5列及以后）
            time_columns = self.excel_data.columns[4:]
            x_axis = time_columns.astype(float)  # 假设时间步长为1，或者可以转换为ns
            
            # 按半径层分组绘制
            colors = ['darkred', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'cyan', 'magenta']
            
            for i, (_, row) in enumerate(self.excel_data.iterrows()):
                color = colors[i % len(colors)]
                radius_range = row['Radius_Range']
                
                # 获取该半径层在所有时间点的密度值
                density_values = []
                for col in time_columns:  # 从第5列开始是时间列
                    try:
                        density_values.append(float(row[col]))
                    except (ValueError, TypeError):
                        density_values.append(0.0)  # 如果转换失败，使用0
                
                plt.plot(x_axis, density_values,
                        marker=self.figure_settings.get('marker_shape', 'o'),
                        markersize=self.figure_settings.get('marker_size', 0),
                        color=color,
                        label=f'Layer {i+1}: {radius_range}')
            
            plt.xlabel(self.figure_settings.get('x_title', 'Time(ns)'), fontsize=self.figure_settings.get('axis_text', 12))
            plt.ylabel(self.figure_settings.get('y_title', 'Density Value'), fontsize=self.figure_settings.get('axis_text', 12))
            plt.title('Multi-Radius Density Analysis - Line Chart', fontsize=14)
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.show()
            
        else:
            # 如果没有标准列，尝试其他方式
            print(f"警告：未找到Radius_Layer列或数据格式不正确，尝试使用其他列")
            print(f"可用列：{list(self.excel_data.columns)}")
            
            # 尝试使用前两列创建简单的折线图
            if len(self.excel_data.columns) >= 2:
                x_axis = self.excel_data.iloc[:, 0].astype(float).to_numpy()
                y_values = self.excel_data.iloc[:, 1].astype(float).to_numpy()
                
                plt.plot(x_axis, y_values,
                        marker=self.figure_settings.get('marker_shape', 'o'),
                        linestyle=self.figure_settings.get('line_style', '-'),
                        linewidth=self.figure_settings.get('line_width', 2),
                        markersize=self.figure_settings.get('marker_size', 5),
                        color=self.figure_settings.get('density_color', '#FF6B6B'),
                        label=f'{self.description} Density')
                
                plt.xlabel(self.figure_settings.get('x_title', 'X Axis'))
                plt.ylabel(self.figure_settings.get('y_title', 'Density'))
                plt.title(self.figure_settings.get('title', f'{self.description} Density'))
                plt.legend()
                plt.grid(True, alpha=0.3)
                plt.tight_layout()
                plt.show()
            else:
                print(f"错误：数据格式不正确，无法绘制折线图")
                return
    
    def _plot_time_heatmap(self):
        """绘制DensityTime的热力图"""
        print(f"DEBUG: _plot_time_heatmap called")
        print(f"DEBUG: excel_data shape: {self.excel_data.shape}")
        print(f"DEBUG: excel_data columns: {list(self.excel_data.columns)}")
        print(f"DEBUG: excel_data head:\n{self.excel_data.head()}")
        
        # 准备数据
        if len(self.excel_data.columns) >= 3:
            print(f"DEBUG: Using 3-column format (Time, Radius, Density)")
            # 假设有Time, Radius, Density三列
            time_data = self.excel_data.iloc[:, 0].astype(float)
            radius_data = self.excel_data.iloc[:, 1].astype(float)
            density_data = self.excel_data.iloc[:, 2].astype(float)
            
            print(f"DEBUG: time_data range: [{time_data.min():.3f}, {time_data.max():.3f}]")
            print(f"DEBUG: radius_data range: [{radius_data.min():.3f}, {radius_data.max():.3f}]")
            print(f"DEBUG: density_data range: [{density_data.min():.3f}, {density_data.max():.3f}]")
            
            # 创建网格
            time_unique = np.sort(time_data.unique())
            radius_unique = np.sort(radius_data.unique())
            
            print(f"DEBUG: time_unique count: {len(time_unique)}")
            print(f"DEBUG: radius_unique count: {len(radius_unique)}")
            
            # 创建网格数据
            time_grid, radius_grid = np.meshgrid(time_unique, radius_unique)
            density_grid = griddata((time_data, radius_data), density_data, 
                                  (time_grid, radius_grid), method='linear')
            
            print(f"DEBUG: density_grid shape: {density_grid.shape}")
            print(f"DEBUG: density_grid range: [{np.nanmin(density_grid):.3f}, {np.nanmax(density_grid):.3f}]")
            print(f"DEBUG: density_grid has NaN: {np.isnan(density_grid).any()}")
            
            # 绘制热力图
            plt.figure(figsize=(10, 8))
            
            # 使用自定义颜色映射
            cmap = self.figure_settings.get('cmap', 'viridis')
            print(f"DEBUG: Using colormap: {cmap}")
            if cmap == 'custom_density':
                colors = ['#000080', '#0000FF', '#00FFFF', '#00FF00', '#FFFF00', '#FF8000', '#FF0000']
                cmap = LinearSegmentedColormap.from_list('custom', colors)
                print(f"DEBUG: Using custom colormap")
            
            im = plt.imshow(density_grid, cmap=cmap, aspect='auto', 
                           extent=[time_unique.min(), time_unique.max(), 
                                 radius_unique.min(), radius_unique.max()],
                           origin='lower')
            
            print(f"DEBUG: imshow created with extent: [{time_unique.min():.3f}, {time_unique.max():.3f}, {radius_unique.min():.3f}, {radius_unique.max():.3f}]")
            
            plt.colorbar(im, label='Density')
            plt.xlabel(self.figure_settings.get('x_title', 'Time (ns)'))
            plt.ylabel(self.figure_settings.get('y_title', 'Radius (nm)'))
            plt.title(self.figure_settings.get('title', f'{self.description} Density Heatmap'))
            
            # 添加等高线
            if self.figure_settings.get('show_contours', True):
                print(f"DEBUG: Adding contours")
                contours = plt.contour(time_grid, radius_grid, density_grid, 
                                      levels=self.figure_settings.get('contour_levels', 10),
                                      colors='white', alpha=0.6, linewidths=0.5)
                plt.clabel(contours, inline=True, fontsize=8)
        
        elif len(self.excel_data.columns) == 2:
            print(f"DEBUG: Using 2-column format (Time, Values) - creating time series heatmap")
            # 2列格式：时间序列数据，创建一维热图
            time_data = self.excel_data.iloc[:, 0].astype(float)
            values_data = self.excel_data.iloc[:, 1].astype(float)
            
            print(f"DEBUG: time_data range: [{time_data.min():.3f}, {time_data.max():.3f}]")
            print(f"DEBUG: values_data range: [{values_data.min():.3f}, {values_data.max():.3f}]")
            
            # 创建一维热图：时间作为x轴，密度值作为颜色
            plt.figure(figsize=(12, 6))
            
            # 创建网格数据用于热图显示
            x_grid = np.linspace(time_data.min(), time_data.max(), len(time_data))
            y_grid = np.array([0, 1])  # 只有一行
            
            # 将密度值扩展为2D数组（重复行）
            values_grid = np.tile(values_data.values, (2, 1))
            
            print(f"DEBUG: values_grid shape: {values_grid.shape}")
            print(f"DEBUG: values_grid range: [{np.nanmin(values_grid):.3f}, {np.nanmax(values_grid):.3f}]")
            
            # 使用自定义颜色映射
            cmap = self.figure_settings.get('cmap', 'viridis')
            print(f"DEBUG: Using colormap: {cmap}")
            
            im = plt.imshow(values_grid, cmap=cmap, aspect='auto', 
                           extent=[time_data.min(), time_data.max(), 0, 1],
                           origin='lower')
            
            print(f"DEBUG: imshow created with extent: [{time_data.min():.3f}, {time_data.max():.3f}, 0, 1]")
            
            plt.colorbar(im, label='Density Value')
            plt.xlabel(self.figure_settings.get('x_title', 'Time (ns)'))
            plt.ylabel(self.figure_settings.get('y_title', 'Density'))
            plt.title(self.figure_settings.get('title', f'{self.description} Density Time Series Heatmap'))
            
            # 隐藏y轴刻度
            plt.yticks([])
            
            # 添加数据点
        
        else:
            print(f"DEBUG: Not enough columns for heatmap, columns: {len(self.excel_data.columns)}")
        
        print(f"DEBUG: Time heatmap plot created")
        plt.tight_layout()
        plt.show()
    
    def _plot_radius_heatmap(self):
        """绘制DensityRadius的热图：每个半径层横着占一层，每层根据时间变化颜色"""
        # DensityRadius数据格式：Radius_Layer, Inner_Radius_A, Outer_Radius_A, Radius_Range, 0, 1, 2, 3, ...
        
        # 检查数据格式
        if 'Radius_Layer' in self.excel_data.columns and len(self.excel_data.columns) > 4:
            # 使用DensityRadius的CSV格式：前4列是基本信息，后面是时间帧数据
            plt.figure(figsize=(12, 8))
            
            # 获取时间列（第5列及以后）
            time_columns = self.excel_data.columns[4:]
            x_axis = time_columns.astype(float)  # 假设时间步长为1，或者可以转换为ns
            
            # 创建数据矩阵：每行是一个半径层，每列是一个时间点
            data_matrix = []
            radius_labels = []
            
            for i, (_, row) in enumerate(self.excel_data.iterrows()):
                radius_labels.append(f'Layer {i+1}: {row["Radius_Range"]}')
                
                # 获取时间序列数据
                density_values = []
                for col in time_columns:  # 从第5列开始是时间列
                    try:
                        density_values.append(float(row[col]))
                    except (ValueError, TypeError):
                        density_values.append(0.0)  # 如果转换失败，使用0
                data_matrix.append(density_values)
            
            data_matrix = np.array(data_matrix, dtype=float)
            
            # 为了确保每个半径层占相等的垂直空间，扩展数据矩阵
            # 每个半径层扩展为多行，这样看起来像横条
            expanded_data_matrix = []
            rows_per_layer = 10  # 每个半径层扩展为10行
            
            for i in range(data_matrix.shape[0]):  # 对每个半径层
                layer_data = data_matrix[i, :]  # 获取该层的所有时间数据
                for _ in range(rows_per_layer):  # 扩展为多行
                    expanded_data_matrix.append(layer_data)
            
            expanded_data_matrix = np.array(expanded_data_matrix)
            
            # 绘制热图
            cmap = self.figure_settings.get('cmap', 'viridis')
            im = plt.imshow(expanded_data_matrix, cmap=cmap, aspect='auto', interpolation='nearest')
            
            # 设置坐标轴
            plt.xlabel(self.figure_settings.get('x_title', 'Time(ns)'), fontsize=self.figure_settings.get('axis_text', 12))
            plt.ylabel(self.figure_settings.get('y_title', 'Radius Layer'), fontsize=self.figure_settings.get('axis_text', 12))
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
            cbar.ax.tick_params(labelsize=self.figure_settings.get('axis_text', 12))
            print('heatmap绘制完成')
            plt.tight_layout()
            plt.show()
            
        else:
            # 如果没有标准列，尝试其他方式
            print(f"警告：未找到Radius_Layer列或数据格式不正确，尝试使用其他列")
            print(f"可用列：{list(self.excel_data.columns)}")
            
            # 尝试使用前两列创建简单的热力图
            if len(self.excel_data.columns) >= 2:
                x_data = self.excel_data.iloc[:, 0].astype(float).to_numpy()
                y_data = self.excel_data.iloc[:, 1].astype(float).to_numpy()
                
                # 创建简单的热力图
                data_matrix = y_data.reshape(1, -1)
                
                plt.figure(figsize=(12, 6))
                im = plt.imshow(data_matrix, cmap='viridis', aspect='auto', interpolation='nearest')
                
                plt.xticks(range(len(x_data)), [f'{x:.1f}' for x in x_data])
                plt.yticks([0], ['Value'])
                
                plt.colorbar(im, label='Value')
                plt.xlabel(self.figure_settings.get('x_title', 'X Axis'))
                plt.ylabel(self.figure_settings.get('y_title', 'Y Axis'))
                plt.title(self.figure_settings.get('title', f'{self.description} Heatmap'))
                plt.tight_layout()
                plt.show()
            else:
                print(f"错误：数据格式不正确，无法绘制热力图")
                return
    
    def _set_axes(self):
        """设置坐标轴范围"""
        if ('x_min' in self.figure_settings and 'x_max' in self.figure_settings and 
            self.figure_settings['x_min'] < self.figure_settings['x_max']):
            plt.xlim(self.figure_settings['x_min'], self.figure_settings['x_max'])
        if ('y_min' in self.figure_settings and 'y_max' in self.figure_settings and 
            self.figure_settings['y_min'] < self.figure_settings['y_max']):
            plt.ylim(self.figure_settings['y_min'], self.figure_settings['y_max'])


class DensityParameterReader:
    """Density参数读取器 - 独立设计，与其他代码解耦"""
    
    def __init__(self, ui_instance):
        self.ui = ui_instance
    
    def read_parameters(self):
        """读取Density绘图参数"""
        parameters = {}
        
        # 基础参数
        parameters['title'] = getattr(self.ui, 'editTitle', None) and self.ui.editTitle.text() or f"Density Analysis"
        parameters['x_label'] = getattr(self.ui, 'editXLabel', None) and self.ui.editXLabel.text() or "X Axis"
        parameters['y_label'] = getattr(self.ui, 'editYLabel', None) and self.ui.editYLabel.text() or "Y Axis"
        
        # 线条参数
        parameters['line_width'] = getattr(self.ui, 'spinLineWidth', None) and self.ui.spinLineWidth.value() or 2
        parameters['line_style'] = getattr(self.ui, 'comboLineStyle', None) and self.ui.comboLineStyle.currentText() or '-'
        parameters['marker_shape'] = getattr(self.ui, 'comboMarkerShape', None) and self.ui.comboMarkerShape.currentText() or 'o'
        parameters['marker_size'] = getattr(self.ui, 'spinMarkerSize', None) and self.ui.spinMarkerSize.value() or 6
        
        # 颜色参数
        parameters['density_color'] = getattr(self.ui, 'comboDensityColor', None) and self.ui.comboDensityColor.currentText() or '#2E8B57'
        
        # 热力图参数
        parameters['colormap'] = getattr(self.ui, 'comboColormap', None) and self.ui.comboColormap.currentText() or 'viridis'
        parameters['show_contours'] = getattr(self.ui, 'checkShowContours', None) and self.ui.checkShowContours.isChecked() or True
        parameters['contour_levels'] = getattr(self.ui, 'spinContourLevels', None) and self.ui.spinContourLevels.value() or 10
        
        # 坐标轴范围
        parameters['x_min'] = getattr(self.ui, 'spinXMin', None) and self.ui.spinXMin.value() or None
        parameters['x_max'] = getattr(self.ui, 'spinXMax', None) and self.ui.spinXMax.value() or None
        parameters['y_min'] = getattr(self.ui, 'spinYMin', None) and self.ui.spinYMin.value() or None
        parameters['y_max'] = getattr(self.ui, 'spinYMax', None) and self.ui.spinYMax.value() or None
        
        return parameters


def create_density_figure(description, excel_data, figure_settings):
    """创建Density图表的工厂函数"""
    return DensityFigure(description, excel_data, figure_settings)

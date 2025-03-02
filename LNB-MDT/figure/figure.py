from abc import abstractmethod, ABC

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.interpolate import griddata
# import cartopy.crs as ccrs

__all__ = ['FigureBar', 'FigureLine', 'FigureScatter', 'read_excel']

TYPE = {'Height(nm)': 0
               , 'SZ': 0
               , 'Mean Curvature(nm -1)': 0
               , 'Area(nm^2)': 0
               , 'Anisotropy': 1
               , 'Gyration(nm)': 1
               , 'Cluster': 1
               , 'PCA': 1
               , 'RadialDistribution': 2
               , 'merge': 2
               , 'Bubble Height(nm)': 1
               , 'Bubble SZ': 1
               , 'Bubble Mean Curvature(nm -1)': 1
               , 'Bubble Area(nm^2)': 1
           }
class Figure(ABC):

    def __init__(self, description, excel_data, figure_settings):
        """
        :param description: # 分析的性质
        :param excel_data: # 读取excel数据结果
        :param figure_settings: # 参数信息
        """
        self.data_type = TYPE[description]
        self.excel_data = excel_data
        self.figure_settings = figure_settings

    @abstractmethod
    def plot(self):
        pass

    @staticmethod
    def _unsupported_plot():
        print('该结果不支持绘制当前的类型图！')


class FigureBar(Figure):
    __slots__ = ('data_type', 'excel_data', 'figure_settings')

    def plot(self):
        plot_strategies = {
            0: self._plot_lipids # lipids
            , 1: self._plot_bubble # bubble
            , 2: self._plot_rdf # rdf
            , -1: self._unsupported_plot # no supply format
        }

        plot_action = plot_strategies.get(self.data_type, -1)
        plot_action()

        self._set_axes()
        plt.show()

    def _plot_lipids(self):
        residue_groups = self.excel_data.groupby('Resname')
        for type_name, df_group in residue_groups:
            data_columns = df_group.iloc[:, 3:]
            result_single = data_columns.mean(axis=1)
            result = result_single.mean(axis=0)

            bars = self._create_bar(type_name, result)
            self._plot_error_bars(bars, result_single, result)
            self._plot_bar_values(bars)
            self._plot_trend_lines(bars)

    def _plot_bubble(self):
        result = self.excel_data.iloc[:, 1]
        bars = self._create_bar(1, result.mean())
        plt.xticks([])
        self._plot_error_bars(bars, result)
        self._plot_bar_values(bars)

    def _plot_rdf(self):
        print('功能正在完善中。。。')


    def _create_bar(self, type_name, height):
        return plt.bar(type_name, height, color=self.figure_settings['color'][type_name])

    def _plot_error_bars(self, bars, data, mean_value=None):
        if not self.figure_settings['error_deci']:
            return
        std_dev = np.std(data)

        for bar in bars:
            bar_height = mean_value if mean_value else bar.get_height()
            plt.errorbar(bar.get_x() + bar.get_width() / 2,
                         bar_height,
                         yerr=std_dev,
                         color='black',
                         capsize=5)

    def _plot_bar_values(self, bars):
        if self.figure_settings['up_bar_value'] <= 0:
            return
        for bar in bars:
            bar_height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width() / 2,
                     bar_height,
                     f'{bar_height:.2f}',
                     ha='center',
                     va='bottom',
                     fontsize=self.figure_settings['up_bar_value'])

    def _plot_trend_lines(self, bars):
        if self.figure_settings['trend_size'] <= 0:
            return
        x, y = [], []
        for bar in bars:
            x.append(bar.get_x() + bar.get_width() * 0.5)
            y.append(bar.get_height())
        plt.plot(x, y
                 , color=self.figure_settings['trend_color']
                 , linestyle='--'
                 , linewidth=self.figure_settings['trend_size'])

    def _set_axes(self):
        plt.xlabel(self.figure_settings['x_title'], fontsize=self.figure_settings['axis_text'])
        plt.ylabel(self.figure_settings['y_title'], fontsize=self.figure_settings['axis_text'])
        if self.figure_settings['y_min'] < self.figure_settings['y_max']:
            plt.ylim(self.figure_settings['y_min'], self.figure_settings['y_max'])


class FigureLine(Figure):
    def plot(self):
        plot_strategy = {
            0: self._plot_lipid
            , 1: self._plot_bubble
            , 2: self._plot_rdf
            , -1: self._unsupported_plot
        }
        plot_method = plot_strategy.get(self.data_type, self._unsupported_plot)
        plot_method()
        self._set_axes()
        plt.show()

    def _plot_lipid(self):
        residue_groups = self.excel_data.groupby('Resname')
        for type_name, df_group in residue_groups:
            result = df_group.iloc[:, 3:].mean(axis=0)
            self._plot_line(result, type_name, self.figure_settings['color'].get(type_name))

    def _plot_bubble(self):
        result = self.excel_data.iloc[:, 1]
        self._plot_line(result, None, self.figure_settings['bubble_color'])

    def _plot_rdf(self):
        x_axis = self.excel_data.iloc[:, 0].to_numpy()
        # plt.xlim(x_axis[0], x_axis[-1])
        type_names = self.excel_data.columns.tolist()[1:]
        for i, type_name in enumerate(type_names, start=1):
            result = self.excel_data.iloc[:, i].to_numpy()
            plt.plot(x_axis, result
                     , marker=self.figure_settings['marker_shape']
                     , markersize=float(self.figure_settings['marker_size'])
                     , label=type_name if self.figure_settings['grid_size'] > 0 else None
                     )

    def _plot_line(self, data, label, color):
        plt.plot(data,
                 marker=self.figure_settings['marker_shape'],
                 markersize=float(self.figure_settings['marker_size']),
                 color=color,
                 label=label if self.figure_settings['grid_size'] > 0 else None)

    def _set_axes(self):
        plt.tick_params(axis='both', labelsize=self.figure_settings['axis_scale'])
        plt.xlabel(self.figure_settings['x_title'], fontsize=self.figure_settings['axis_text'])
        plt.ylabel(self.figure_settings['y_title'], fontsize=self.figure_settings['axis_text'])
        if self.figure_settings['grid_size'] > 0:
            plt.legend(fontsize=self.figure_settings['grid_size'])
        if self.figure_settings['x_min'] < self.figure_settings['x_max']:
            plt.xlim(self.figure_settings['x_min'], self.figure_settings['x_max'])
        if self.figure_settings['y_min'] < self.figure_settings['y_max']:
            plt.ylim(self.figure_settings['y_min'], self.figure_settings['y_max'])


class FigureScatter(Figure):
    def plot(self):
        if self.data_type == 0:
            self._plot_3d_scatter()
        else:
            self._unsupported_plot()

    def _plot_3d_scatter(self):
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
                            s=self.figure_settings['shape_size'],
                            marker=self.figure_settings['shape'].get(type_name, 'o'),
                            label=type_name)
            legend_handles.append(sc)

        ax.legend(handles=legend_handles, loc='best',
                  fontsize=self.figure_settings['grid_size'] if self.figure_settings['grid_size'] > 0 else None)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical', boundaries=boundaries)
        plt.show()

    def _normalize_and_color_map(self, result):
        if self.figure_settings['bar_min'] < self.figure_settings['bar_max']:
            norm = Normalize(self.figure_settings['bar_min'], self.figure_settings['bar_max'])
            boundaries = np.linspace(self.figure_settings['bar_min'], self.figure_settings['bar_max'], 256)
        else:
            norm = Normalize(result.min(), result.max())
            boundaries = np.linspace(result.min(), result.max(), 256)
        cmap = plt.get_cmap(self.figure_settings['bar_color'])
        return norm, cmap, boundaries

# class FigureMap(Figure):
#     def plot(self):
#         if self.type == 0:
#             positions = self.results.iloc[1:, 2:5].to_numpy()
#             center_pos = np.mean(positions, axis=0)
#             r = np.mean(np.linalg.norm(positions - center_pos, axis=1))
#
#
#
#
# class Mocartoo:
#     def __init__(self, positions, values, method: str, v1=None):
#         """
#
#         :param v1: 选中的居中的原子
#         :param positions: 所有需要进行转换的原子
#         :param values: 每个原子的值的大小
#         :param method: 插值方法的选择
#         """
#         self.v1 = v1 if v1.any() else positions[0, ]
#         self.positions = positions
#         self.values = values
#         self._method = method
#         self.center_pos = np.mean(self.positions, axis=0)
#         distance = np.linalg.norm(self.center_pos - self.v1)  # 计算中心原子到选中原子的距离
#         self.v2 = np.array([self.center_pos[0], self.center_pos[1], self.center_pos[2] + distance])
#
#     def _compute_rotation_matrix(self):
#         axis = np.cross(self.v1 - self.center_pos, self.v2 - self.center_pos)
#         axis_norm = np.linalg.norm(axis)
#         if axis_norm == 0:
#             raise ValueError("v1 and v2 are collinear, cannot compute rotation axis.")
#         axis = axis / axis_norm  # 计算轴的单位向量
#         angle = np.arccos(
#             np.dot(self.v1 - self.center_pos, self.v2 - self.center_pos)
#             / (np.linalg.norm(self.v1 - self.center_pos) * np.linalg.norm(self.v2 - self.center_pos)))
#         I = np.identity(3)
#         skew = np.array([[0, -axis[2], axis[1]],
#                          [axis[2], 0, -axis[0]],
#                          [-axis[1], axis[0], 0]])
#         rotation_matrix = I + np.sin(angle) * skew + (1 - np.cos(angle)) * np.dot(skew, skew)
#         return rotation_matrix
#
#     def _rotate_points(self):
#         rotation_matrix = self._compute_rotation_matrix()
#         rotated_points = np.dot(self.positions - self.center_pos, rotation_matrix.T) + self.center_pos
#         return rotated_points
#
#     def make_figure(self):
#         rotated_points = self._rotate_points()
#         points = rotated_points
#         # 计算每个点与质心的向量
#         vectors = points - self.center_pos
#         # 将向量转换为球坐标
#         r = np.linalg.norm(vectors, axis=1)  # 这应该接近我们之前计算的半径
#         theta = np.arccos(vectors[:, 1] / r)  # 纬度（与z轴的夹角）
#         phi = np.arctan2(vectors[:, 0], vectors[:, 2])  # 经度（在xy平面上的角度）
#         # 将弧度转换为度
#         lats = np.degrees(theta) - 90
#         lons = np.degrees(phi)
#         index = np.nonzero(self.values)  # ????
#         values = self.values[index]
#         ##############################################
#         # 创建网格点用于插值
#         grid_lon, grid_lat = np.meshgrid(np.linspace(-180, 180, 100), np.linspace(-90, 90, 50))
#         # 使用scipy进行插值
#         grid_values = griddata((lons, lats), values, (grid_lon, grid_lat), method='cubic')
#         grid_values = np.nan_to_num(grid_values)
#         fig = plt.figure(figsize=(10, 5))
#         ax = plt.axes(projection=ccrs.Mollweide())
#         if self._method == 'scatter':
#             sc = ax.scatter(lons, lats, c=values, transform=ccrs.Geodetic(), cmap='RdBu')
#         elif self._method == 'interplot':
#             contour = ax.contourf(grid_lon, grid_lat, grid_values, cmap='RdBu', transform=ccrs.PlateCarree())
#         ax.set_global()
#         plt.show()


# ///////////////////////  Simple Factory   //////////////////////////////////

def read_excel(file_path):
    try:
        # 读取文件的第一行作为描述信息
        comments = []
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    comments.append(line.strip()[2:])  # 去掉#号
                else:
                    break  # 遇到非注释行停

        # 从第二行开始读取文件，第一行（文件中的第二行）作为表头
        df = pd.read_csv(file_path, skiprows=len(comments), header=0)  # 跳过第一行，直接读取表头和数据
        valid_comments = comments[1]
        return valid_comments, df
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return None, None






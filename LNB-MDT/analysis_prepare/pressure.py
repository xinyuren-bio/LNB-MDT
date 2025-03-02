import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from mayavi import mlab
import plotly.graph_objects as go
from matplotlib.colors import ListedColormap

try:
    from .analysis_base import AnalysisBase
except:
    from analysis_base import AnalysisBase
import warnings

warnings.filterwarnings('ignore')

__all__ = ['Pressure']

u_per_A_s2_to_Pa = 1.66053906660e7


class Pressure(AnalysisBase):
    def __init__(self, u, residuesGroup: dict, gas: list, n_bins: int=None, n_circles: int = None, file_path: str=None) -> None:
        super().__init__(u.trajectory)
        self.u = u
        # if n_bins:
        self.n_bins = n_bins
        # if n_circles:
        self.n_circles = n_circles
        self.file_path = file_path
        self.ball_residues: list[str] = [i for i in residuesGroup]
        self.headSp = {sp :residuesGroup[sp][0] for sp in self.ball_residues}
        self.gas_list = gas # resname N2 and name N2
        self.headAtoms = self.u.atoms[[]]
        self.numResidues = {}
        for i, sp in enumerate(self.ball_residues):
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]), updating=False)

        self.results.Pressure = None

    @property
    def Pressure(self):
        return self.results.Pressure

    def _prepare(self):
        """
        准备结果数组，用于存储计算结果。
        """
        try:
            if self.n_bins:
                self.results.Pressure = np.full([self.n_frames, self.n_bins, self.n_bins, self.n_bins], fill_value=np.NaN)
            else:
                self.results.Pressure = np.full([self.n_frames, self.n_circles], fill_value=np.NaN)
        except MemoryError:
            logging.error("内存不足，无法创建结果数组。")
            raise MemoryError("内存不足")

    def _single_frame(self):
        if self.n_bins:
            self.get_bins_pressure()
        else:
            self.get_circles_pressure()

    # def _conclude(self):
    #     if self.n_bins:
    #         self.results.Pressure = self.results.Pressure.mean(axis=0)
    #     else:
    #         presure_mean = self.results.Pressure.mean(axis=0)

    def get_bins_pressure(self):
        min_value = np.min(np.min(self.headAtoms.positions, axis=0))
        max_value = np.max(np.max(self.headAtoms.positions, axis=0))
        # 计算球心和半径
        center = np.array([(min_value + max_value) * 0.5,
                           (min_value + max_value) * 0.5, (min_value + max_value) * 0.5])
        # 计算半径为立方体对角线的一半
        radius = (max_value - min_value) * 0.5
        selected = self.u.select_atoms(f'resname {self.gas_list[0]} and name {self.gas_list[1]}')
        positions = selected.positions - center
        radial_distances = np.linalg.norm(positions, axis=1)
        within_sphere = radial_distances <= radius
        positions = positions[within_sphere]
        masses = selected.masses[within_sphere]

        # 计算每个原子的动能
        try:
            velocities = selected.velocities[within_sphere]
            v_squared = np.sum(velocities ** 2, axis=1)
            kinetic_energy = 0.5 * masses * v_squared
        except:
            kinetic_energy = 0.5 * masses
        # 为每个坐标轴创建bins
        bins_ = np.linspace(-radius, radius, self.n_bins + 1)
        # 计算每个坐标轴的bin中心
        centers_ = (bins_[:-1] + bins_[1:]) * 0.5

        # 计算每个原子所在的bin索引
        digitize_x = np.digitize(positions[:, 0], bins_) - 1
        digitize_y = np.digitize(positions[:, 1], bins_) - 1
        digitize_z = np.digitize(positions[:, 2], bins_) - 1

        # 修正索引，确保在[0, n_bins-1]范围内
        digitize_x = np.clip(digitize_x, 0, self.n_bins - 1)
        digitize_y = np.clip(digitize_y, 0, self.n_bins - 1)
        digitize_z = np.clip(digitize_z, 0, self.n_bins - 1)

        # 计算每个原子的3D网格索引
        bin_indices = digitize_x * self.n_bins * self.n_bins + digitize_y * self.n_bins + digitize_z
        # 总网格数
        total_bins = self.n_bins ** 3

        # 计算每个bin的总动能
        total_kinetic_energy_per_bin = np.bincount(bin_indices, weights=kinetic_energy, minlength=total_bins)

        # 计算每个bin的体积
        V = ((bins_[1] - bins_[0]) ** 3)
        # 生成网格
        XX, YY, ZZ = np.meshgrid(centers_, centers_, centers_, indexing='ij')

        # 计算边界
        x_min, y_min, z_min = np.min(self.headAtoms.positions - center, axis=0)
        x_max, y_max, z_max = np.max(self.headAtoms.positions - center, axis=0)

        # 条件检查
        cond_x = (XX >= x_min) & (XX <= x_max)
        cond_y = (YY >= y_min) & (YY <= y_max)
        cond_z = (ZZ >= z_min) & (ZZ <= z_max)

        # 找到所有条件都满足的布尔数组
        cond_all = np.all([cond_x, cond_y, cond_z], axis=0)
        XX = XX[cond_all]
        YY = YY[cond_all]
        ZZ = ZZ[cond_all]
        cond_all = cond_all.ravel()
        # 使用np.where找到满足所有条件的三维索引
        # 计算压力：P = (2 * K) / (3 * V)
        pressure_per_bin = (2 * total_kinetic_energy_per_bin) / (3 * V)
        pressure_per_bin *= u_per_A_s2_to_Pa
        pressure_per_bin = pressure_per_bin[cond_all]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax_scatter = ax.scatter(XX, YY, ZZ, c=pressure_per_bin, cmap='viridis', alpha=0.5)
        # ax_scatter = ax.scatter(XX, YY, ZZ)
        cbar = plt.colorbar(ax_scatter)
        plt.show()

    def get_circles_pressure(self):
        min_value = np.min(np.min(self.headAtoms.positions, axis=0))
        max_value = np.max(np.max(self.headAtoms.positions, axis=0))

        center = np.array([(min_value + max_value) * 0.5,
                           (min_value + max_value) * 0.5, (min_value + max_value) * 0.5])

        radius = (max_value - min_value) * 0.5

        selected = self.u.select_atoms(f'resname {self.gas_list[0]} and name {self.gas_list[1]}', updating=True)
        self.n_atoms = selected.n_atoms
        positions = selected.positions - center
        radial_distances = np.linalg.norm(positions, axis=1)
        within_sphere = radial_distances <= radius
        masses = selected.masses[within_sphere]
        radial_distances = radial_distances[within_sphere]
        try:
            velocities = selected.velocities[within_sphere]
            kinetic_energy = 0.5 * masses * np.sum(velocities ** 2, axis=1)
        except:
            kinetic_energy = 0.5 * masses
        circles = np.linspace(0, radius, self.n_circles + 1)
        # 计算每个原子的动能
        volumes = (4 / 3) * np.pi * (circles[1:] ** 3 - circles[:-1] ** 3)  # 计算壳层体积
        ids = np.digitize(radial_distances, circles) - 1
        ids = np.clip(ids, 0, self.n_circles - 1)
        # 计算压力：P = (2 * K) / (3 * V)
        # 计算每个 bin 的总动能
        total_kinetic_energy_per_bin = np.bincount(ids, weights=kinetic_energy, minlength=self.n_circles)

        # pressure_per_circle = np.where(volumes > 0, (2 * total_kinetic_energy_per_bin) / (3 * volumes), 0)
        pressure_per_circle = (2 * total_kinetic_energy_per_bin) / (3 * volumes)
        pressure_per_circle *= u_per_A_s2_to_Pa
        centers = (circles[:-1] + circles[1:]) / 2  # 计算 bin 的中心
        # 存储结果
        self.results.Pressure[self._frame_index, :] = pressure_per_circle
        self.plot_concentric_circles_pressure(centers, pressure_per_circle)

    def plot_concentric_circles_pressure(self, centers, pressure_per_bin):
        """
        绘制同心圆压力分布图。

        参数：
        - centers: 径向 bin 的中心位置
        - pressure_per_bin: 每个径向 bin 的压力值
        """
        grid_size = 500  # 定义网格分辨率
        x = np.linspace(-centers[-1], centers[-1], grid_size)
        y = np.linspace(-centers[-1], centers[-1], grid_size)
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X ** 2 + Y ** 2)
        # 初始化压力网格
        pressure_grid = np.zeros_like(R)

        # 分配压力值到网格
        for i in range(self.n_circles):
            inner_radius = centers[i] - (centers[1] - centers[0]) / 2
            outer_radius = centers[i] + (centers[1] - centers[0]) / 2
            mask = (R >= inner_radius) & (R < outer_radius)
            pressure_grid[mask] = pressure_per_bin[i]

        # 定义自定义颜色映射，从红色到紫色
        colors = plt.cm.RdPu(np.linspace(0, 1, self.n_circles))
        custom_cmap = ListedColormap(colors)
        # 绘制等高线图
        plt.figure(figsize=(8, 8))
        contour = plt.contourf(X, Y, pressure_grid, levels=self.n_circles, cmap=custom_cmap, vmin=0, vmax=450000)
        plt.colorbar(contour, label='Pa')
        plt.xlabel('X (Angstrom)')
        plt.ylabel('Y (Angstrom)')
        # plt.title('压力分布等高线图（红至紫）')
        plt.axis('equal')  # 确保x和y轴的比例相同
        plt.text(-40, -70, f'num of n2 :{self.n_atoms}, mean_pressure:{np.mean(pressure_grid)}')
        plt.show()

    def writeExcel(self):
        arr1 = self.r.reshape([-1, 1]).T
        combined_array = np.row_stack((arr1, self.results.Rad))
        columnHead = ['Presure(Pa)']
        df = pd.DataFrame(combined_array.T, columns=columnHead)
        explanation_df = pd.DataFrame(['Radial Distribution'])
        with pd.CsvWriter(self.file_path, engine='openpyxl') as writer:
            explanation_df.to_excel(writer, sheet_name='Sheet1', index=False, header=False)
            df.to_excel(writer, sheet_name='Sheet1', index=False, startrow=1, header=True)


if __name__ == "__main__":
    import MDAnalysis as mda
    # u = mda.Universe("E:/ach.gro")
    # for i in [120, 140, 160, 180, 200]:
    #     u = mda.Universe(fr"E:\test\{i}.gro")
    #     u = mda.Universe(fr"E:\test\{i}.gro")
    u = mda.Universe('E:/ach.gro', 'E:/ach.xtc')
    cls = Pressure(u
                   , {'DPPC': ['NC3'], 'DAPC': ['NC3']}
                   , ['N2', 'N2']
                   , n_circles=10
                   # , n_bins=5
                 # , filePath='E:/excel/qq.csv')
                 )
        # cls = CalRad(u, {'DPPC':['NC3'],'D3PC':['NC3'],'CHOL':['ROH']}, ncircle=50, filePath='E:/excel/rad1.xlsx')
    cls.run(100, 2100, 400)

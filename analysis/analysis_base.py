from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.spatial import KDTree
from MDAnalysis.analysis.base import Results
from MDAnalysis.lib.log import ProgressBar


class AnalysisBase(ABC):
    def __init__(self, trajectory, verbose=False, **kwargs):
        self._trajectory = trajectory
        self._verbose = verbose
        self.results = Results()
        # 支持的图表类型列表，子类可以重写此属性
        self.supported_figure_types = ['Line Chart', 'Bar Chart']

    def _setup_frames(self, trajectory, start=None, stop=None, step=None,
                      frames=None):
        self._trajectory = trajectory
        if frames is not None:
            if not all(opt is None for opt in [start, stop, step]):
                raise ValueError("start/stop/step cannot be combined with "
                                 "frames")
            slicer = frames
        else:
            start, stop, step = trajectory.check_slice_indices(start, stop,
                                                               step)
            slicer = slice(start, stop, step)
        self._sliced_trajectory = trajectory[slicer]
        self.start = start
        self.stop = stop
        self.step = step
        self.n_frames = len(self._sliced_trajectory)
        self.frames = np.zeros(self.n_frames, dtype=int)
        self.times = np.zeros(self.n_frames)

    @abstractmethod
    def _single_frame(self):
        raise NotImplementedError("Only implemented in child classes")

    def _prepare(self):
        pass

    def _conclude(self):
        pass

    def run(self, start=None, stop=None, step=None, frames=None,
            verbose=None, *, progressbar_kwargs={}, callBack=None):

        verbose = getattr(self, '_verbose',
                          False) if verbose is None else verbose

        self._setup_frames(self._trajectory, start=start, stop=stop,
                           step=step, frames=frames)
        self._prepare()

        for i, ts in enumerate(ProgressBar(
                self._sliced_trajectory,
                verbose=verbose,
                **progressbar_kwargs)):
            self._frame_index = i
            self._ts = ts
            self.frames[i] = ts.frame
            self.times[i] = ts.time
            self._single_frame()
            if callBack:
                callBack((i * self.step/(self.stop - self.start) - 0.01) * 100)
        self._conclude()
        # if callBack:
        #     callBack(100)
        return self


def get_normals(k, particles):
    """
    计算粒子集合中每个粒子的法向量。

    通过查找 k 个最近邻，计算邻居点云的协方差矩阵，
    并取最小特征值对应的特征向量作为法向量。

    Args:
        k (int): 用于计算法向量的最近邻数量。
        particles (np.ndarray): 粒子坐标数组，形状为 (N, 3)。

    Returns:
        np.ndarray: 每个粒子的法向量数组，形状为 (N, 3)。
                    如果特征值分解失败，对应法向量可能为 NaN。
    """
    n_particles = particles.shape[0]
    if n_particles == 0:
        return np.empty((0, 3), dtype=particles.dtype)
    if n_particles < k:
        print(f"警告: 粒子数 ({n_particles}) 小于 k ({k})。可能无法可靠计算法线。")
        # 可以选择返回 NaN，或者减少 k 值等处理方式
        k = n_particles # 例如，使用所有可用粒子

    kdtree = KDTree(particles)
    try:
        # 查询最近邻，包含自身
        _, idxs = kdtree.query(particles, k=k)
    except ValueError as e:
         print(f"错误: KDTree 查询失败，可能是 k 值 ({k}) 或粒子数问题: {e}")
         return np.full((n_particles, 3), np.nan)

    # 获取邻居坐标，形状 (N, k, 3)
    neighbor_points = particles[idxs]

    # 计算质心，形状 (N, 1, 3)
    centroids = np.mean(neighbor_points, axis=1, keepdims=True)

    # 中心化邻居点，形状 (N, k, 3)
    centered_points = neighbor_points - centroids

    # 计算协方差矩阵，形状 (N, 3, 3)
    # cov = sum[(p - centroid).T @ (p - centroid)] / k
    # 使用 einsum 高效计算批量协方差
    cov_matrices = np.einsum('nik,nil->nkl', centered_points, centered_points) / k

    # --- 优化的特征值分解 ---
    try:
        # 使用 numpy.linalg.eigh 直接处理堆叠的矩阵 (N, 3, 3)
        # 它利用底层的 BLAS/LAPACK 库，通常是多线程且高度优化的。
        # eigh 返回排序后的特征值（升序）和对应的特征向量。
        eigenvalues, eigenvectors = np.linalg.eigh(cov_matrices)
        # eigenvalues 形状 (N, 3), eigenvectors 形状 (N, 3, 3)

        # 法向量是对应最小特征值的特征向量
        # 由于 eigh 按升序排序，最小特征值在索引 0 处
        normals = eigenvectors[:, :, 0] # 形状 (N, 3)

    except np.linalg.LinAlgError as e:
        # 如果某个矩阵是奇异的或计算中出现数值问题，会抛出 LinAlgError
        print(f"警告: 特征值分解中出现线性代数错误: {e}")
        print("部分法向量可能无法计算，将填充 NaN。")
        # 创建一个全为 NaN 的数组，后续可以根据需要处理
        # 或者尝试找出哪些计算失败了，只填充那些
        normals = np.full((n_particles, 3), np.nan)
        # 更精细的处理：可以尝试逐个计算失败的，但这会慢很多
        # Alternatively, could loop ONLY through failed ones if identifiable

    return normals


@dataclass
class WriteExcel(ABC):

    @abstractmethod
    def run(self):
        pass
    
    def _convert_frames_to_time(self, frames, trajectory=None):
        """
        将帧数转换为时间单位（ns或ps）
        
        Args:
            frames: 帧数列表
            trajectory: MDAnalysis轨迹对象
            
        Returns:
            tuple: (转换后的时间值列表, 时间单位字符串)
        """
        if trajectory is None:
            # 如果没有轨迹对象，返回原始帧数
            return frames, "frame"
        
        try:
            # 获取时间步长
            trajectory[0]  # 跳转到第一帧
            dt = trajectory.ts.dt  # 获取时间步长（ps）
            
            # 计算总时间
            total_time_ps = len(frames) * dt
            
            # 根据总时间选择单位
            if total_time_ps >= 1000:  # >= 10ns
                # 使用ns单位
                time_values = [(frame - frames[0]) * dt / 1000 for frame in frames]
                time_unit = "ns"
            else:
                # 使用ps单位
                time_values = [(frame - frames[0]) * dt for frame in frames]
                time_unit = "ps"
            
            return time_values, time_unit
            
        except Exception as e:
            print(f"Warning: Could not convert frames to time: {e}")
            return frames, "frame"

    @staticmethod
    def _write_to_csv(file_path: str, comments: list, df: pd.DataFrame, type_: str = None):
        try:
            # 修改文件路径以包含类型
            if type_:
                file_path = file_path.replace('.csv', f'_{type_}.csv')

            # 写入注释
            with open(file_path, mode='w', newline='') as f:
                for comment in comments:
                    f.write(f"# {comment}\n")

            # 追加数据框内容
            df.to_csv(file_path, mode='a', index=False, header=True)
            print(f"Data successfully written to {file_path}")
        except PermissionError:
            print(f"Permission denied: {file_path}")

            raise PermissionDeniedError(f"Permission denied: {file_path}")
        except Exception as e:
            print(f"An error occurred: {e}")
            raise


@dataclass
class WriteExcelLipids(WriteExcel):
    parameters: str
    frames: list
    resids: np.ndarray
    resnames: np.ndarray
    positions: np.ndarray[float]
    results: np.ndarray[float]
    file_path: str
    description: str
    lipids_type: dict
    trajectory: object = None  # 添加轨迹对象用于获取时间信息

    def run(self):
        # column_frame = [i * self.step for i in range(self.n_frames)]
        # LIPID
        lipids_ratio = ':'.join(self.lipids_type) + '=' + ':'.join(map(str, self.lipids_type.values()))
        comments = ['Created by LNB-MDT v1.0', self.description, lipids_ratio, 'TYPE:Lipids', 'Parameters:' + self.parameters]
        
        # 转换帧数为时间单位
        time_values, time_unit = self._convert_frames_to_time(self.frames, self.trajectory)
        
        # 添加时间单位信息到注释中
        comments_with_time = comments + [f'TIME_UNIT:{time_unit}']
        
        df_lipid = pd.DataFrame({
            'Resid': self.resids.astype(int),
            'Resname': self.resnames,
            'Coordinates': [f"{x:.2f},{y:.2f},{z:.2f}" for x, y, z in self.positions],
            **{f"{time_val:.2f}": self.results[:, i].round(3)
               for i, time_val in enumerate(time_values)},
        })
        # FRAME
        df_frames = pd.DataFrame(
            {
                f'Time ({time_unit})': time_values,
                'Values': np.mean(self.results, axis=0).round(3)
            }
        )

        # 储存表头的信息
        # 储存每一帧的数据信息，残基号，残基名称等
        self._write_to_csv(self.file_path, comments_with_time, df_lipid)
        comments_frames = ['Created by LNB-MDT v1.0', 'Bubble ' + self.description, lipids_ratio, 'TYPE:Bubble', 'Parameters:' + self.parameters, f'TIME_UNIT:{time_unit}']
        self._write_to_csv(self.file_path, comments_frames, df_frames, 'frames')


@dataclass
class WriteExcelBubble(WriteExcel):
    frames: int
    results: np.ndarray
    file_path: str
    description: str
    parameters: str
    trajectory: object = None  # 添加轨迹对象用于获取时间信息
    def run(self):
        # column_frame = [i * self.step for i in range(self.n_frames)]
        
        # 转换帧数为时间单位
        time_values, time_unit = self._convert_frames_to_time(self.frames, self.trajectory)
        
        df = pd.DataFrame(
            {
                f'Time ({time_unit})': time_values,
                'Values': self.results.round(3)
            }
        )
        comments = ['Created by LNB-MDT v1.0', self.description, 'TYPE:Bubble', 'Parameters:' + self.parameters]
        self._write_to_csv(self.file_path, comments, df)



class PermissionDeniedError(Exception):
    """Custom exception for permission denied errors."""
    pass
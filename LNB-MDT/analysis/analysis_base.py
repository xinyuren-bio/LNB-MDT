from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.spatial import KDTree
from scipy.linalg import eigh
from MDAnalysis.analysis.base import Results
from MDAnalysis.lib.log import ProgressBar


class AnalysisBase(ABC):
    def __init__(self, trajectory, verbose=False, **kwargs):
        self._trajectory = trajectory
        self._verbose = verbose
        self.results = Results()

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
    kdtree = KDTree(particles)
    _, idxs = kdtree.query(particles, k=k)
    neighbor_points = particles[idxs]
    centroids = np.mean(neighbor_points, axis=1, keepdims=True)
    centered_points = neighbor_points - centroids
    # 计算协方差矩阵
    # 使用爱因斯坦求和约定进行批量协方差计算
    # cov_matrix[i] = (centered_points[i].T @ centered_points[i]) / (k)
    cov_matrices = np.einsum('ijk,ijl->ikl', centered_points, centered_points) / k  # (N, 3, 3)
    # 执行批量特征值分解
    _, eigenvectors = zip(*[eigh(cov) for cov in cov_matrices])
    eigenvectors = np.array(eigenvectors)  # (N, 3, 3)
    normals = eigenvectors[:, :, 0]
    return normals


@dataclass
class WriteExcel(ABC):

    @abstractmethod
    def run(self):
        pass

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
    step: int
    n_frames: int
    resids: np.ndarray
    resnames: np.ndarray
    positions: np.ndarray[float]
    results: np.ndarray[float]
    file_path: str
    description: str
    lipids_type: dict

    def run(self):
        column_frame = [i * self.step for i in range(self.n_frames)]
        # LIPID
        lipids_ratio = ':'.join(self.lipids_type) + '=' + ':'.join(map(str, self.lipids_type.values()))
        comments = ['Created by LNB-MDT v1.0', self.description, lipids_ratio, 'TYPE:Lipids', 'Parameters:' + self.parameters]
        df_lipid = pd.DataFrame({
            'Resid': self.resids.astype(int),
            'Resname': self.resnames,
            'Coordinates': [f"{x:.2f},{y:.2f},{z:.2f}" for x, y, z in self.positions],
            **{str(frame): self.results[:, i].round(3)
               for i, frame in enumerate(column_frame)}
        })
        # FRAME
        df_frames = pd.DataFrame(
            {
                'Frames': column_frame
                , 'Values': np.mean(self.results, axis=0).round(3)
            }
        )

        # 储存表头的信息
        # 储存每一帧的数据信息，残基号，残基名称等
        self._write_to_csv(self.file_path, comments, df_lipid)
        comments_frames = ['Created by LNB-MDT v1.0', 'Bubble ' + self.description, lipids_ratio, 'TYPE:Bubble', 'Parameters:' + self.parameters]
        self._write_to_csv(self.file_path, comments_frames, df_frames, 'frames')


@dataclass
class WriteExcelBubble(WriteExcel):
    step: int
    n_frames: int
    results: np.ndarray
    file_path: str
    description: str
    parameters: str
    def run(self):
        column_frame = [i * self.step for i in range(self.n_frames)]
        df = pd.DataFrame(
            {
                'Frames': column_frame
                , 'Values': self.results.round(3)
            }
        )
        comments = ['Created by LNB-MDT v1.0', self.description, 'TYPE:Bubble', 'Parameters:' + self.parameters]
        self._write_to_csv(self.file_path, comments, df)



class PermissionDeniedError(Exception):
    """Custom exception for permission denied errors."""
    pass
from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
import pandas as pd

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


@dataclass
class WriteExcel(ABC):

    @abstractmethod
    def run(self):
        pass

    @staticmethod
    def _write_to_excel(file_path: str, explanation_df: pd.DataFrame, df: pd.DataFrame):
        try:
            with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
                explanation_df.to_excel(writer, sheet_name='Sheet1', index=False, header=False)
                df.to_excel(writer, sheet_name='Sheet1', index=False, startrow=1, header=True)
        except Exception as e:
            print(f"Error writing to Excel: {e}")


@dataclass
class WriteExcelLipids(WriteExcel):
    step: int
    n_frames: int
    resids: np.ndarray
    resnames: np.ndarray
    positions: np.ndarray
    results: np.ndarray
    file_path: str
    description: str
    value_divition: float
    lipids_type: dict

    def run(self):
        column_frame = [i * self.step for i in range(self.n_frames)]
        column_head = ['resid', 'resname', 'X', 'Y', 'Z'] + column_frame

        res_data = np.column_stack(
            (self.resids, self.resnames, self.positions[:, 0], self.positions[:, 1], self.positions[:, 2], self.results/self.value_divition))

        lipids_ratio = ':'.join(self.lipids_type) + '=' + ':'.join(map(str, self.lipids_type.values()))

        # 储存表头的信息
        explanation_df = pd.DataFrame([[self.description] + [lipids_ratio]])
        # 储存每一帧的数据信息，残基号，残基名称等
        df = pd.DataFrame(res_data, columns=column_head)
        self._write_to_excel(self.file_path, explanation_df, df)


@dataclass
class WriteExcelBubble(WriteExcel):
    step: int
    n_frames: int
    results: np.ndarray
    file_path: str
    description: str
    value_divition: float

    def run(self):
        column_frame = [i * self.step for i in range(self.n_frames)]
        column_head = ['Frames', 'Values']
        df = pd.DataFrame(np.column_stack((column_frame
                                           , self.results / self.value_divition))
                          , columns=column_head)
        explanation_df = pd.DataFrame([self.description])
        self._write_to_excel(self.file_path, explanation_df, df)

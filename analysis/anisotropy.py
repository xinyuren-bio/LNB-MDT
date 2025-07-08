import warnings
import os
import sys
warnings.filterwarnings('ignore')

import pandas as pd
from scipy.linalg import eigh
import numpy as np
import MDAnalysis as mda
from joblib import Parallel, delayed

if __name__ == '__main__':
    current_file_path = os.path.abspath(__file__)
    current_dir = os.path.dirname(current_file_path)
    package_root = os.path.abspath(os.path.join(current_dir, '..'))
    if package_root not in sys.path:
        sys.path.insert(0, package_root)
from analysis.analysis_base import *

__all__ = ['Anisotropy']


class Anisotropy(AnalysisBase):
    """
    计算脂质分子在双分子层中的各向异性。

    该类通过计算序参数张量及其特征值来分析脂质分子的取向序。
    """

    def __init__(self, universe: mda.Universe, residueGroup: dict, file_path: str = None, parallel: bool = False, n_jobs: int = -1):
        """
        初始化各向异性分析类。

        Parameters
        ----------
        universe : MDAnalysis.Universe
            包含分子动力学轨迹的 MDAnalysis Universe 对象。
            
        residueGroup : dict
            一个指定每种脂质头部原子名称的字典。
            格式: {'lipid_name': ['head_atom_name1', 'head_atom_name2']}
            示例: {'DPPC': ['PO4'], 'CHOL': ['ROH']}
            
        file_path : str, optional
            用于保存分析结果 (CSV 文件) 的路径。
            如果为 None，结果将不会保存到磁盘。
            
        parallel : bool, optional
            是否开启并行计算。默认为 False。
            
        n_jobs : int, optional
            并行计算时使用的工作核心数。默认为 -1，表示使用所有可用的核心。
        """
        super().__init__(universe.trajectory)
        self.u = universe
        self.file_path = file_path
        self.parallel = parallel
        self.n_jobs = n_jobs
        
        # --- 原子选择部分：已恢复为您的原始方法 ---
        self.residues = list(residueGroup)
        # 将头部原子名称转换为空格分隔的字符串
        self.headSp = {sp: ' '.join(residueGroup[sp]) for sp in residueGroup}
        print("Head atoms:", self.headSp)

        # 初始化一个空的 AtomGroup
        self.headAtoms = self.u.atoms[[]]

        # 循环选择指定脂质类型的头部原子并添加到 AtomGroup
        for i in range(len(self.residues)):
            self.headAtoms += self.u.select_atoms(
                f'resname {self.residues[i]} and name {self.headSp[self.residues[i]]}', 
                updating=False
            )
        # --- 原子选择部分结束 ---
        
        if self.headAtoms.n_atoms == 0:
            raise ValueError("Atom selection is empty. Please check your `residueGroup` dictionary and atomic names.")

        self._n_residues = self.headAtoms.n_residues
        self.resids = self.headAtoms.resids
        self.resnames = self.headAtoms.resnames
        self.results.Anisotropy = None

        self.parameters = str(residueGroup)

    @property
    def Anisotropy(self):
        """返回计算出的各向异性结果。"""
        return self.results.Anisotropy

    def _prepare(self):
        """准备用于存储结果的数组。"""
        self.results.Anisotropy = np.full(self.n_frames, fill_value=np.nan)

    @staticmethod
    def _calculate_anisotropy(positions: np.ndarray) -> float: 
        covariance_matrix = np.cov(positions.T)
        eigenvalues, _ = eigh(covariance_matrix)
        anisotropy = np.sqrt(1 - eigenvalues[0]**2 / eigenvalues[2]**2)
        return anisotropy

    def _single_frame(self):
        positions = self.headAtoms.positions
        self.results.Anisotropy[self._frame_index] = self._calculate_anisotropy(positions)

    def _get_positions_generator(self):
        for ts in self.u.trajectory[self.start:self.stop:self.step]:
            yield self.headAtoms.positions

    def run(self, start=None, stop=None, step=None, verbose=None):
        """
        执行分析。

        根据 `self.parallel` 的值选择串行或并行模式。
        """
        if self.parallel:
            self.start = start if start is not None else 0
            self.stop = stop if stop is not None and stop < self._trajectory.n_frames else self._trajectory.n_frames
            self.step = step if step is not None else 1
            self.n_frames = len(range(self.start, self.stop, self.step))
            self._prepare()
            print(f"Running in parallel on {self.n_jobs} jobs...")
            positions_generator = self._get_positions_generator()
            verbose_level = 10 if verbose else 0
            results_list = Parallel(n_jobs=self.n_jobs, verbose=verbose_level)(
                delayed(self._calculate_anisotropy)(pos) for pos in positions_generator
            )
            self.results.Anisotropy = np.array(results_list)
            self._conclude()
        else:
            # --- 串行执行 (调用父类方法，这部分是正确的) ---
            print("Running in serial mode...")
            super().run(start=start, stop=stop, step=step, verbose=verbose)
    def _conclude(self):
        """分析结束后的收尾工作，如保存文件。"""
        if self.file_path:
            dict_parameter = {
                'frames': list(range(self.start, self.stop, self.step)),
                'results': self.results.Anisotropy,
                'file_path': self.file_path,
                'description': 'Anisotropy',
                'parameters': self.parameters
            }
            WriteExcelBubble(**dict_parameter).run()
            print(f"Results saved to {self.file_path}")


if __name__ == "__main__":
    # --- 示例使用 ---
    # 确保文件路径正确
    gro_file = "cases/lnb.gro"
    xtc_file = "cases/md.xtc"
    csv_file_serial = "cases/csv/anisotropy_serial.csv"
    csv_file_parallel = "cases/csv/anisotropy_parallel.csv"

    u = mda.Universe(gro_file, xtc_file)
    
    residue_group = {'DPPC':['PO4'], 'DUPC':['PO4'], 'CHOL':['ROH']}

    # 1. 串行
    analysis_serial = Anisotropy(u, residue_group, file_path=csv_file_serial, parallel=False)
    analysis_serial.run(start=0, stop=100, verbose=True)
    print("Serial calculation finished.")
    print("Results (first 5 frames):", analysis_serial.Anisotropy[:5])
# 
    print("\n" + "="*50 + "\n")

    # 2. 并行
    print("--- Testing Parallel Calculation ---")
    analysis_parallel = Anisotropy(u, residue_group, file_path=csv_file_parallel, parallel=True, n_jobs=2)
    analysis_parallel.run(start=0, stop=500 , verbose=True)
    print("Parallel calculation finished.")

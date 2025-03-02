import numpy as np
from MDAnalysis.analysis.base import AnalysisBase
import MDAnalysis as mda
import warnings
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
"""
原理：选取CHOL附近任意A以内的分子,
      并且统计这些分子归属于哪些residue,并统计residue数目
使用方法：
    需要参数：
    1.u
    eg：u = mda.Universe("ach.tpr", "ach.xtc")
    注：可以只为单独的gro文件
    2.chol_name->str
    eg:chol_name = "CHOL"
    注：输入选定的残基名称，用来获取他附近任意截断半径内的原子
    3.other_name->list
    eg:other_name = ["DDPC", "DAPC"]
    注：要进行分析的残基组分
    4.cut_off_value->float
    注：截断半径为多少A
    return：
    np.array([len(other_name), nframes])
    用法：
    eg：    
    cp = ChlPf(u, "CHOL", ["DPPC", 'DAPC'], 6)
    cp.run(1000, 2000, 2, verbose=True)
    cp.make_figure()
"""


class ChlPf(AnalysisBase):

    def __init__(self, universe, chol_name:str, other_name:list, cut_off_value:float):
        super().__init__(universe.trajectory)
        self.u = universe
        self._chol = chol_name
        self._other_residue = other_name
        self.results.CPF = None
        self._cut_value = cut_off_value

    def _prepare(self):
        self.results.CPF = np.zeros([2, self.n_frames])
        self.arr = np.zeros([len(self._other_residue)])

    def _single_frame(self):
        select_atoms = self.u.select_atoms('(around %d resname %s) and (resname %s)'
                                           % (self._cut_value, self._chol, " ".join(self._other_residue)))

        for i in self._other_residue:
            n_residue = select_atoms.select_atoms('resname %s' % i).n_residues
            self.arr[self._other_residue.index(i)] = n_residue

        num_all = np.sum(self.arr, axis=0)
        self.arr = self.arr / num_all
        self.results.CPF[:, self._frame_index] = self.arr

    def make_figure(self):
        data = self.results.CPF  # 用随机数填充数组作为示例

        # 绘制第一条曲线，即数组的第一行
        for i in self._other_residue:
            plt.plot(data[self._other_residue.index(i), :], label='%s' % i)
        # 设置图例
        plt.legend()
        # 设置x轴标签和y轴标签
        plt.xlabel('Frames')
        plt.ylabel('Cholesterol Preference')
        # 显示网格
        plt.grid(True)
        # 显示
        plt.show()


# u = mda.Universe("D:/awork/lnb/march/113/1500DA/.gro")
u = mda.Universe("D:/awork/lnb/march/113/1500DA/ach.tpr", "D:/awork/lnb/march/113/1500DA/ach.xtc")
cp = ChlPf(u, "CHOL", ["DPPC", 'DAPC'], 6)
cp.run(1000, 2000, 2, verbose=True)
cp.make_figure()



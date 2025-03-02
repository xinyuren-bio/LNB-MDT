import warnings
warnings.filterwarnings('ignore')

import MDAnalysis as mda
import numpy as np
from scipy.spatial import KDTree
from scipy.linalg import eigh

try:
    from .analysis_base import *
except:
    from analysis_base import *

__all__ = ['SZ']




# suppress some MDAnalysis warnings when writing PDB/GRO files


class SZ(AnalysisBase):
    def __init__(self, universe, residuesGroup: dict, chain:str, k: int = None, path=None):
        """
        format:
        residuesGroup - > { 'DPPC': ['PO4'], 'DAPC' : ['PO4']}
        chain - > 'Chain A'
        """
        super().__init__(universe.trajectory)
        self.u = universe
        self.k = k
        self.file_path = path
        self.residues = list(residuesGroup)
        self.chain = chain

        self.tailSp = {}
        self.headSp = {}

        self.headAtoms = self.u.atoms[[]]
        self.tailAtoms = self.u.atoms[[]]

        for sp in residuesGroup:
            self.headSp[sp] = residuesGroup[sp][0]

            if self.chain == 'Chain A':
                self.tailSp[sp] = '??A'
            elif self.chain == 'Chain B':
                self.tailSp[sp] = '??B'
            else:
                self.tailSp[sp] = '??A ??B'

        # 将选择的头部原子和尾部原子按照残基名称进行排序
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]), updating=False)

            self.tailAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.tailSp[sp]), updating=False)

        # 计算每个残基尾部原子的平均位置
        self._atomTailMean = np.array([np.mean(i.positions, axis=0) for i in self.tailAtoms.split('residue')])
        # 得到每个残基头部原子和尾部平均位置的向量，用于后续进行同一区域的区分
        self._atomTailToHead = self._atomTailMean - self.headAtoms.positions
        # 计算分析的原子中所包含的残疾数目
        # 用来写入EXCEL
        self.resids = self.headAtoms.resids
        self.resnames = self.headAtoms.resnames

        self._n_residues = self.headAtoms.n_residues

        self._headMask = {
            sp :self.resnames == sp for sp in self.residues
        }

        self._tailMask = {
            sp: self.tailAtoms.resnames == sp for sp in self.residues
        }

        self._allAtomsMask = {
            sp: self.u.atoms.resnames == sp for sp in self.residues
        }
        # 获得每种类型的残基数量
        self.numSp = {
            sp: i.n_residues for sp, i in self.headAtoms.groupby('resnames').items()
        }


        self.results.SZ = None

    def _prepare(self):
        self.results.SZ = np.full([self._n_residues, self.n_frames],
                                  fill_value=np.NaN)

    def _single_frame(self):
        head_positions = self.headAtoms.positions
        normals = get_normals(self.k, head_positions)
        def compute(tail_atm, residue):
            tail_positions = tail_atm.positions.reshape([self.numSp[residue], -1, 3])
            # 获得每个残基每条尾链的最后一个原子的位置
            tail_tail = tail_positions[:, -1, :]
            head_to_tail = head_positions[self._headMask[sp]] - tail_tail
            consistent_normals = ensure_consistent_normals(sp_normal, head_to_tail)

            # /////////////////////////////////////////////////////////////////
            # 待优化，不需要每帧都计算链的数目
            chain_num = tail_positions.shape[1]

            vectors = tail_positions[:, :chain_num - 2, :] - tail_positions[:, 2:, :]
            vectors_norm = np.linalg.norm(vectors, axis=2, keepdims=True)  # 保留维度以进行广播
            arr = np.sum(vectors * consistent_normals[:, np.newaxis, :], axis=2)[:,:,np.newaxis] / vectors_norm
            # 计算theta和sz
            theta = np.mean(arr, axis=1)
            sz = (3 * theta ** 2 - 1) / 2
            return sz.reshape(-1)

        for sp in self.residues:
            # 获得每个残基他的法向量
            sp_normal = normals[self._headMask[sp]]
            # 得到需要进行分析的原子
            sp_tail = self.tailAtoms[self._tailMask[sp]]
            # 判断该种类型的残基是否选取了两条链进行分析
            if self.tailSp[sp] == '??A ??B':
                chainList = ['??A', '??B']
                spALL = np.full([self.numSp[sp], 2], fill_value=np.NaN)
                for i, chain in enumerate(chainList):
                    sp_chain = sp_tail.select_atoms('name %s' % chain)
                    spALL[:, i] = compute(sp_chain, sp)
                sz = np.mean(spALL, axis=1)
            else:
                sz = compute(sp_tail, sp)
            self.results.SZ[self._headMask[sp], self._frame_index] = sz

    @property
    def SZ(self):
        return self.results.SZ

    def _conclude(self):
        if self.file_path:
            lipids_ratio = {sp: self.u.select_atoms(f'resname {sp}').n_residues for sp in self.residues}
            dict_parameter = {'step': self.step, 'n_frames': self.n_frames, 'resids': self.resids,
                              'resnames': self.resnames,
                              'positions': self.headAtoms.positions, 'results': self.results.SZ,
                              'file_path': self.file_path, 'description': 'SZ',
                              'value_divition': 1, 'lipids_type': lipids_ratio}
            WriteExcelLipids(**dict_parameter).run()



def ensure_consistent_normals(normals, vector):
    multi = normals * vector
    multi_naga = multi[:, ] < 0
    normals[multi_naga] = -normals[multi_naga]
    return normals


if __name__ == "__main__":
    u = mda.Universe("E:/ach.gro", 'E:/ach.xtc')
    cls2 = SZ(u, {'DPPC': ['PO4'], 'DAPC': ['PO4']}, 'Chain A', k=14, path='E:/untitled2.csv')
    cls2.run(0, 100, verbose=True)


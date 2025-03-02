import numpy as np
from MDAnalysis.analysis.base import AnalysisBase
import MDAnalysis as mda
import warnings
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')


class Flip(AnalysisBase):
    def __init__(self, u, LNBSheel: list, CHOLName: str, CHOLAtom:list):
        super().__init__(u.trajectory)
        self.u = u
        self.LNBSheel = self.u.select_atoms('resname %s' % ' '.join(LNBSheel))
        self.cholName = CHOLName
        self.headCHOL = self.u.select_atoms('resname %s and name %s' % (self.cholName, CHOLAtom[0]))
        self.tailCHOL = self.u.select_atoms('resname %s and name %s' % (self.cholName, CHOLAtom[1]))
        self.numCHOL = self.u.select_atoms('resname %s' % self.cholName).n_residues
        self.filpRatio = None

    def _prepare(self):
        self.filpRatio = np.zeros([self.n_frames])
        self.flipYes = np.zeros([self.n_frames])
        self.flipNo = np.zeros([self.n_frames])
        self.flipArray = np.full([self.numCHOL], 1)

    @property
    def flipRatio(self):
        return self.filpRatio

    def _single_frame(self):
        head = self.headCHOL.positions
        tail = self.tailCHOL.positions
        center = self.LNBSheel.center_of_mass()
        headToTail = head - tail
        headToCenter = head - center
        dotProduct = (headToCenter * headToTail).sum(axis=1)
        originFlipArray = self.flipArray.copy()
        self.flipArray[dotProduct < 0] = -1
        self.flipArray[dotProduct > 0] = 1
        self.flipYes[self._frame_index] = np.sum(self.flipArray[dotProduct > 0])
        self.flipNo[self._frame_index] = np.sum(self.flipArray == -1)
        if np.all(originFlipArray == self.flipArray):
            self.filpRatio[self._frame_index] = 0
        else:
            self.filpRatio[self._frame_index] = np.sum(originFlipArray != self.flipArray)

    def _conclude(self):
        print(self.filpRatio)

if __name__ == "__main__":
    u = mda.Universe("E:/awork/lnb/june/01/B/ach.gro", "E:/awork/lnb/june/01/B/lnb_nojump.xtc")
    cls2 = Flip(u, ['DPPC', 'D3PC', 'CHOL'], 'CHOL', ['ROH','R5'])
    cls2.run(100,-1,1)
    import matplotlib.pyplot as plt
    plt.plot(cls2.filpRatio)
    # plt.plot(cls2.flipYes)
    plt.plot(cls2.flipNo)
    plt.show()
    # cls2.run()
    # cls2.make_figure_2d('bar')




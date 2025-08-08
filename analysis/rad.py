import warnings
warnings.filterwarnings('ignore')
import logging
import argparse  # 添加argparse导入
import ast       # 用于安全地解析residues_group字典字符串
import sys       # 用于sys.exit

import numpy as np
import pandas as pd

try:
    from .analysis_base import *
except:
    from analysis_base import *

__all__ = ['CalRad']


class CalRad(AnalysisBase):
    def __init__(self, u, residuesGroup: dict, ncircle: int=50, filePath: str=None) -> None:
        super().__init__(u.trajectory)
        self.u = u
        self.ncircle: int = ncircle
        self.filePath = filePath
        self.residues: list[str] = list(residuesGroup)
        self.headSp = {sp :residuesGroup[sp][0] for sp in self.residues}
        self.boxHalf = self.u.dimensions[:3] / 2
        self.boxMinHalf = np.min(self.u.dimensions[:3] / 2)
        self.dr = self.boxMinHalf / self.ncircle / 10
        self.r = np.linspace(0, self.boxMinHalf, self.ncircle) / 10
        self.headAtoms = self.u.atoms[[]]
        self.numResidues = {}

        for i, sp in enumerate(self.residues):
            self.headAtoms += self.u.select_atoms('resname %s and name %s'
                                                  % (sp, self.headSp[sp]), updating=False)

        self._residueMask = {sp: self.headAtoms.atoms.resnames == sp for sp in self.residues}
        self.numResidues = {sp: self.u.select_atoms('resname %s' % sp, updating=False).n_residues for sp in self.residues}
        self.results.Rad = None

    @property
    def RAD(self):
        return self.results.Rad

    def _prepare(self):
        """
        准备结果数组，用于存储计算结果。
        """
        try:
            # 创建结果数组
            self.results.Rad = np.zeros([len(self.residues), self.ncircle, self.n_frames])
        except MemoryError:
            logging.error("内存不足，无法创建结果数组。")
            raise MemoryError("内存不足")

    def _single_frame(self):
        center_of_mass = self.u.select_atoms('all').center_of_mass()
        head_positions = self.headAtoms.positions
        try:
            distances_to_center = np.linalg.norm(head_positions - center_of_mass, axis=1)
        except ValueError as e:
            print(f"Error: {e}")
            return
        if head_positions.shape[1] != 3 or center_of_mass.shape[0] != 3:
            raise ValueError("Dimensions of positions and center must be consistent.")

        for residue_index, residue in enumerate(self.residues):
            residue_mask = self._residueMask[residue]
            distances = distances_to_center[residue_mask] / 10

            for circle_index in range(self.ncircle):
                lower_bound = circle_index * self.dr
                upper_bound = (circle_index + 1) * self.dr
                num_atoms = np.sum((distances >= lower_bound) & (distances < upper_bound))
                circle_volume = 4 / 3 * np.pi * (upper_bound ** 3 - lower_bound ** 3)
                self.results.Rad[residue_index, circle_index, self._frame_index] = num_atoms / circle_volume / self.numResidues[residue]

    def _conclude(self):
        self.results.Rad = np.mean(self.results.Rad, axis=-1)
        self.results.Rad.reshape([len(self.residues), -1])
        self.writeExcel()

    def writeExcel(self):
        arr1 = self.r.reshape([-1, 1]).T
        combined_array = np.row_stack((arr1, self.results.Rad))
        columnHead = ['Distance'] + [sp for sp in self.residues]
        df = pd.DataFrame(combined_array.T, columns=columnHead)
        explanation_df = pd.DataFrame(['Radial Distribution'])
        with pd.ExcelWriter(self.filePath, engine='openpyxl') as writer:
            explanation_df.to_excel(writer, sheet_name='Sheet1', index=False, header=False)
            df.to_excel(writer, sheet_name='Sheet1', index=False, startrow=1, header=True)


# --- Command-line Argument Parsing ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform Radial Distribution analysis on molecular dynamics trajectories."
    )

    parser.add_argument(
        "--gro-file",
        type=str,
        default="cases/lnb.gro",
        help="Path to the GRO file (topology file)."
    )
    parser.add_argument(
        "--output-excel",
        type=str,
        default="cases/csv/radial_distribution.xlsx",
        help="Path to the output Excel file for radial distribution results."
    )
    parser.add_argument(
        "--residues",
        type=str,
        default="{'DPPC': ['NC3'], 'D3PC': ['NC3'], 'CHOL': ['ROH']}",
        help="A dictionary string defining residue groups for analysis. E.g., \"{'DPPC': ['NC3'], 'CHOL': ['ROH']}\""
    )
    parser.add_argument(
        "--n-circle",
        type=int,
        default=50,
        help="Number of concentric circles for radial analysis."
    )
    parser.add_argument(
        "--start-frame",
        type=int,
        default=0,
        help="Starting frame for analysis (0-indexed)."
    )
    parser.add_argument(
        "--stop-frame",
        type=int,
        help="Stopping frame for analysis (exclusive). Defaults to end of trajectory."
    )
    parser.add_argument(
        "--step-frame",
        type=int,
        default=1,
        help="Step size for frames during analysis."
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Parse residues_group from string
    try:
        residues_group_parsed = ast.literal_eval(args.residues)
        if not isinstance(residues_group_parsed, dict):
            raise ValueError("Residues argument must be a dictionary string.")
    except (ValueError, SyntaxError) as e:
        print(f"Error: Could not parse residues argument: {e}")
        print("Please ensure it's a valid dictionary string, e.g., \"{'DPPC': ['NC3'], 'CHOL': ['ROH']}\"")
        sys.exit(1)

    print("\n--- Initializing MDAnalysis Universe ---")
    try:
        import MDAnalysis as mda
        u = mda.Universe(args.gro_file)
    except Exception as e:
        print(f"Error loading MDAnalysis Universe: {e}")
        print("Please check if GRO file exists and is valid.")
        sys.exit(1)

    print("\n--- Running Radial Distribution Analysis ---")
    rad_analysis = CalRad(
        u,
        residues_group_parsed,
        ncircle=args.n_circle,
        filePath=args.output_excel
    )
    rad_analysis.run(
        start=args.start_frame,
        stop=args.stop_frame,
        step=args.step_frame
    )

    print("\n--- Analysis Finished ---")
    print(f"Results saved to: {args.output_excel}")

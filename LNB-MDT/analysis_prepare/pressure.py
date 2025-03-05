import numpy as np
from MDAnalysis import Universe
import matplotlib.pyplot as plt

u = Universe("E:/ach.gro", "E:/ach.xtc")  # 替换为你的文件路径
lipids = u.select_atoms("name P04")  # 假设 P 是脂质头基团原子

# 2. 设置参数
box = u.dimensions[:3]  # 模拟盒子尺寸 [x, y, z]
cube_size = 0.4  # 立方体边长 (nm)
N_segments = 100  # IK 轮廓分段数
n_cubes = np.floor(box / cube_size).astype(int)  # 每个维度立方体数
V_cube = cube_size ** 3  # 立方体体积

# 初始化压力张量数组 (n_frames, n_cubes_x, n_cubes_y, n_cubes_z, 3, 3)
n_frames = len(u.trajectory)
P_V = np.zeros((len(range(100, 200, 5)), n_cubes[0], n_cubes[1], n_cubes[2], 3, 3))

# 3. 计算每个时间步的压力张量
for t in range(100, 200, 5):
    ts = u.trajectory[t]
    frame_idx = ts.frame

    # 获取位置、速度、质量
    positions = u.atoms.positions / 10  # 转换为 nm
    velocities = u.atoms.velocities / 10  # 假设单位正确
    masses = u.atoms.masses

    # (1) 动能部分
    for i, (pos, vel, mass) in enumerate(zip(positions, velocities, masses)):
        cube_idx = np.floor(pos / cube_size).astype(int)
        if np.all(cube_idx < n_cubes) and np.all(cube_idx >= 0):
            for alpha in range(3):
                for beta in range(3):
                    P_V[frame_idx, cube_idx[0], cube_idx[1], cube_idx[2], alpha, beta] += (
                            mass * vel[alpha] * vel[beta] / V_cube
                    )

    # (2) 构型部分（示例：双体相互作用）
    pairs = [(i, j) for i in range(len(u.atoms)) for j in range(i + 1, len(u.atoms))]  # 所有粒子对
    for i, j in pairs:
        r_i = positions[i]
        r_j = positions[j]
        r_ij = r_j - r_i
        # 计算力 (假设简单 LJ 势，需替换为实际力场计算)
        dist = np.linalg.norm(r_ij)
        if dist < 1.0:  # 假设截断距离 1 nm
            F_ij = 24 * (2 * dist ** (-13) - dist ** (-7)) * r_ij / dist  # LJ 力示例
            for lambda_ in range(N_segments + 1):
                s = lambda_ / N_segments
                r_lambda = r_i + s * r_ij
                cube_idx = np.floor(r_lambda / cube_size).astype(int)
                if np.all(cube_idx < n_cubes) and np.all(cube_idx >= 0):
                    for alpha in range(3):
                        for beta in range(3):
                            P_V[frame_idx, cube_idx[0], cube_idx[1], cube_idx[2], alpha, beta] -= (
                                    F_ij[alpha] * r_ij[beta] / (N_segments * V_cube)
                            )

# 4. 时间平均
P_V_avg = np.mean(P_V, axis=0)

# 5. 球对称平均
center = np.mean(lipids.positions / 10, axis=0)  # 囊泡中心
r_bins = np.arange(0, 10, 0.1)  # 半径 bins (0-10 nm)
P_rr = np.zeros(len(r_bins) - 1)
P_T = np.zeros(len(r_bins) - 1)
counts = np.zeros(len(r_bins) - 1)

for i in range(n_cubes[0]):
    for j in range(n_cubes[1]):
        for k in range(n_cubes[2]):
            cube_pos = np.array([i, j, k]) * cube_size + cube_size / 2
            r_vec = cube_pos - center
            r = np.linalg.norm(r_vec)
            bin_idx = np.searchsorted(r_bins, r) - 1
            if 0 <= bin_idx < len(r_bins) - 1:
                # 简化为径向和切向分量（需实现完整变换矩阵 T）
                P = P_V_avg[i, j, k]
                P_rr[bin_idx] += P[2, 2]  # 假设 z 为径向，简化处理
                P_T[bin_idx] += (P[0, 0] + P[1, 1]) / 2
                counts[bin_idx] += 1

P_rr /= counts
P_T /= counts

# 6. 绘图
plt.plot(r_bins[:-1], P_rr, label="P_rr(r)")
plt.plot(r_bins[:-1], P_T, label="P_T(r)")
plt.xlabel("Radius (nm)")
plt.ylabel("Pressure (bar)")  # 单位需根据力场调整
plt.legend()
plt.show()
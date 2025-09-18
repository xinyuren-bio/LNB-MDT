import MDAnalysis as mda
import numpy as np

# 1. 加载你的模拟体系
# 你需要同时提供拓扑文件和轨迹文件
try:
    # 替换成你自己的文件名
    u = mda.Universe("/Users/renxinyu/LNB/LNB-MDT_v1.0/cases/lnb.gro", "/Users/renxinyu/LNB/LNB-MDT_v1.0/cases/md.xtc") 
except Exception as e:
    print(f"无法加载文件: {e}")
    print("请确保你的拓扑文件 (如 .gro, .pdb) 和轨迹文件 (.xtc) 路径正确。")
    # 为了演示，我们加载一个测试文件 (如果你的文件加载失败)
    print("\n--- 下面是使用 MDAnalysis 测试文件的演示 ---")
    from MDAnalysis.tests.datafiles import GRO, XTC
    u = mda.Universe(GRO, XTC)


# 2. 检查轨迹中有多少帧
print(f"轨迹总帧数: {len(u.trajectory)}")


# 3. 提取两帧之间的时间间隔 (你最关心的问题)

if len(u.trajectory) > 1:
    # 方法一：直接访问第二帧的 .dt 属性
    # (第0帧的 dt 可能是 0 或未定义，所以看第1帧最稳妥)
    u.trajectory[1]  # 跳转到第二帧 (索引为 1)
    time_per_frame = u.trajectory.ts.dt
    print(f"\n[方法一] 每一帧代表的时长 (ts.dt): {time_per_frame:.2f} ps")

    # 方法二：手动计算两帧的 .time 差值
    u.trajectory[0] # 跳转到第一帧 (索引为 0)
    time_0 = u.trajectory.ts.time
    
    u.trajectory[1] # 跳转到第二帧 (索引为 1)
    time_1 = u.trajectory.ts.time
    
    manual_diff = time_1 - time_0
    print(f"[方法二] 手动计算时间差: ({time_1:.2f} ps - {time_0:.2f} ps) = {manual_diff:.2f} ps")

else:
    print("轨迹中只有一帧，无法计算时间间隔。")


# 4. (可选) 遍历所有帧并打印它们的时间信息
print("\n--- 遍历并打印所有帧的时间 ---")
for ts in u.trajectory:
    # ts 是代表当前帧的 "timestep" 对象
    print(f"帧: {ts.frame}, 绝对时间 (ts.time): {ts.time:.2f} ps, 帧间隔 (ts.dt): {ts.dt:.2f} ps")


headatoms = u.select_atoms('name PO4 ROH')

gasatoms = u.select_atoms('name N2')

center_of_mass = headatoms.center_of_mass()

print(type(center_of_mass))
print(type(center_of_mass[0]))
print(center_of_mass[0])

# 检查质心是否有效
if center_of_mass is not None and not np.any(np.isnan(center_of_mass)):
    print(f"质心坐标: {center_of_mass}")
    # 使用逻辑运算符选择环形区域内的原子
    try:
        gas_radius = gasatoms.select_atoms(f'sphlayer 0.0 10.0 point {center_of_mass[0]:.3f} {center_of_mass[1]:.3f} {center_of_mass[2]:.3f}')
        print(f"选择的原子数量: {gas_radius.n_atoms}")
    except Exception as e:
        print(f"逻辑运算符语法错误: {e}")
        # 备用方案：使用简单的point选择
        try:
            gas_radius = gasatoms.select_atoms(f'point {center_of_mass[0]:.3f} {center_of_mass[1]:.3f} {center_of_mass[2]:.3f} 10.0')
            print(f"使用简单point选择语句，选择的原子数量: {gas_radius.n_atoms}")
        except Exception as e2:
            print(f"point语法也错误: {e2}")
else:
    print("质心无效，无法进行选择")

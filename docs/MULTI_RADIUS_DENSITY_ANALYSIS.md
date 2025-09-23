# 多半径密度分析可视化功能

## 概述

这个新功能允许您分析气泡中不同半径环形区域内气体密度随时间的变化，并提供多种可视化方式来展示三维数据（时间、半径范围、密度）。

## 计算逻辑

- **输入参数**: 最大半径 (max_radius) 和分段数量 (number_segments)
- **半径分段**: 将0到最大半径等分为指定数量的段
- **密度计算**: 计算每个环形区域内的密度，而不是累积密度
- **示例**: max_radius=50, number_segments=5 → 半径分段: [10, 20, 30, 40, 50]
  - 第1段: 0-10 Å 环形区域
  - 第2段: 10-20 Å 环形区域  
  - 第3段: 20-30 Å 环形区域
  - 第4段: 30-40 Å 环形区域
  - 第5段: 40-50 Å 环形区域

## 功能特点

- **多半径分析**: 同时分析多个半径的密度变化
- **多种可视化**: 支持折线图、热力图、3D表面图
- **并行处理**: 支持多核并行计算，提高分析效率
- **灵活配置**: 通过命令行参数或Python API灵活配置

## 可视化方案

### 1. 折线图 (推荐)
- **优点**: 最直观，容易比较不同半径的密度变化趋势
- **适用场景**: 展示密度随时间的变化模式，比较不同半径的差异
- **图表特点**: X轴为时间，Y轴为密度，不同半径用不同颜色线条表示

### 2. 热力图
- **优点**: 直观展示密度分布模式，颜色深浅表示密度大小
- **适用场景**: 识别密度分布的空间-时间模式
- **图表特点**: X轴为时间，Y轴为半径，颜色表示密度

### 3. 3D表面图
- **优点**: 立体展示三维关系，视觉效果丰富
- **适用场景**: 展示密度在时间和半径维度上的连续变化
- **图表特点**: X轴为时间，Y轴为半径，Z轴为密度

## 使用方法

### 命令行使用

```bash
# 基本使用 - 生成所有类型的图表
python analysis/density_multi_radius.py \
    --gro-file cases/lnb.gro \
    --xtc-file cases/md.xtc \
    --max-radius 50 \
    --number-segments 5 \
    --output-csv cases/csv/density_multi_radius_results.csv \
    --plot-dir cases/plots

# 只生成折线图
python analysis/density_multi_radius.py \
    --plot-type line \
    --max-radius 60 \
    --number-segments 6

# 使用并行处理
python analysis/density_multi_radius.py \
    --parallel \
    --n-jobs 4 \
    --max-radius 80 \
    --number-segments 8
```

### Python API使用

```python
from analysis.density import DensityRadius, DensityVisualizer
import MDAnalysis as mda

# 1. 加载轨迹
u = mda.Universe('cases/lnb.gro', 'cases/md.xtc')

# 2. 定义分析参数
residues_group = {'DPPC': ['PO4'], 'DUPC': ['PO4'], 'CHOL': ['ROH']}
gas_group = {'N2': ['N2']}
max_radius = 50
number_segments = 5

# 3. 运行分析
density_analysis = DensityRadius(
    u,
    ResiudeGroup=residues_group,
    GasGroup=gas_group,
    max_radius=max_radius,
    number_segments=number_segments,
    MW=14,
    filePath='results.csv',
    parallel=True,
    n_jobs=4
)

density_analysis.run()

# 4. 生成可视化
visualizer = DensityVisualizer(data_file='results.csv')

# 生成所有图表
visualizer.plot_all(save_dir='plots/')

# 或者单独生成特定类型的图表
visualizer.plot_line_chart(save_path='plots/line_chart.png')
visualizer.plot_heatmap(save_path='plots/heatmap.png')
visualizer.plot_3d_surface(save_path='plots/3d_surface.png')
```

## 参数说明

### 分析参数
- `--max-radius`: 最大半径 (Å)，默认50
- `--number-segments`: 分段数量，默认5
- `--MW`: 分子量，默认14
- `--residues`: 残基组定义
- `--gas-group`: 气体组定义

### 可视化参数
- `--plot-type`: 图表类型 (`line`, `heatmap`, `3d`, `all`)
- `--plot-dir`: 图表保存目录

### 性能参数
- `--parallel`: 启用并行处理
- `--n-jobs`: 并行任务数，-1表示使用所有CPU核心

## 输出文件

### 数据文件
- **CSV文件**: 包含frame、radius、density三列的数据文件
- **格式**: 每行代表一个时间帧-半径组合的密度值

### 图表文件
- **折线图**: `density_line_chart.png`
- **热力图**: `density_heatmap.png`
- **3D表面图**: `density_3d_surface.png`

## 示例数据

运行示例程序查看可视化效果：

```bash
python examples/density_visualization_example.py
```

## 性能优化建议

1. **并行处理**: 对于大轨迹文件，使用 `--parallel` 参数
2. **帧采样**: 使用 `--step-frame` 参数减少分析帧数
3. **半径选择**: 根据实际需求选择合适的半径范围
4. **内存管理**: 对于超大轨迹，考虑分段分析

## 常见问题

### Q: 如何选择合适的最大半径和分段数量？
A: 建议根据气泡大小选择最大半径，分段数量5-10个比较合适。例如：气泡直径100Å，可选择max_radius=50, number_segments=5，得到10Å间隔的环形区域。

### Q: 密度单位是什么？
A: 输出密度单位为 kg/m³。

### Q: 如何解释热力图中的颜色？
A: 颜色越深（接近紫色）表示密度越高，颜色越浅（接近黄色）表示密度越低。

### Q: 3D图在什么情况下最有用？
A: 当您需要同时观察密度在时间和半径两个维度上的连续变化时，3D图最有价值。

## 扩展功能

这个模块设计为可扩展的，您可以：

1. **添加新的可视化类型**: 在 `DensityVisualizer` 类中添加新方法
2. **自定义图表样式**: 修改matplotlib参数
3. **添加统计分析**: 计算密度统计指标
4. **导出其他格式**: 支持导出为Excel、JSON等格式

## 依赖包

确保安装以下Python包：

```bash
pip install matplotlib seaborn pandas numpy MDAnalysis joblib scipy
```

# 分析模块命令行使用指南

本指南说明如何使用命令行参数运行各个分析模块。

## 通用参数

所有分析模块都支持以下通用参数：

- `--gro-file`: GRO文件路径（拓扑文件）
- `--xtc-file`: XTC文件路径（轨迹文件）
- `--output-csv`: 输出CSV文件路径
- `--residues`: 残基组字典字符串
- `--parallel`: 启用并行处理
- `--n-jobs`: 并行作业数量（-1表示使用所有可用CPU核心）
- `--start-frame`: 分析起始帧（0索引）
- `--stop-frame`: 分析结束帧（独占）
- `--step-frame`: 帧步长
- `--verbose`: 启用详细输出

## 各模块特定参数

### 1. PCA分析 (pca.py)

```bash
python analysis/pca.py --gro-file cases/lnb.gro --xtc-file cases/md.xtc --output-csv cases/csv/pca_results.csv --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" --parallel --verbose
```

### 2. 面积分析 (area.py)

```bash
python analysis/area.py --gro-file cases/lnb.gro --xtc-file cases/md.xtc --output-csv cases/csv/area_results.csv --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" --k-value 20 --max-normal-angle 140 --parallel --verbose
```

特定参数：
- `--k-value`: Voronoi镶嵌的K值
- `--max-normal-angle`: 最大法线角度（度）

### 3. 曲率分析 (curvature.py)

```bash
python analysis/curvature.py --gro-file cases/lnb.gro --xtc-file cases/md.xtc --output-csv cases/csv/curvature_results.csv --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" --k-value 20 --method mean --parallel --verbose
```

特定参数：
- `--k-value`: 曲率计算的K值
- `--method`: 曲率类型（'mean' 或 'gaussian'）

### 4. 高度分析 (height.py)

```bash
python analysis/height.py --gro-file cases/lnb.gro --xtc-file cases/md.xtc --output-csv cases/csv/height_results.csv --residues "{'DPPC': (['PO4'], ['C4B', 'C4A']), 'CHOL':(['ROH'], ['R5'])}" --k-value 20 --parallel --verbose
```

特定参数：
- `--k-value`: 高度计算的K值

### 5. 聚类分析 (cluster.py)

```bash
python analysis/cluster.py --gro-file cases/lnb.gro --xtc-file cases/md.xtc --output-csv cases/csv/cluster_results.csv --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" --cutoff 8.0 --parallel --verbose
```

特定参数：
- `--cutoff`: 聚类截止距离（埃）

### 6. 各向异性分析 (anisotropy.py)

```bash
python analysis/anisotropy.py --gro-file cases/lnb.gro --xtc-file cases/md.xtc --output-csv cases/csv/anisotropy_results.csv --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" --parallel --verbose
```

### 7. 回转半径分析 (gyration.py)

```bash
python analysis/gyration.py --gro-file cases/lnb.gro --xtc-file cases/md.xtc --output-csv cases/csv/gyration_results.csv --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" --parallel --verbose
```

### 8. Sz序参数分析 (sz.py)

```bash
python analysis/sz.py --gro-file cases/lnb.gro --xtc-file cases/md.xtc --output-csv cases/csv/sz_results.csv --residues "{'DPPC': ['PO4'], 'DUPC': ['PO4']}" --chain sn1 --k-value 15 --parallel --verbose
```

特定参数：
- `--chain`: 链类型（'sn1', 'sn2', 或 'both'）
- `--k-value`: Sz计算的K值

### 9. N-聚类分析 (n_cluster.py)

```bash
python analysis/n_cluster.py --gro-file cases/lnb.gro --xtc-file cases/md.xtc --output-csv cases/csv/ncluster_results.csv --residues "{'DAPC': ['GL1', 'GL2'], 'DPPC': ['PO4']}" --cutoff 12.0 --n-cutoff 10 --parallel --verbose
```

特定参数：
- `--cutoff`: 聚类截止距离（埃）
- `--n-cutoff`: 最小聚类大小阈值

### 10. 径向分布分析 (rad.py)

```bash
python analysis/rad.py --gro-file cases/lnb.gro --output-excel cases/csv/radial_distribution.xlsx --residues "{'DPPC': ['NC3'], 'CHOL': ['ROH']}" --n-circle 50
```

特定参数：
- `--output-excel`: 输出Excel文件路径
- `--n-circle`: 径向分析的同心圆数量

## 使用示例

### 基本使用（串行处理）

```bash
python analysis/pca.py --gro-file my_system.gro --xtc-file my_trajectory.xtc --residues "{'DPPC': ['PO4']}" --output-csv results.csv
```

### 并行处理

```bash
python analysis/pca.py --gro-file my_system.gro --xtc-file my_trajectory.xtc --residues "{'DPPC': ['PO4']}" --output-csv results.csv --parallel --n-jobs 4 --verbose
```

### 指定帧范围

```bash
python analysis/pca.py --gro-file my_system.gro --xtc-file my_trajectory.xtc --residues "{'DPPC': ['PO4']}" --output-csv results.csv --start-frame 100 --stop-frame 1000 --step-frame 10
```

## 注意事项

1. **residues参数格式**: 必须使用有效的Python字典字符串格式
2. **文件路径**: 确保GRO和XTC文件存在且可访问
3. **并行处理**: 使用`--parallel`标志启用并行处理，使用`--n-jobs`指定作业数量
4. **输出目录**: 确保输出目录存在，否则程序会尝试创建

## 错误处理

如果遇到参数解析错误，程序会显示详细的错误信息和正确的使用格式。常见的错误包括：

- 无效的residues字典格式
- 文件不存在
- 参数类型错误

## 与UI的兼容性

这些命令行参数功能不会影响原有的UI调用。UI仍然可以通过直接调用分析类的构造函数来使用这些模块。

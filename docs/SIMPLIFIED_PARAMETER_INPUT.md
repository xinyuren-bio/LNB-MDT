# 简化参数输入指南

## 概述

为了解决输入复杂字典字符串的麻烦，我们为所有analysis模块添加了简化的参数输入格式支持。现在您可以使用更简单、更直观的方式来指定residues和gas-group参数。

## 支持的输入格式

### 1. 简单格式 (推荐)

**Residues参数:**
```bash
# 基本格式: RESIDUE:ATOM
-r DPPC:PO4,CHOL:ROH

# 多原子格式: RESIDUE:ATOM1+ATOM2
-r DPPC:PO4+GLY,CHOL:ROH

# 多个残基
-r DPPC:PO4,DUPC:PO4,CHOL:ROH
```

**Gas-group参数:**
```bash
# 基本格式
-a N2:N2

# 多个气体
-a N2:N2,O2:O2

# 只有气体名（原子名与气体名相同）
-a N2
```

### 2. 配置文件格式

**JSON配置文件:**
```bash
# 使用JSON配置文件
-r @cases/config/residues_config.json
-a @cases/config/gas_config.json
```

**文本配置文件:**
```bash
# 使用文本配置文件
-r @cases/config/residues.txt
-a @cases/config/gas.txt
```

### 3. 传统字典格式 (仍然支持)

```bash
# 字典字符串格式
-r "{'DPPC': ['PO4'], 'CHOL': ['ROH']}"
-a "{'N2': ['N2']}"
```

## 使用示例

### 命令行使用

**原来的复杂方式:**
```bash
python analysis/densitywithframe.py \
    --gro-file cases/lnb.gro \
    --xtc-file cases/md.xtc \
    --residues "{'DPPC': ['PO4'], 'DUPC': ['PO4'], 'CHOL': ['ROH']}" \
    --gas-group "{'N2': ['N2']}" \
    --output-csv results.csv
```

**现在的简单方式:**
```bash
# 使用短参数 + 简单格式
python analysis/densitywithframe.py \
    -g cases/lnb.gro \
    -x cases/md.xtc \
    -r DPPC:PO4,DUPC:PO4,CHOL:ROH \
    -a N2:N2 \
    -o results.csv
```

**使用配置文件:**
```bash
python analysis/densitywithframe.py \
    -g cases/lnb.gro \
    -x cases/md.xtc \
    -r @cases/config/residues_config.json \
    -a @cases/config/gas_config.json \
    -o results.csv
```

### Python API使用

```python
from analysis.parameter_utils import parse_residues_simple, parse_gas_group_simple

# 简单格式
residues = parse_residues_simple('DPPC:PO4,CHOL:ROH')
gas_group = parse_gas_group_simple('N2:N2')

# 多原子格式
residues = parse_residues_simple('DPPC:PO4+GLY,CHOL:ROH')

# 配置文件
residues = parse_residues_simple('@config/residues.json')
gas_group = parse_gas_group_simple('@config/gas.json')
```

## 配置文件格式

### JSON配置文件

**residues_config.json:**
```json
{
  "DPPC": ["PO4"],
  "DUPC": ["PO4"],
  "CHOL": ["ROH"]
}
```

**gas_config.json:**
```json
{
  "N2": ["N2"],
  "O2": ["O2"]
}
```

### 文本配置文件

**residues.txt:**
```
# Residues configuration file
# Format: RESIDUE:ATOM1+ATOM2
DPPC:PO4
DUPC:PO4
CHOL:ROH
```

**gas.txt:**
```
# Gas configuration file
# Format: GAS:ATOM
N2:N2
O2:O2
```

## 参数映射表

### 短参数别名

| 长参数 | 短参数 | 说明 |
|---------|---------|------|
| `--gro-file` | `-g` | GRO文件路径 |
| `--xtc-file` | `-x` | XTC文件路径 |
| `--output-csv` | `-o` | 输出CSV文件路径 |
| `--residues` | `-r` | 残基组定义 |
| `--gas-group` | `-a` | 气体组定义 |
| `--MW` | `-m` | 分子量 |
| `--radius` | `-R` | 半径 |
| `--parallel` | `-p` | 启用并行处理 |
| `--n-jobs` | `-j` | 并行任务数 |
| `--start-frame` | `-s` | 起始帧 |
| `--stop-frame` | `-e` | 结束帧 |
| `--step-frame` | `-t` | 帧步长 |
| `--verbose` | `-v` | 详细输出 |

## 实际使用对比

### 密度分析示例

**原来的方式:**
```bash
python analysis/densitywithframe.py \
    --gro-file cases/lnb.gro \
    --xtc-file cases/md.xtc \
    --residues "{'DPPC': ['PO4'], 'DUPC': ['PO4'], 'CHOL': ['ROH']}" \
    --gas-group "{'N2': ['N2']}" \
    --MW 14 \
    --radius 50 \
    --output-csv cases/csv/density_results.csv \
    --parallel \
    --n-jobs 4
```

**现在的方式:**
```bash
python analysis/densitywithframe.py \
    -g cases/lnb.gro \
    -x cases/md.xtc \
    -r DPPC:PO4,DUPC:PO4,CHOL:ROH \
    -a N2:N2 \
    -m 14 \
    -R 50 \
    -o cases/csv/density_results.csv \
    -p \
    -j 4
```

**使用配置文件:**
```bash
python analysis/densitywithframe.py \
    -g cases/lnb.gro \
    -x cases/md.xtc \
    -r @cases/config/residues_config.json \
    -a @cases/config/gas_config.json \
    -m 14 \
    -R 50 \
    -o cases/csv/density_results.csv \
    -p \
    -j 4
```

## 优势

1. **输入更简单**: 不需要输入复杂的引号和括号
2. **更直观**: 格式更接近自然语言
3. **支持配置文件**: 可以保存常用配置，避免重复输入
4. **向后兼容**: 仍然支持原有的字典格式
5. **短参数**: 所有参数都有简短的别名

## 注意事项

1. **空格处理**: 参数中的空格会被自动处理
2. **大小写敏感**: 残基名和原子名区分大小写
3. **配置文件路径**: 使用`@`前缀指定配置文件路径
4. **错误处理**: 如果格式不正确，会显示详细的错误信息和格式说明

## 创建配置文件

您可以使用以下命令创建示例配置文件：

```bash
python analysis/parameter_utils.py
```

这会在 `cases/config/` 目录下创建示例配置文件。

## 支持的模块

所有analysis模块都支持简化的参数输入：

- `densitywithframe.py`
- `densitywithradius.py`
- `area.py`
- `height.py`
- `curvature.py`
- `pca.py`
- `cluster.py`
- `anisotropy.py`
- `gyration.py`
- `sz.py`
- `n_cluster.py`
- `rad.py`

现在您可以享受更简单、更直观的命令行体验了！

分析模块详解
============

LNB-MDT提供了丰富的分子动力学分析模块，每个模块都针对特定的物理性质进行分析。

模块概览
--------

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">📐 PCA分析</h3>
   <p>主成分分析，研究分子构象变化</p>
   <ul style="margin-bottom: 0;">
   <li>降维分析</li>
   <li>构象变化</li>
   <li>运动模式</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">📏 面积分析</h3>
   <p>APL面积计算</p>
   <ul style="margin-bottom: 0;">
   <li>分子面积</li>
   <li>密度分布</li>
   <li>包装效率</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">🌊 曲率分析</h3>
   <p>膜曲率计算（平均/高斯）</p>
   <ul style="margin-bottom: 0;">
   <li>平均曲率</li>
   <li>高斯曲率</li>
   <li>膜形变</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">📊 高度分析</h3>
   <p>分子高度分布分析</p>
   <ul style="margin-bottom: 0;">
   <li>Z坐标分布</li>
   <li>膜厚度</li>
   <li>表面粗糙度</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #fa709a 0%, #fee140 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">🔗 聚类分析</h3>
   <p>分子聚集行为分析</p>
   <ul style="margin-bottom: 0;">
   <li>聚集模式</li>
   <li>聚类大小</li>
   <li>相互作用</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #a8edea 0%, #fed6e3 100%); color: #333; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">🎯 各向异性分析</h3>
   <p>分子取向各向异性计算</p>
   <ul style="margin-bottom: 0;">
   <li>取向分布</li>
   <li>有序参数</li>
   <li>分子排列</li>
   </ul>
   </div>

   </div>

详细模块说明
------------

PCA分析 (pca.py)
~~~~~~~~~~~~~~~~

**功能描述**
主成分分析用于研究脂质分子的构象变化和运动模式。

**算法原理**
- 对分子坐标进行主成分分析
- 提取主要的运动模式
- 降维到主要成分空间

**关键参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- **residues**: 残基组字典，指定要分析的分子类型
- **n_components**: 主成分数量（默认：3）
- **start_frame**: 起始帧（默认：0）
- **stop_frame**: 结束帧（默认：-1，表示到最后）
- **step_frame**: 帧步长（默认：1）

   </div>

**使用示例**

.. code:: bash

   python analysis/pca.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/pca_results.csv \
     --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" \
     --n-components 3 \
     --parallel \
     --verbose

**输出结果**
- CSV文件包含每个帧的主成分值
- 可用于分析分子构象变化趋势
- 支持可视化分析结果

面积分析 (area.py)
~~~~~~~~~~~~~~~~~~

**功能描述**
使用Voronoi镶嵌方法计算脂质分子的面积分布。

**算法原理**
- 构建Voronoi图
- 计算每个分子的Voronoi面积
- 分析面积分布和变化

**关键参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- **k-value**: Voronoi镶嵌的k值（默认：20）
- **max-normal-angle**: 最大法线角度（默认：140度）
- **residues**: 残基组字典

   </div>

**使用示例**

.. code:: bash

   python analysis/area.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/area_results.csv \
     --residues "{'DPPC': ['PO4']}" \
     --k-value 20 \
     --max-normal-angle 140 \
     --parallel \
     --verbose

**输出结果**
- 每个分子的Voronoi面积
- 面积分布统计信息
- 可用于分析膜密度和包装

曲率分析 (curvature.py)
~~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
计算脂质膜的平均曲率和高斯曲率。

**算法原理**
- 基于局部表面拟合
- 计算曲率张量
- 提取平均曲率和高斯曲率

**关键参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- **method**: 曲率类型（'mean' 或 'gaussian'）
- **k-value**: 曲率计算的k值（默认：20）
- **residues**: 残基组字典

   </div>

**使用示例**

.. code:: bash

   python analysis/curvature.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/curvature_results.csv \
     --residues "{'DPPC': ['PO4']}" \
     --k-value 20 \
     --method mean \
     --parallel \
     --verbose

**输出结果**
- 每个分子的曲率值
- 曲率分布统计
- 可用于分析膜形变和稳定性

高度分析 (height.py)
~~~~~~~~~~~~~~~~~~~~~

**功能描述**
分析脂质分子的高度分布和膜厚度。

**算法原理**
- 计算分子在Z方向的位置
- 分析高度分布
- 计算膜厚度和表面粗糙度

**关键参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- **k-value**: 高度计算的k值（默认：20）
- **residues**: 残基组字典，支持多组原子

   </div>

**使用示例**

.. code:: bash

   python analysis/height.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/height_results.csv \
     --residues "{'DPPC': (['PO4'], ['C4B', 'C4A']), 'CHOL': (['ROH'], ['R5'])}" \
     --k-value 20 \
     --parallel \
     --verbose

**输出结果**
- 每个分子的高度值
- 高度分布统计
- 膜厚度分析

聚类分析 (cluster.py)
~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
分析脂质分子的聚集行为和聚类模式。

**算法原理**
- 基于距离的聚类算法
- 识别分子聚集
- 分析聚类大小和分布

**关键参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- **cutoff**: 聚类截止距离（默认：8.0埃）
- **residues**: 残基组字典

   </div>

**使用示例**

.. code:: bash

   python analysis/cluster.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/cluster_results.csv \
     --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" \
     --cutoff 8.0 \
     --parallel \
     --verbose

**输出结果**
- 聚类大小分布
- 聚类数量统计
- 聚集行为分析

各向异性分析 (anisotropy.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
计算分子取向的各向异性参数。

**算法原理**
- 计算分子取向向量
- 分析取向分布
- 计算各向异性参数

**关键参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- **residues**: 残基组字典

   </div>

**使用示例**

.. code:: bash

   python analysis/anisotropy.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/anisotropy_results.csv \
     --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" \
     --parallel \
     --verbose

**输出结果**
- 各向异性参数
- 取向分布统计
- 分子排列分析

回转半径分析 (gyration.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
计算分子的回转半径，反映分子的紧凑程度。

**算法原理**
- 计算分子质心
- 计算回转半径
- 分析分子形状变化

**关键参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- **residues**: 残基组字典

   </div>

**使用示例**

.. code:: bash

   python analysis/gyration.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/gyration_results.csv \
     --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" \
     --parallel \
     --verbose

**输出结果**
- 回转半径值
- 形状变化分析
- 分子紧凑度

Sz序参数分析 (sz.py)
~~~~~~~~~~~~~~~~~~~~~

**功能描述**
计算脂质链的Sz序参数，反映链的有序程度。

**算法原理**
- 计算链取向向量
- 计算Sz序参数
- 分析链有序性

**关键参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- **chain**: 链类型（'sn1', 'sn2', 或 'both'）
- **k-value**: Sz计算的k值（默认：15）
- **residues**: 残基组字典

   </div>

**使用示例**

.. code:: bash

   python analysis/sz.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/sz_results.csv \
     --residues "{'DPPC': ['PO4'], 'DUPC': ['PO4']}" \
     --chain sn1 \
     --k-value 15 \
     --parallel \
     --verbose

**输出结果**
- Sz序参数值
- 链有序性分析
- 相变行为

N-聚类分析 (n_cluster.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
统计聚类数量，分析聚集模式。

**算法原理**
- 基于距离的聚类
- 统计聚类数量
- 分析聚集模式

**关键参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- **cutoff**: 聚类截止距离（默认：12.0埃）
- **n-cutoff**: 最小聚类大小阈值（默认：10）
- **residues**: 残基组字典

   </div>

**使用示例**

.. code:: bash

   python analysis/n_cluster.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/ncluster_results.csv \
     --residues "{'DAPC': ['GL1', 'GL2'], 'DPPC': ['PO4']}" \
     --cutoff 12.0 \
     --n-cutoff 10 \
     --parallel \
     --verbose

**输出结果**
- 聚类数量统计
- 聚集模式分析
- 相互作用强度

径向分布分析 (rad.py)
~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
计算径向分布函数，分析分子间的距离分布。

**算法原理**
- 计算分子间距离
- 构建径向分布函数
- 分析相互作用

**关键参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- **n-circle**: 径向分析的同心圆数量（默认：50）
- **residues**: 残基组字典

   </div>

**使用示例**

.. code:: bash

   python analysis/rad.py \
     --gro-file cases/lnb.gro \
     --output-excel results/radial_distribution.xlsx \
     --residues "{'DPPC': ['NC3'], 'CHOL': ['ROH']}" \
     --n-circle 50

**输出结果**
- Excel文件包含径向分布数据
- 距离分布统计
- 相互作用分析

参数优化建议
------------

k值选择
~~~~~~~~

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 15px; border-radius: 8px; border-left: 4px solid #2196f3;">

**k值选择原则：**

- **小系统**: k = 10-15
- **中等系统**: k = 15-25  
- **大系统**: k = 25-35
- **高密度**: 增加k值
- **低密度**: 减少k值

**优化方法：**
使用机器学习模块的k值优化器：

.. code:: python

   from machine_learning import KValueOptimizer
   optimizer = KValueOptimizer('area')
   best_k = optimizer.optimize()

   </div>

截止距离选择
~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**截止距离选择：**

- **聚类分析**: 8-12埃
- **相互作用**: 5-8埃
- **长程相互作用**: 12-20埃

**选择依据：**
- 分子大小
- 相互作用强度
- 系统密度

   </div>

并行处理优化
~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**并行处理建议：**

- **CPU核心数**: 使用 `--n-jobs -1` 自动检测
- **内存考虑**: 大系统减少并行数
- **I/O限制**: SSD硬盘可增加并行数

**性能优化：**
- 使用SSD存储轨迹文件
- 增加系统内存
- 优化网络文件系统

   </div>

结果解读指南
------------

CSV文件格式
~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

**标准CSV格式：**

.. code:: csv

   # Created by LNB-MDT v1.0
   # PCA Analysis
   # TYPE:Bubble
   # Parameters:{'DPPC': ['PO4']}
   Frames,Values
   0,0.787
   1,0.801
   2,0.800
   ...

**列说明：**
- **Frames**: 帧编号
- **Values**: 分析结果值
- **注释行**: 包含分析参数和类型

   </div>

数据可视化
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**可视化建议：**

1. **时间序列图**: 显示分析值随时间变化
2. **分布直方图**: 显示数值分布
3. **相关性分析**: 不同分析间的相关性
4. **统计摘要**: 均值、标准差、分位数

**使用图表模块：**
- 导入CSV数据
- 选择图表类型
- 自定义样式
- 导出高质量图片

   </div>

统计分析
~~~~~~~~

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**统计指标：**

- **均值**: 反映平均水平
- **标准差**: 反映变异性
- **分位数**: 反映分布特征
- **趋势**: 线性回归分析
- **周期性**: 傅里叶变换

**Python分析示例：**

.. code:: python

   import pandas as pd
   import numpy as np
   
   # 读取结果
   data = pd.read_csv('results/pca_results.csv')
   
   # 基本统计
   stats = data['Values'].describe()
   
   # 趋势分析
   from scipy import stats
   slope, intercept, r_value, p_value, std_err = stats.linregress(data['Frames'], data['Values'])

   </div>

最佳实践
--------

分析流程建议
~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**推荐分析流程：**

1. **预处理**: 检查轨迹质量，去除异常帧
2. **参数优化**: 使用ML模块优化关键参数
3. **并行分析**: 启用并行处理提高效率
4. **结果验证**: 检查结果的合理性
5. **可视化**: 使用图表模块可视化结果
6. **统计分析**: 进行深入的统计分析

   </div>

质量控制
~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**质量控制检查：**

- **数据完整性**: 检查是否有缺失值
- **数值合理性**: 检查结果是否在合理范围内
- **收敛性**: 检查分析是否收敛
- **重现性**: 多次运行验证结果一致性

   </div>

性能优化
~~~~~~~~

.. raw:: html

   <div style="background-color: #fce4ec; padding: 15px; border-radius: 8px; border-left: 4px solid #e91e63;">

**性能优化技巧：**

- **分段处理**: 大轨迹文件分段分析
- **参数调优**: 根据系统特点调整参数
- **硬件优化**: 使用SSD和充足内存
- **算法选择**: 选择最适合的分析算法

   </div>

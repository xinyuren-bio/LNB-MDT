快速开始指南
============

本指南将帮助您在5分钟内开始使用LNB-MDT进行脂质纳米泡分析。

配置VMD路径
----------

首次使用LNB-MDT需要配置VMD路径。VMD用于分子可视化和轨迹分析。

1. **编辑配置文件**
   
   打开项目根目录的 `config.ini` 文件，修改 `vmd_path` 为您的VMD实际安装路径：

.. code:: text

   # 常见路径示例
   Windows: C:/Program Files/VMD/vmd.exe
   macOS:   /Applications/VMD.app/Contents/vmd/vmd_MACOSXARM64
   Linux:   /usr/local/bin/vmd

2. **保存并重启**
   
   保存配置文件后重新启动程序。

启动程序
--------

图形界面启动
~~~~~~~~~~~~

使用图形界面是开始使用LNB-MDT最简单的方式：

.. code:: python

   # 激活环境
   conda activate LNB-MDT
   
   # 启动主程序
   python main.py

启动后您将看到LNB-MDT的主界面，包含以下功能模块：

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); color: white; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-center: 0;">🧬 Generation Module</h4>
   </div>

   <div style="background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); color: white; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-center: 0;">📊 Analysis Module</h4>
   </div>

   <div style="background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%); color: white; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-center: 0;">📈 Figure Module</h4>
   </div>

   <div style="background: linear-gradient(135deg, #fa709a 0%, #fee140 100%); color: white; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-center: 0;">🔧 VMD Module</h4>
   </div>

   </div>

命令行启动
~~~~~~~~~~

对于批量处理和自动化分析，可以使用命令行工具：

.. code:: python

   # 激活环境
   conda activate LNB-MDT
   
   # 查看帮助信息
   python analysis/pca.py --help

基本分析流程
------------

准备数据文件
~~~~~~~~~~~~

LNB-MDT需要以下文件进行分析：

- **GRO文件**: 分子拓扑结构文件
- **XTC文件**: 分子动力学轨迹文件

项目包含示例数据文件：
- `cases/lnb.gro` - 示例拓扑文件  
- `cases/md.xtc` - 示例轨迹文件

选择分析类型
~~~~~~~~~~~~

LNB-MDT提供多种分析类型：

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #7b1fa2;">📐 Anisotropy</h4>
   <p style="margin-bottom: 0;">主成分分析，研究分子构象变化</p>
   </div>

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #388e3c;">📏 APL</h4>
   <p style="margin-bottom: 0;">Voronoi镶嵌面积计算</p>
   </div>

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #f57c00;">🌊 SZ</h4>
   <p style="margin-bottom: 0;">膜曲率计算（平均/高斯）</p>
   </div>

   <div style="background-color: #fce4ec; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #c2185b;">📊 Cluster</h4>
   <p style="margin-bottom: 0;">分子聚集行为分析</p>
   </div>

   </div>

配置参数
~~~~~~~~

关键参数包括：

- **残基组**: 指定要分析的分子类型和原子
- **帧范围**: 选择分析的时间范围  
- **计算参数**: k值、截止距离等
- **并行处理**: 启用多核加速

简化参数输入
^^^^^^^^^^^^

LNB-MDT现在支持更简单的参数输入方式，让命令行使用更加便捷：

**短参数别名:**
.. code:: text

   -g  --gro-file      GRO文件路径
   -x  --xtc-file      XTC文件路径  
   -o  --output-csv    输出CSV文件路径
   -r  --residues      残基组定义
   -a  --gas-group     气体组定义
   -m  --MW           分子量
   -R  --radius       半径
   -p  --parallel     启用并行处理
   -j  --n-jobs       并行任务数
   -s  --start-frame   起始帧
   -e  --stop-frame    结束帧
   -t  --step-frame    帧步长
   -v  --verbose       详细输出

**简化的residues和gas-group格式:**
.. code:: text

   # 简单格式（推荐）
   -r DPPC:PO4,CHOL:ROH
   -a N2:N2
   
   # 多原子格式
   -r DPPC:PO4+GLY,CHOL:ROH
   
   
   # 传统字典格式（仍然支持）
   -r "{'DPPC': ['PO4'], 'CHOL': ['ROH']}"

运行分析
~~~~~~~~

图形界面运行
^^^^^^^^^^^^

1. 在界面中加载GRO和XTC文件
2. 选择分析类型
3. 配置参数
4. 点击"运行"按钮
5. 查看结果

命令行运行
^^^^^^^^^^

LNB-MDT支持简化的命令行参数输入，让您更轻松地使用命令行工具：

**传统方式（仍然支持）:**
.. code-block:: python

   # PCA分析示例
   python analysis/pca.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/pca_results.csv \
     --residues "{'DPPC': ['PO4']}" \
     --parallel \
     --verbose

**新的简化方式（推荐）:**
.. code-block:: python

   # 使用短参数和简单格式
   python analysis/pca.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/pca_results.csv \
     -r DPPC:PO4 \
     -p \
     -v

查看结果
~~~~~~~~

分析完成后，LNB-MDT会生成以下输出：

- **CSV文件**: 包含分析结果的数值数据
- **图表**: 可视化分析结果  
- **日志**: 分析过程的详细信息

结果解读：

- 查看CSV文件中的数值结果
- 使用图表模块可视化数据
- 结合VMD进行分子可视化

实际示例
--------

PCA分析
~~~~~~~

分析脂质分子的构象变化：

**传统方式:**
.. code-block:: python

   python analysis/pca.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/pca_test.csv \
     --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" \
     --start-frame 0 \
     --stop-frame 100 \
     --parallel \
     --verbose

**简化方式:**
.. code-block:: python

   python analysis/pca.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/pca_test.csv \
     -r DPPC:PO4,CHOL:ROH \
     -s 0 \
     -e 100 \
     -p \
     -v

面积分析
~~~~~~~~

计算脂质分子的Voronoi镶嵌面积：

**传统方式:**
.. code-block:: python

   python analysis/area.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/area_test.csv \
     --residues "{'DPPC': ['PO4']}" \
     --k-value 20 \
     --max-normal-angle 140 \
     --parallel \
     --verbose

**简化方式:**
.. code-block:: python

   python analysis/area.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/area_test.csv \
     -r DPPC:PO4 \
     -k 20 \
     --max-normal-angle 140 \
     -p \
     -v

曲率分析
~~~~~~~~

计算脂质膜的曲率特性：

**传统方式:**
.. code-block:: python

   python analysis/curvature.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/curvature_test.csv \
     --residues "{'DPPC': ['PO4']}" \
     --k-value 20 \
     --method mean \
     --parallel \
     --verbose

**简化方式:**
.. code-block:: python

   python analysis/curvature.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/curvature_test.csv \
     -r DPPC:PO4 \
     -k 20 \
     -M mean \
     -p \
     -v

密度分析
~~~~~~~~

分析气泡中气体密度随时间的变化：

**简化方式（推荐）:**
.. code-block:: python

   python analysis/densitywithframe.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -o results/density_test.csv \
     -r DPPC:PO4,CHOL:ROH \
     -a N2:N2 \
     -m 14 \
     -R 50 \
     -p \
     -v

机器学习功能
------------

LNB-MDT集成了强大的机器学习功能，包括参数优化、异常检测和属性预测。

参数优化
~~~~~~~~

自动参数优化功能可以帮助找到最佳的分析参数：

.. code:: python

   from machine_learning import AnalysisParameterOptimizer
   
   # 创建优化器
   optimizer = AnalysisParameterOptimizer('area')
   
   # 运行优化
   results = optimizer.optimize()
   print(f"最佳参数: {results['best_parameters']}")

异常检测
~~~~~~~~

异常模式检测可以识别轨迹中的异常行为：

.. code:: python

   from machine_learning import MDAnomalyDetector
   
   # 创建检测器
   detector = MDAnomalyDetector(method='isolation_forest')
   
   # 分析轨迹
   results = detector.analyze_trajectory(
       gro_file="cases/lnb.gro",
       xtc_file="cases/md.xtc",
       residues={'DPPC': ['PO4']}
   )

属性预测
~~~~~~~~

分子属性预测可以基于轨迹数据预测分子的物理化学性质：

.. code:: python

   from machine_learning import MDPropertyPredictor
   
   # 创建预测器
   predictor = MDPropertyPredictor(
       model_type='random_forest',
       target_property='diffusion_coefficient'
   )
   
   # 训练模型
   results = predictor.fit(X_train, y_train)

VMD集成
--------

LNB-MDT支持与VMD的无缝集成，用于分子可视化和轨迹分析。

VMD路径配置
~~~~~~~~~~~

首次使用需要配置VMD路径：

1. **找到VMD安装路径**

.. code:: text

   Windows: 通常在 C:/Program Files/VMD/vmd.exe
   macOS:   通常在 /Applications/VMD.app/Contents/vmd/vmd_MACOSXARM64
   Linux:   通常在 /usr/local/bin/vmd

2. **编辑配置文件**
   
   打开项目根目录的 `config.ini` 文件，修改 `vmd_path` 为您的VMD实际安装路径：

.. code:: ini

   [VMD]
   vmd_path = /Applications/VMD.app/Contents/vmd/vmd_MACOSXARM64

3. **验证配置**
   
   保存配置文件后重新启动LNB-MDT程序。

启动VMD
~~~~~~~

图形界面启动：

1. 点击"Start VMD"按钮
2. 等待VMD启动
3. 拖拽CSV文件到VMD窗口
4. 选择分子进行可视化

命令行启动：

.. code:: python

   # 启动VMD
   python -c "from modules.vmd_control import VMDTcp; vmd = VMDTcp(); vmd.start()"

可视化操作
~~~~~~~~~~

操作步骤：

1. 在LNB-MDT中加载分析结果
2. 选择要可视化的帧和分子
3. VMD自动跳转到对应帧
4. 高亮显示选中的分子
5. 调整可视化参数

下一步
------

恭喜！您已经成功完成了LNB-MDT的快速开始！

接下来可以：

- 学习 :doc:`analysis_modules` 的深度使用  
- 探索 :doc:`machine_learning` 功能
- 查看 :doc:`api_reference` 了解API详情

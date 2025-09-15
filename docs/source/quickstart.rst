快速开始指南
============

本指南将帮助您在5分钟内开始使用LNB-MDT进行脂质纳米泡分析。

启动程序
--------

图形界面启动
~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color:rgb(0, 0, 0); margin-top: 0;">🖥️ 启动图形界面</h3>
   <pre style="background-color:rgb(255, 255, 255); color: #ecf0f1; padding: 15px; border-radius: 5px; overflow-x: auto;">
   <code># 激活环境
   conda activate LNB-MDT
   
   # 启动主程序
   python main.py</code>
   </pre>
   </div>

启动后您将看到LNB-MDT的主界面，包含以下功能模块：

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); color: white; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0;">🧬 Generation Module</h4>
   <p style="margin-bottom: 0;">脂质纳气泡生成</p>
   </div>

   <div style="background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); color: white; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0;">📊 Analysis Module</h4>
   <p style="margin-bottom: 0;">轨迹分析</p>
   </div>

   <div style="background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%); color: white; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0;">📈 Figure Module</h4>
   <p style="margin-bottom: 0;">数据可视化</p>
   </div>

   <div style="background: linear-gradient(135deg, #fa709a 0%, #fee140 100%); color: white; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0;">🔧 VMD Module</h4>
   <p style="margin-bottom: 0;">VMD可视化</p>
   </div>

   </div>

命令行启动
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #e65100; margin-top: 0;">💻 使用命令行工具</h3>
   <p>适合批量处理和自动化分析：</p>
   <pre style="background-color:rgb(255, 255, 255); color: #ecf0f1; padding: 15px; border-radius: 5px; overflow-x: auto;">
   <code># 激活环境
   conda activate LNB-MDT
   
   # 查看帮助信息
   python analysis/pca.py --help</code>
   </pre>
   </div>

基本分析流程
------------

步骤1：准备数据文件
~~~~~~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 15px; border-radius: 8px; border-left: 4px solid #2196f3;">

**必需文件：**

- **GRO文件**: 分子拓扑结构文件
- **XTC文件**: 分子动力学轨迹文件

**示例数据：**
项目包含示例数据文件：
- `cases/lnb.gro` - 示例拓扑文件
- `cases/md.xtc` - 示例轨迹文件


步骤2：选择分析类型
~~~~~~~~~~~~~~~~~~~~

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

步骤3：配置参数
~~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

**关键参数：**

- **残基组**: 指定要分析的分子类型和原子
- **帧范围**: 选择分析的时间范围
- **计算参数**: k值、截止距离等
- **并行处理**: 启用多核加速


步骤4：运行分析
~~~~~~~~~~~~~~~~

图形界面运行
^^^^^^^^^^^^

1. 在界面中加载GRO和XTC文件
2. 选择分析类型
3. 配置参数
4. 点击"运行"按钮
5. 查看结果

命令行运行
^^^^^^^^^^

.. code-block:: bash

   # PCA分析示例
   python analysis/pca.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/pca_results.csv \
     --residues "{'DPPC': ['PO4']}" \
     --parallel \
     --verbose

步骤5：查看结果
~~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**输出文件：**

- **CSV文件**: 包含分析结果的数值数据
- **图表**: 可视化分析结果
- **日志**: 分析过程的详细信息

**结果解读：**
- 查看CSV文件中的数值结果
- 使用图表模块可视化数据
- 结合VMD进行分子可视化

   </div>

实际示例
--------

示例1：PCA分析
~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #7b1fa2; margin-top: 0;">🧬 PCA主成分分析</h3>
   <p>分析脂质分子的构象变化：</p>
   <pre style="background-color: #2c3e50; color: #ecf0f1; padding: 15px; border-radius: 5px; overflow-x: auto;">
   <code>python analysis/pca.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/pca_test.csv \
     --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" \
     --start-frame 0 \
     --stop-frame 100 \
     --parallel \
     --verbose</code>
   </pre>
   </div>

示例2：面积分析
~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #388e3c; margin-top: 0;">📏 Voronoi面积分析</h3>
   <p>计算脂质分子的Voronoi镶嵌面积：</p>
   <pre style="background-color: #2c3e50; color: #ecf0f1; padding: 15px; border-radius: 5px; overflow-x: auto;">
   <code>python analysis/area.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/area_test.csv \
     --residues "{'DPPC': ['PO4']}" \
     --k-value 20 \
     --max-normal-angle 140 \
     --parallel \
     --verbose</code>
   </pre>
   </div>

示例3：曲率分析
~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #f57c00; margin-top: 0;">🌊 膜曲率分析</h3>
   <p>计算脂质膜的曲率特性：</p>
   <pre style="background-color: #2c3e50; color: #ecf0f1; padding: 15px; border-radius: 5px; overflow-x: auto;">
   <code>python analysis/curvature.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/curvature_test.csv \
     --residues "{'DPPC': ['PO4']}" \
     --k-value 20 \
     --method mean \
     --parallel \
     --verbose</code>
   </pre>
   </div>

机器学习功能
------------

LNB-MDT集成了强大的机器学习功能：

参数优化
~~~~~~~~

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**自动参数优化：**

.. code-block:: python

   from machine_learning import AnalysisParameterOptimizer
   
   # 创建优化器
   optimizer = AnalysisParameterOptimizer('area')
   
   # 运行优化
   results = optimizer.optimize()
   print(f"最佳参数: {results['best_parameters']}")

   </div>

异常检测
~~~~~~~~

.. raw:: html

   <div style="background-color: #fce4ec; padding: 15px; border-radius: 8px; border-left: 4px solid #e91e63;">

**异常模式检测：**

.. code-block:: python

   from machine_learning import MDAnomalyDetector
   
   # 创建检测器
   detector = MDAnomalyDetector(method='isolation_forest')
   
   # 分析轨迹
   results = detector.analyze_trajectory(
       gro_file="cases/lnb.gro",
       xtc_file="cases/md.xtc",
       residues={'DPPC': ['PO4']}
   )

   </div>

属性预测
~~~~~~~~

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**分子属性预测：**

.. code-block:: python

   from machine_learning import MDPropertyPredictor
   
   # 创建预测器
   predictor = MDPropertyPredictor(
       model_type='random_forest',
       target_property='diffusion_coefficient'
   )
   
   # 训练模型
   results = predictor.fit(X_train, y_train)

   </div>

VMD集成
--------

LNB-MDT支持与VMD的无缝集成：

启动VMD
~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**图形界面：**
1. 点击"Start VMD"按钮
2. 等待VMD启动
3. 拖拽CSV文件到VMD窗口
4. 选择分子进行可视化

**命令行：**
.. code-block:: bash

   # 启动VMD
   python -c "from modules.vmd_control import VMDTcp; vmd = VMDTcp(); vmd.start()"

   </div>

可视化操作
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**操作步骤：**
1. 在LNB-MDT中加载分析结果
2. 选择要可视化的帧和分子
3. VMD自动跳转到对应帧
4. 高亮显示选中的分子
5. 调整可视化参数

   </div>

Next Step
------

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 20px; border-radius: 8px; margin: 20px 0; text-align: center;">
   <h3 style="color: #1976d2; margin-top: 0;">🎉 恭喜！</h3>
   <p>您已经成功完成了LNB-MDT的快速开始！</p>
   <p><strong>接下来可以：</strong></p>
   <ul style="text-align: left; display: inline-block;">
   <li>📖 查看 <a href="user_guide.html">用户指南</a> 了解详细功能</li>
   <li>🔬 学习 <a href="analysis_modules.html">分析模块</a> 的深度使用</li>
   <li>🤖 探索 <a href="machine_learning.html">机器学习</a> 功能</li>
   <li>💻 掌握 <a href="command_line.html">命令行工具</a> 的高级用法</li>
   </ul>
   </div>

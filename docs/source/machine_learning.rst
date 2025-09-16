机器学习模块
============

LNB-MDT集成了强大的机器学习功能，为分子动力学分析提供智能化支持。

模块概览
--------

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">🔧 参数优化</h3>
   <p>贝叶斯优化自动寻找最佳分析参数</p>
   <ul style="margin-bottom: 0;">
   <li>高斯过程回归</li>
   <li>多目标优化</li>
   <li>并行评估</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">🔍 异常检测</h3>
   <p>识别分子动力学轨迹中的异常模式</p>
   <ul style="margin-bottom: 0;">
   <li>Isolation Forest</li>
   <li>Local Outlier Factor</li>
   <li>Elliptic Envelope</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">📊 属性预测</h3>
   <p>使用ML模型预测分子属性</p>
   <ul style="margin-bottom: 0;">
   <li>随机森林</li>
   <li>梯度提升</li>
   <li>神经网络</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">🎯 模式识别</h3>
   <p>识别分子动力学中的模式</p>
   <ul style="margin-bottom: 0;">
   <li>聚类分析</li>
   <li>分类识别</li>
   <li>特征提取</li>
   </ul>
   </div>

   </div>

安装和配置
----------

依赖安装
~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**机器学习模块需要额外的依赖包：**

.. code:: bash

   # 安装ML依赖
   pip install scikit-learn scipy matplotlib seaborn joblib
   
   # 可选：深度学习支持
   pip install tensorflow torch

   </div>

验证安装
~~~~~~~~

.. code:: python

   # 验证ML模块安装
   from machine_learning import ParameterOptimizer, AnomalyDetector, PropertyPredictor
   print("机器学习模块安装成功！")

参数优化
--------

贝叶斯优化
~~~~~~~~~~~

**功能描述**
使用贝叶斯优化自动寻找最佳分析参数，提高分析效率和准确性。

**算法原理**
- 使用高斯过程回归建模目标函数
- 通过采集函数指导搜索
- 平衡探索和利用

**使用示例**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #2c3e50; margin-top: 0;">🔧 参数优化示例</h3>
   <pre style="background-color: #2c3e50; color: #ecf0f1; padding: 15px; border-radius: 5px; overflow-x: auto;">
   <code>from machine_learning import AnalysisParameterOptimizer
   import time

   # 创建面积分析优化器
   optimizer = AnalysisParameterOptimizer('area')

   # 定义目标函数
   def objective_function(params):
       try:
           from analysis.area import Area
           
           analyzer = Area(
               gro_file="cases/lnb.gro",
               xtc_file="cases/md.xtc",
               residues={'DPPC': ['PO4']},
               **params
           )
           
           start_time = time.time()
           results = analyzer.run()
           computation_time = time.time() - start_time
           
           # 目标函数：计算时间 + 结果质量
           objective = computation_time + len(results) * 0.001
           
           return objective
           
       except Exception as e:
           print(f"错误: {e}")
           return float('inf')

   # 运行优化
   results = optimizer.optimize()

   print(f"最佳参数: {results['best_parameters']}")
   print(f"最佳得分: {results['best_score']}")</code>
   </pre>
   </div>

**优化参数类型**

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**支持的优化参数：**

- **k-value**: Voronoi镶嵌的k值 (5-50)
- **cutoff**: 聚类截止距离 (5.0-15.0)
- **n-components**: PCA主成分数量 (2-10)
- **max-normal-angle**: 最大法线角度 (120-160)

   </div>

k值优化器
~~~~~~~~~

**功能描述**
专门用于优化k值的工具，针对不同分析类型提供最优k值。

**使用示例**

.. code:: python

   from machine_learning import KValueOptimizer

   # 创建k值优化器
   optimizer = KValueOptimizer('area')

   # 运行优化
   best_k = optimizer.optimize(
       gro_file="cases/lnb.gro",
       xtc_file="cases/md.xtc",
       residues={'DPPC': ['PO4']}
   )

   print(f"最佳k值: {best_k}")

异常检测
--------

算法选择
~~~~~~~~

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background-color: #e3f2fd; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #1976d2;">🌲 Isolation Forest</h4>
   <p>适用于一般异常检测，对高维数据效果好</p>
   </div>

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #7b1fa2;">🔍 Local Outlier Factor</h4>
   <p>适用于局部异常检测，能识别密度异常</p>
   </div>

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #388e3c;">📊 Elliptic Envelope</h4>
   <p>适用于高斯分布数据的异常检测</p>
   </div>

   </div>

轨迹异常检测
~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #2c3e50; margin-top: 0;">🔍 轨迹异常检测示例</h3>
   <pre style="background-color: #2c3e50; color: #ecf0f1; padding: 15px; border-radius: 5px; overflow-x: auto;">
   <code>from machine_learning import MDAnomalyDetector
   import matplotlib.pyplot as plt

   # 创建异常检测器
   detector = MDAnomalyDetector(
       method='isolation_forest',
       contamination=0.1
   )

   # 分析轨迹
   results = detector.analyze_trajectory(
       gro_file="cases/lnb.gro",
       xtc_file="cases/md.xtc",
       residues={'DPPC': ['PO4'], 'CHOL': ['ROH']},
       start_frame=0,
       stop_frame=1000,
       step_frame=5
   )

   # 打印结果
   print(f"总帧数: {len(results['predictions'])}")
   print(f"异常帧数: {results['n_anomalies']}")
   print(f"异常比例: {results['anomaly_ratio']:.2%}")

   # 可视化结果
   detector.plot_anomalies(results, save_path="anomaly_analysis.png")

   # 分析特定异常
   anomaly_frames = results['anomaly_indices']
   print(f"异常帧: {anomaly_frames}")</code>
   </pre>
   </div>

**特征提取**

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**自动提取的特征：**

- **结构特征**: 分子间距离、角度、二面角
- **动力学特征**: 速度、加速度、扩散系数
- **热力学特征**: 能量、温度、压力
- **统计特征**: 均值、方差、相关性

   </div>

属性预测
--------

模型类型
~~~~~~~~

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background-color: #e3f2fd; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #1976d2;">🌲 随机森林</h4>
   <p style="margin-bottom: 0;">适用于非线性关系</p>
   </div>

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #7b1fa2;">🚀 梯度提升</h4>
   <p style="margin-bottom: 0;">高精度预测</p>
   </div>

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #388e3c;">🧠 神经网络</h4>
   <p style="margin-bottom: 0;">复杂模式识别</p>
   </div>

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #f57c00;">📈 支持向量机</h4>
   <p style="margin-bottom: 0;">小样本学习</p>
   </div>

   </div>

预测示例
~~~~~~~~

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #2c3e50; margin-top: 0;">📊 属性预测示例</h3>
   <pre style="background-color: #2c3e50; color: #ecf0f1; padding: 15px; border-radius: 5px; overflow-x: auto;">
   <code>from machine_learning import MDPropertyPredictor
   import numpy as np
   import pandas as pd

   # 创建预测器
   predictor = MDPropertyPredictor(
       model_type='random_forest',
       target_property='diffusion_coefficient',
       n_estimators=100,
       max_depth=10
   )

   # 生成训练数据（替换为真实数据）
   np.random.seed(42)
   n_samples = 1000
   n_features = 15

   # 创建合成特征和目标
   X = np.random.randn(n_samples, n_features)
   y = np.random.randn(n_samples) * 0.1 + 1.0  # 合成扩散系数

   # 训练模型
   results = predictor.fit(X, y, test_size=0.2)

   # 打印性能指标
   print(f"训练R²: {results['train_r2']:.4f}")
   print(f"测试R²: {results['test_r2']:.4f}")
   print(f"交叉验证均值: {results['cv_mean']:.4f}")

   # 获取特征重要性
   importance = predictor.get_feature_importance()
   print("前5个重要特征:")
   for feature, score in sorted(importance.items(), key=lambda x: x[1], reverse=True)[:5]:
       print(f"  {feature}: {score:.4f}")

   # 可视化结果
   predictor.plot_results(results, save_path="property_prediction.png")

   # 保存模型
   predictor.save_model('diffusion_predictor.pkl')</code>
   </pre>
   </div>

**可预测的属性**

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**支持的预测属性：**

- **扩散系数**: 分子扩散能力
- **渗透性**: 膜渗透性
- **稳定性**: 系统稳定性
- **相变温度**: 相变行为
- **机械性质**: 弹性模量等

   </div>

模式识别
--------

聚类分析
~~~~~~~~

**功能描述**
识别分子动力学中的聚类模式和结构特征。

**使用示例**

.. code:: python

   from machine_learning import PatternRecognizer

   # 创建模式识别器
   recognizer = PatternRecognizer(method='kmeans')

   # 分析轨迹模式
   patterns = recognizer.analyze_patterns(
       gro_file="cases/lnb.gro",
       xtc_file="cases/md.xtc",
       residues={'DPPC': ['PO4']}
   )

   # 可视化模式
   recognizer.plot_patterns(patterns)

分类识别
~~~~~~~~

**功能描述**
对分子动力学状态进行分类识别。

**使用示例**

.. code:: python

   from machine_learning import StateClassifier

   # 创建状态分类器
   classifier = StateClassifier(
       model_type='random_forest',
       n_classes=3
   )

   # 训练分类器
   classifier.fit(X_train, y_train)

   # 预测状态
   predictions = classifier.predict(X_test)

   # 评估性能
   accuracy = classifier.evaluate(X_test, y_test)
   print(f"分类准确率: {accuracy:.4f}")

数据处理
--------

特征工程
~~~~~~~~

**功能描述**
从分子动力学轨迹中提取有意义的特征。

**特征类型**

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**提取的特征类型：**

- **几何特征**: 距离、角度、体积
- **动力学特征**: 速度、加速度、扩散
- **热力学特征**: 能量、温度、压力
- **统计特征**: 均值、方差、相关性
- **拓扑特征**: 连通性、聚类系数

   </div>

**使用示例**

.. code:: python

   from machine_learning import FeatureExtractor

   # 创建特征提取器
   extractor = FeatureExtractor()

   # 提取特征
   features = extractor.extract_features(
       gro_file="cases/lnb.gro",
       xtc_file="cases/md.xtc",
       residues={'DPPC': ['PO4']}
   )

   print(f"提取特征数: {features.shape[1]}")
   print(f"特征名称: {extractor.feature_names}")

数据预处理
~~~~~~~~~~

**功能描述**
对原始数据进行清洗和预处理。

**预处理步骤**

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**预处理流程：**

1. **数据清洗**: 去除异常值和缺失值
2. **特征缩放**: 标准化和归一化
3. **特征选择**: 选择重要特征
4. **数据分割**: 训练集和测试集分割
5. **交叉验证**: 模型验证

   </div>

**使用示例**

.. code:: python

   from machine_learning import DataProcessor

   # 创建数据处理器
   processor = DataProcessor()

   # 预处理数据
   X_processed, y_processed = processor.preprocess(
       X_raw, y_raw,
       scale=True,
       select_features=True,
       test_size=0.2
   )

   print(f"处理后特征数: {X_processed.shape[1]}")

API参考
-------

参数优化器
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

**AnalysisParameterOptimizer**

.. code:: python

   class AnalysisParameterOptimizer:
       def __init__(self, analysis_type, **kwargs):
           """
           初始化参数优化器
           
           参数:
           - analysis_type: 分析类型 ('pca', 'area', 'curvature', ...)
           - n_initial_points: 初始随机点数 (默认: 10)
           - n_iterations: 优化迭代次数 (默认: 50)
           - random_state: 随机种子 (默认: 42)
           """
       
       def optimize(self, objective_function=None):
           """
           运行优化过程
           
           返回:
           - Dict包含最佳参数、最佳得分和优化历史
           """
       
       def save_model(self, filepath):
           """保存优化模型"""
       
       def load_model(self, filepath):
           """加载优化模型"""

   </div>

异常检测器
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

**MDAnomalyDetector**

.. code:: python

   class MDAnomalyDetector:
       def __init__(self, method='isolation_forest', **kwargs):
           """
           初始化异常检测器
           
           参数:
           - method: 检测方法 ('isolation_forest', 'lof', 'elliptic_envelope')
           - contamination: 预期异常比例 (默认: 0.1)
           """
       
       def analyze_trajectory(self, gro_file, xtc_file, **kwargs):
           """
           分析轨迹中的异常
           
           返回:
           - Dict包含预测结果、异常索引和统计信息
           """
       
       def plot_anomalies(self, results, save_path=None):
           """可视化异常检测结果"""

   </div>

属性预测器
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

**MDPropertyPredictor**

.. code:: python

   class MDPropertyPredictor:
       def __init__(self, model_type='random_forest', **kwargs):
           """
           初始化属性预测器
           
           参数:
           - model_type: 模型类型 ('random_forest', 'gradient_boosting', 'neural_network', 'svr')
           - target_property: 预测目标属性
           """
       
       def fit(self, X, y, test_size=0.2):
           """
           训练模型
           
           返回:
           - Dict包含训练结果和性能指标
           """
       
       def predict(self, X):
           """预测新数据"""
       
       def get_feature_importance(self):
           """获取特征重要性"""
       
       def plot_results(self, results, save_path=None):
           """可视化预测结果"""

   </div>

最佳实践
--------

参数优化建议
~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**优化建议：**

1. **定义清晰目标**: 确保目标函数同时考虑准确性和效率
2. **设置合理边界**: 使用领域知识设置参数范围
3. **监控进度**: 使用日志跟踪优化进度
4. **验证结果**: 在测试数据上验证优化参数

   </div>

异常检测建议
~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**检测建议：**

1. **选择合适方法**: Isolation Forest用于一般异常，LOF用于局部异常
2. **调整污染率**: 根据预期异常比例设置contamination
3. **特征选择**: 使用与分析相关的特征
4. **结果解释**: 在系统背景下分析检测到的异常

   </div>

属性预测建议
~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**预测建议：**

1. **特征工程**: 从轨迹数据中提取有意义的特征
2. **模型选择**: 尝试多种模型并比较性能
3. **交叉验证**: 使用交叉验证评估模型泛化能力
4. **特征重要性**: 分析特征重要性理解预测

   </div>

故障排除
--------

常见问题
~~~~~~~~

问题1：导入错误
^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #ffebee; padding: 15px; border-radius: 8px; border-left: 4px solid #f44336;">

**解决方案：**
.. code:: bash

   pip install scikit-learn scipy matplotlib seaborn joblib

   </div>

问题2：内存不足
^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**解决方案：**
- 减少批处理大小
- 使用特征选择降低维度
- 分块处理数据

   </div>

问题3：性能不佳
^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**解决方案：**
- 检查特征缩放
- 尝试不同模型类型
- 调整超参数
- 增加训练数据

   </div>

问题4：收敛问题
^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**解决方案：**
- 检查参数边界
- 调整优化参数
- 验证目标函数

   </div>

性能优化
--------

计算优化
~~~~~~~~

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**优化策略：**

- **并行处理**: 使用多核CPU加速
- **内存管理**: 优化内存使用
- **算法选择**: 选择高效算法
- **数据预处理**: 减少计算复杂度

   </div>

模型优化
~~~~~~~~

.. raw:: html

   <div style="background-color: #fce4ec; padding: 15px; border-radius: 8px; border-left: 4px solid #e91e63;">

**优化方法：**

- **超参数调优**: 使用网格搜索或随机搜索
- **特征选择**: 选择最重要的特征
- **模型集成**: 使用集成方法提高性能
- **正则化**: 防止过拟合

   </div>

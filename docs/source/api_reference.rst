API参考
========

LNB-MDT的完整API参考文档。

分析模块API
-----------

分析基类
~~~~~~~~

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

**AnalysisBase**

所有分析模块的基类，提供通用的分析功能。

   </div>

.. code-block:: python

   class AnalysisBase:
       """
       分析基类
       
       所有分析模块都继承自此类，提供通用的分析功能。
       """
       
       def __init__(self, gro_file, xtc_file, residues, **kwargs):
           """
           初始化分析器
           
           参数:
           - gro_file (str): GRO拓扑文件路径
           - xtc_file (str): XTC轨迹文件路径
           - residues (dict): 残基组字典
           - **kwargs: 其他参数
           """
       
       def run(self, start_frame=0, stop_frame=-1, step_frame=1, **kwargs):
           """
           运行分析
           
           参数:
           - start_frame (int): 起始帧
           - stop_frame (int): 结束帧
           - step_frame (int): 帧步长
           
           返回:
           - dict: 分析结果
           """
       
       def save_results(self, output_file, results):
           """
           保存分析结果
           
           参数:
           - output_file (str): 输出文件路径
           - results (dict): 分析结果
           """

PCA分析
~~~~~~~

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 15px; border-radius: 8px; border-left: 4px solid #2196f3;">

**PCA类**

主成分分析类，用于分析分子构象变化。

   </div>

.. code-block:: python

   class PCA(AnalysisBase):
       """
       主成分分析类
       
       用于分析脂质分子的构象变化和运动模式。
       """
       
       def __init__(self, gro_file, xtc_file, residues, n_components=3, **kwargs):
           """
           初始化PCA分析器
           
           参数:
           - gro_file (str): GRO拓扑文件路径
           - xtc_file (str): XTC轨迹文件路径
           - residues (dict): 残基组字典
           - n_components (int): 主成分数量，默认3
           """
       
       def run(self, start_frame=0, stop_frame=-1, step_frame=1, **kwargs):
           """
           运行PCA分析
           
           返回:
           - dict: 包含主成分值的字典
           """
       
       def get_explained_variance_ratio(self):
           """
           获取解释方差比
           
           返回:
           - numpy.ndarray: 解释方差比数组
           """

面积分析
~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**Area类**

Voronoi镶嵌面积分析类。

   </div>

.. code-block:: python

   class Area(AnalysisBase):
       """
       Voronoi镶嵌面积分析类
       
       计算脂质分子的Voronoi镶嵌面积。
       """
       
       def __init__(self, gro_file, xtc_file, residues, k_value=20, 
                    max_normal_angle=140, **kwargs):
           """
           初始化面积分析器
           
           参数:
           - gro_file (str): GRO拓扑文件路径
           - xtc_file (str): XTC轨迹文件路径
           - residues (dict): 残基组字典
           - k_value (int): Voronoi镶嵌的k值，默认20
           - max_normal_angle (float): 最大法线角度，默认140度
           """
       
       def run(self, start_frame=0, stop_frame=-1, step_frame=1, **kwargs):
           """
           运行面积分析
           
           返回:
           - dict: 包含面积值的字典
           """
       
       def calculate_voronoi_area(self, coordinates, k_value):
           """
           计算Voronoi面积
           
           参数:
           - coordinates (numpy.ndarray): 分子坐标
           - k_value (int): k值
           
           返回:
           - numpy.ndarray: Voronoi面积数组
           """

曲率分析
~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**Curvature类**

膜曲率分析类。

   </div>

.. code-block:: python

   class Curvature(AnalysisBase):
       """
       膜曲率分析类
       
       计算脂质膜的平均曲率和高斯曲率。
       """
       
       def __init__(self, gro_file, xtc_file, residues, k_value=20, 
                    method='mean', **kwargs):
           """
           初始化曲率分析器
           
           参数:
           - gro_file (str): GRO拓扑文件路径
           - xtc_file (str): XTC轨迹文件路径
           - residues (dict): 残基组字典
           - k_value (int): 曲率计算的k值，默认20
           - method (str): 曲率类型，'mean'或'gaussian'，默认'mean'
           """
       
       def run(self, start_frame=0, stop_frame=-1, step_frame=1, **kwargs):
           """
           运行曲率分析
           
           返回:
           - dict: 包含曲率值的字典
           """
       
       def calculate_mean_curvature(self, coordinates, k_value):
           """
           计算平均曲率
           
           参数:
           - coordinates (numpy.ndarray): 分子坐标
           - k_value (int): k值
           
           返回:
           - numpy.ndarray: 平均曲率数组
           """
       
       def calculate_gaussian_curvature(self, coordinates, k_value):
           """
           计算高斯曲率
           
           参数:
           - coordinates (numpy.ndarray): 分子坐标
           - k_value (int): k值
           
           返回:
           - numpy.ndarray: 高斯曲率数组
           """

高度分析
~~~~~~~~

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**Height类**

分子高度分析类。

   </div>

.. code-block:: python

   class Height(AnalysisBase):
       """
       分子高度分析类
       
       分析脂质分子的高度分布和膜厚度。
       """
       
       def __init__(self, gro_file, xtc_file, residues, k_value=20, **kwargs):
           """
           初始化高度分析器
           
           参数:
           - gro_file (str): GRO拓扑文件路径
           - xtc_file (str): XTC轨迹文件路径
           - residues (dict): 残基组字典，支持多组原子
           - k_value (int): 高度计算的k值，默认20
           """
       
       def run(self, start_frame=0, stop_frame=-1, step_frame=1, **kwargs):
           """
           运行高度分析
           
           返回:
           - dict: 包含高度值的字典
           """
       
       def calculate_height(self, coordinates, reference_atoms, k_value):
           """
           计算分子高度
           
           参数:
           - coordinates (numpy.ndarray): 分子坐标
           - reference_atoms (list): 参考原子列表
           - k_value (int): k值
           
           返回:
           - numpy.ndarray: 高度值数组
           """

聚类分析
~~~~~~~~

.. raw:: html

   <div style="background-color: #fce4ec; padding: 15px; border-radius: 8px; border-left: 4px solid #e91e63;">

**Cluster类**

分子聚类分析类。

   </div>

.. code-block:: python

   class Cluster(AnalysisBase):
       """
       分子聚类分析类
       
       分析脂质分子的聚集行为和聚类模式。
       """
       
       def __init__(self, gro_file, xtc_file, residues, cutoff=8.0, **kwargs):
           """
           初始化聚类分析器
           
           参数:
           - gro_file (str): GRO拓扑文件路径
           - xtc_file (str): XTC轨迹文件路径
           - residues (dict): 残基组字典
           - cutoff (float): 聚类截止距离，默认8.0埃
           """
       
       def run(self, start_frame=0, stop_frame=-1, step_frame=1, **kwargs):
           """
           运行聚类分析
           
           返回:
           - dict: 包含聚类信息的字典
           """
       
       def find_clusters(self, coordinates, cutoff):
           """
           寻找聚类
           
           参数:
           - coordinates (numpy.ndarray): 分子坐标
           - cutoff (float): 截止距离
           
           返回:
           - list: 聚类列表
           """

机器学习模块API
---------------

参数优化器
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**ParameterOptimizer类**

贝叶斯参数优化器。

   </div>

.. code-block:: python

   class ParameterOptimizer:
       """
       参数优化器基类
       
       使用贝叶斯优化自动寻找最佳分析参数。
       """
       
       def __init__(self, parameter_bounds, objective_function, 
                    n_initial_points=10, n_iterations=50, random_state=42):
           """
           初始化参数优化器
           
           参数:
           - parameter_bounds (dict): 参数边界字典
           - objective_function (callable): 目标函数
           - n_initial_points (int): 初始随机点数，默认10
           - n_iterations (int): 优化迭代次数，默认50
           - random_state (int): 随机种子，默认42
           """
       
       def optimize(self):
           """
           运行优化过程
           
           返回:
           - dict: 包含最佳参数、最佳得分和优化历史
           """
       
       def save_model(self, filepath):
           """
           保存优化模型
           
           参数:
           - filepath (str): 保存路径
           """
       
       def load_model(self, filepath):
           """
           加载优化模型
           
           参数:
           - filepath (str): 加载路径
           """

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

**AnalysisParameterOptimizer类**

分析参数优化器。

   </div>

.. code-block:: python

   class AnalysisParameterOptimizer(ParameterOptimizer):
       """
       分析参数优化器
       
       针对特定分析类型的参数优化器。
       """
       
       def __init__(self, analysis_type, **kwargs):
           """
           初始化分析参数优化器
           
           参数:
           - analysis_type (str): 分析类型 ('pca', 'area', 'curvature', ...)
           - **kwargs: 其他参数
           """
       
       def get_parameter_bounds(self, analysis_type):
           """
           获取参数边界
           
           参数:
           - analysis_type (str): 分析类型
           
           返回:
           - dict: 参数边界字典
           """

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**KValueOptimizer类**

k值优化器。

   </div>

.. code-block:: python

   class KValueOptimizer:
       """
       k值优化器
       
       专门用于优化k值的工具。
       """
       
       def __init__(self, analysis_type, **kwargs):
           """
           初始化k值优化器
           
           参数:
           - analysis_type (str): 分析类型
           - **kwargs: 其他参数
           """
       
       def optimize(self, gro_file, xtc_file, residues, **kwargs):
           """
           优化k值
           
           参数:
           - gro_file (str): GRO文件路径
           - xtc_file (str): XTC文件路径
           - residues (dict): 残基组字典
           
           返回:
           - int: 最佳k值
           """

异常检测器
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**AnomalyDetector类**

异常检测器基类。

   </div>

.. code-block:: python

   class AnomalyDetector:
       """
       异常检测器基类
       
       检测数据中的异常模式。
       """
       
       def __init__(self, method='isolation_forest', contamination=0.1, **kwargs):
           """
           初始化异常检测器
           
           参数:
           - method (str): 检测方法 ('isolation_forest', 'lof', 'elliptic_envelope')
           - contamination (float): 预期异常比例，默认0.1
           - **kwargs: 其他参数
           """
       
       def fit(self, data):
           """
           拟合异常检测模型
           
           参数:
           - data (numpy.ndarray): 训练数据
           """
       
       def predict(self, data):
           """
           预测异常
           
           参数:
           - data (numpy.ndarray): 测试数据
           
           返回:
           - numpy.ndarray: 预测结果 (-1: 异常, 1: 正常)
           """
       
       def predict_proba(self, data):
           """
           预测异常概率
           
           参数:
           - data (numpy.ndarray): 测试数据
           
           返回:
           - numpy.ndarray: 异常概率数组
           """

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**MDAnomalyDetector类**

分子动力学异常检测器。

   </div>

.. code-block:: python

   class MDAnomalyDetector(AnomalyDetector):
       """
       分子动力学异常检测器
       
       检测分子动力学轨迹中的异常。
       """
       
       def __init__(self, method='isolation_forest', contamination=0.1, **kwargs):
           """
           初始化MD异常检测器
           
           参数:
           - method (str): 检测方法
           - contamination (float): 预期异常比例
           - **kwargs: 其他参数
           """
       
       def analyze_trajectory(self, gro_file, xtc_file, residues, 
                             start_frame=0, stop_frame=-1, step_frame=1, **kwargs):
           """
           分析轨迹中的异常
           
           参数:
           - gro_file (str): GRO文件路径
           - xtc_file (str): XTC文件路径
           - residues (dict): 残基组字典
           - start_frame (int): 起始帧
           - stop_frame (int): 结束帧
           - step_frame (int): 帧步长
           
           返回:
           - dict: 异常检测结果
           """
       
       def plot_anomalies(self, results, save_path=None):
           """
           可视化异常检测结果
           
           参数:
           - results (dict): 异常检测结果
           - save_path (str): 保存路径
           """

属性预测器
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**PropertyPredictor类**

属性预测器基类。

   </div>

.. code-block:: python

   class PropertyPredictor:
       """
       属性预测器基类
       
       使用机器学习模型预测分子属性。
       """
       
       def __init__(self, model_type='random_forest', target_property='diffusion_coefficient', **kwargs):
           """
           初始化属性预测器
           
           参数:
           - model_type (str): 模型类型 ('random_forest', 'gradient_boosting', 'neural_network', 'svr')
           - target_property (str): 预测目标属性
           - **kwargs: 其他参数
           """
       
       def fit(self, X, y, test_size=0.2, random_state=42):
           """
           训练模型
           
           参数:
           - X (numpy.ndarray): 特征数据
           - y (numpy.ndarray): 目标数据
           - test_size (float): 测试集比例，默认0.2
           - random_state (int): 随机种子，默认42
           
           返回:
           - dict: 训练结果和性能指标
           """
       
       def predict(self, X):
           """
           预测新数据
           
           参数:
           - X (numpy.ndarray): 特征数据
           
           返回:
           - numpy.ndarray: 预测结果
           """
       
       def get_feature_importance(self):
           """
           获取特征重要性
           
           返回:
           - dict: 特征重要性字典
           """

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

**MDPropertyPredictor类**

分子动力学属性预测器。

   </div>

.. code-block:: python

   class MDPropertyPredictor(PropertyPredictor):
       """
       分子动力学属性预测器
       
       从轨迹数据预测分子属性。
       """
       
       def __init__(self, model_type='random_forest', target_property='diffusion_coefficient', **kwargs):
           """
           初始化MD属性预测器
           
           参数:
           - model_type (str): 模型类型
           - target_property (str): 预测目标属性
           - **kwargs: 其他参数
           """
       
       def extract_features(self, gro_file, xtc_file, residues, **kwargs):
           """
           从轨迹中提取特征
           
           参数:
           - gro_file (str): GRO文件路径
           - xtc_file (str): XTC文件路径
           - residues (dict): 残基组字典
           
           返回:
           - numpy.ndarray: 提取的特征
           """
       
       def plot_results(self, results, save_path=None):
           """
           可视化预测结果
           
           参数:
           - results (dict): 预测结果
           - save_path (str): 保存路径
           """

模式识别器
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**PatternRecognizer类**

模式识别器。

   </div>

.. code-block:: python

   class PatternRecognizer:
       """
       模式识别器
       
       识别分子动力学中的模式和结构特征。
       """
       
       def __init__(self, method='kmeans', **kwargs):
           """
           初始化模式识别器
           
           参数:
           - method (str): 识别方法 ('kmeans', 'dbscan', 'hierarchical')
           - **kwargs: 其他参数
           """
       
       def analyze_patterns(self, gro_file, xtc_file, residues, **kwargs):
           """
           分析轨迹模式
           
           参数:
           - gro_file (str): GRO文件路径
           - xtc_file (str): XTC文件路径
           - residues (dict): 残基组字典
           
           返回:
           - dict: 模式分析结果
           """
       
       def plot_patterns(self, patterns, save_path=None):
           """
           可视化模式
           
           参数:
           - patterns (dict): 模式分析结果
           - save_path (str): 保存路径
           """

数据处理模块API
---------------

特征提取器
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**FeatureExtractor类**

特征提取器。

   </div>

.. code-block:: python

   class FeatureExtractor:
       """
       特征提取器
       
       从分子动力学轨迹中提取有意义的特征。
       """
       
       def __init__(self, **kwargs):
           """
           初始化特征提取器
           
           参数:
           - **kwargs: 其他参数
           """
       
       def extract_features(self, gro_file, xtc_file, residues, **kwargs):
           """
           提取特征
           
           参数:
           - gro_file (str): GRO文件路径
           - xtc_file (str): XTC文件路径
           - residues (dict): 残基组字典
           
           返回:
           - numpy.ndarray: 提取的特征
           """
       
       def extract_geometric_features(self, coordinates):
           """
           提取几何特征
           
           参数:
           - coordinates (numpy.ndarray): 分子坐标
           
           返回:
           - numpy.ndarray: 几何特征
           """
       
       def extract_dynamic_features(self, coordinates, velocities):
           """
           提取动力学特征
           
           参数:
           - coordinates (numpy.ndarray): 分子坐标
           - velocities (numpy.ndarray): 分子速度
           
           返回:
           - numpy.ndarray: 动力学特征
           """

数据处理器
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**DataProcessor类**

数据处理器。

   </div>

.. code-block:: python

   class DataProcessor:
       """
       数据处理器
       
       对原始数据进行清洗和预处理。
       """
       
       def __init__(self, **kwargs):
           """
           初始化数据处理器
           
           参数:
           - **kwargs: 其他参数
           """
       
       def preprocess(self, X, y, scale=True, select_features=True, test_size=0.2):
           """
           预处理数据
           
           参数:
           - X (numpy.ndarray): 特征数据
           - y (numpy.ndarray): 目标数据
           - scale (bool): 是否缩放特征，默认True
           - select_features (bool): 是否选择特征，默认True
           - test_size (float): 测试集比例，默认0.2
           
           返回:
           - tuple: 处理后的特征和目标数据
           """
       
       def clean_data(self, data):
           """
           清洗数据
           
           参数:
           - data (numpy.ndarray): 原始数据
           
           返回:
           - numpy.ndarray: 清洗后的数据
           """
       
       def scale_features(self, X):
           """
           缩放特征
           
           参数:
           - X (numpy.ndarray): 特征数据
           
           返回:
           - numpy.ndarray: 缩放后的特征
           """

VMD控制模块API
-------------

VMD连接器
~~~~~~~~~

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**VMDTcp类**

VMD TCP连接器。

   </div>

.. code-block:: python

   class VMDTcp:
       """
       VMD TCP连接器
       
       通过TCP连接控制VMD程序。
       """
       
       def __init__(self, rctl_path, vmd_path):
           """
           初始化VMD连接器
           
           参数:
           - rctl_path (str): 远程控制脚本路径
           - vmd_path (str): VMD程序路径
           """
       
       def start(self):
           """
           启动VMD
           
           返回:
           - int: 连接状态码
           """
       
       def stop(self):
           """
           停止VMD
           """
       
       def send_command(self, command):
           """
           发送命令到VMD
           
           参数:
           - command (str): VMD命令
           
           返回:
           - str: VMD响应
           """

VMD命令
~~~~~~~

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

**VMDCommands类**

VMD命令集合。

   </div>

.. code-block:: python

   class VMDCommands:
       """
       VMD命令集合
       
       提供常用的VMD命令。
       """
       
       @staticmethod
       def loadTopology(gro_file):
           """
           加载拓扑文件
           
           参数:
           - gro_file (str): GRO文件路径
           
           返回:
           - str: VMD命令
           """
       
       @staticmethod
       def loadTrajectory(xtc_file):
           """
           加载轨迹文件
           
           参数:
           - xtc_file (str): XTC文件路径
           
           返回:
           - str: VMD命令
           """
       
       @staticmethod
       def gotoFrame(frame):
           """
           跳转到指定帧
           
           参数:
           - frame (str): 帧号
           
           返回:
           - str: VMD命令
           """
       
       @staticmethod
       def highlightResid(resids):
           """
           高亮指定残基
           
           参数:
           - resids (list): 残基ID列表
           
           返回:
           - str: VMD命令
           """
       
       @staticmethod
       def setRepresentation(style):
           """
           设置显示样式
           
           参数:
           - style (str): 显示样式
           
           返回:
           - str: VMD命令
           """
       
       @staticmethod
       def setColoringMethod(method):
           """
           设置着色方法
           
           参数:
           - method (str): 着色方法
           
           返回:
           - str: VMD命令
           """

图表模块API
----------

图表生成器
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 15px; border-radius: 8px; border-left: 4px solid #2196f3;">

**FigurePage类**

图表页面类。

   </div>

.. code-block:: python

   class FigurePage:
       """
       图表页面类
       
       提供图表生成和可视化功能。
       """
       
       @staticmethod
       def figureBtnMakeFigure(ui):
           """
           生成图表按钮处理
           
           参数:
           - ui: 用户界面对象
           """
       
       @staticmethod
       def figureBtnColor(ui):
           """
           颜色选择按钮处理
           
           参数:
           - ui: 用户界面对象
           """
       
       @staticmethod
       def figureBtnShape(ui):
           """
           形状选择按钮处理
           
           参数:
           - ui: 用户界面对象
           """

工具函数API
----------

文件工具
~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**BtnGetPath类**

文件路径获取工具。

   </div>

.. code-block:: python

   class BtnGetPath:
       """
       文件路径获取工具
       
       提供文件选择对话框功能。
       """
       
       @staticmethod
       def run(edit_widget, file_type):
           """
           运行文件选择对话框
           
           参数:
           - edit_widget: 编辑控件
           - file_type (str): 文件类型
           """

分析工具
~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**NextClick函数**

下一步按钮处理函数。

   </div>

.. code-block:: python

   def NextClick(ui):
       """
       下一步按钮处理函数
       
       参数:
       - ui: 用户界面对象
       """

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**BtnGeneClick函数**

生成按钮处理函数。

   </div>

.. code-block:: python

   def BtnGeneClick(ui):
       """
       生成按钮处理函数
       
       参数:
       - ui: 用户界面对象
       """

.. raw:: html

   <div style="background-color: #fce4ec; padding: 15px; border-radius: 8px; border-left: 4px solid #e91e63;">

**lipidsSelect函数**

脂质选择函数。

   </div>

.. code-block:: python

   def lipidsSelect(ui):
       """
       脂质选择函数
       
       参数:
       - ui: 用户界面对象
       """

错误处理API
----------

异常类
~~~~~~

.. raw:: html

   <div style="background-color: #ffebee; padding: 15px; border-radius: 8px; border-left: 4px solid #f44336;">

**LNBException类**

LNB-MDT自定义异常类。

   </div>

.. code-block:: python

   class LNBException(Exception):
       """
       LNB-MDT自定义异常类
       
       用于处理LNB-MDT特定的错误。
       """
       
       def __init__(self, message, error_code=None):
           """
           初始化异常
           
           参数:
           - message (str): 错误消息
           - error_code (str): 错误代码
           """
           super().__init__(message)
           self.error_code = error_code

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**AnalysisError类**

分析错误异常类。

   </div>

.. code-block:: python

   class AnalysisError(LNBException):
       """
       分析错误异常类
       
       用于处理分析过程中的错误。
       """
       
       def __init__(self, message, analysis_type=None):
           """
           初始化分析错误
           
           参数:
           - message (str): 错误消息
           - analysis_type (str): 分析类型
           """
           super().__init__(message)
           self.analysis_type = analysis_type

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**VMDError类**

VMD错误异常类。

   </div>

.. code-block:: python

   class VMDError(LNBException):
       """
       VMD错误异常类
       
       用于处理VMD相关的错误。
       """
       
       def __init__(self, message, vmd_command=None):
           """
           初始化VMD错误
           
           参数:
           - message (str): 错误消息
           - vmd_command (str): VMD命令
           """
           super().__init__(message)
           self.vmd_command = vmd_command

使用示例
--------

基本使用
~~~~~~~~

.. code-block:: python

   # 基本分析使用示例
   from analysis.pca import PCA
   
   # 创建PCA分析器
   analyzer = PCA(
       gro_file="cases/lnb.gro",
       xtc_file="cases/md.xtc",
       residues={'DPPC': ['PO4']},
       n_components=3
   )
   
   # 运行分析
   results = analyzer.run(start_frame=0, stop_frame=100)
   
   # 保存结果
   analyzer.save_results("results/pca_results.csv", results)

机器学习使用
~~~~~~~~~~~~

.. code-block:: python

   # 机器学习使用示例
   from machine_learning import AnalysisParameterOptimizer, MDAnomalyDetector
   
   # 参数优化
   optimizer = AnalysisParameterOptimizer('area')
   results = optimizer.optimize()
   
   # 异常检测
   detector = MDAnomalyDetector(method='isolation_forest')
   anomalies = detector.analyze_trajectory(
       gro_file="cases/lnb.gro",
       xtc_file="cases/md.xtc",
       residues={'DPPC': ['PO4']}
   )

VMD集成使用
~~~~~~~~~~~

.. code-block:: python

   # VMD集成使用示例
   from modules.vmd_control import VMDTcp, VMDCommands
   
   # 创建VMD连接
   vmd = VMDTcp("./remote_ctl.tcl", "C:/Program Files/VMD/vmd.exe")
   
   # 启动VMD
   vmd.start()
   
   # 发送命令
   vmd.send_command(VMDCommands.loadTopology("cases/lnb.gro"))
   vmd.send_command(VMDCommands.loadTrajectory("cases/md.xtc"))
   vmd.send_command(VMDCommands.gotoFrame("100"))
   
   # 停止VMD
   vmd.stop()

错误处理使用
~~~~~~~~~~~~

.. code-block:: python

   # 错误处理使用示例
   from _exception import AnalysisError, VMDError
   
   try:
       # 运行分析
       analyzer = PCA(gro_file="invalid.gro", xtc_file="invalid.xtc", residues={})
       results = analyzer.run()
   except AnalysisError as e:
       print(f"分析错误: {e}")
       print(f"分析类型: {e.analysis_type}")
   except VMDError as e:
       print(f"VMD错误: {e}")
       print(f"VMD命令: {e.vmd_command}")
   except Exception as e:
       print(f"未知错误: {e}")

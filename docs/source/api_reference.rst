API Reference
=============

Complete API reference documentation for LNB-MDT.

Analysis Module API
-------------------

Analysis Base Class
~~~~~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

**AnalysisBase**

Base class for all analysis modules, providing common analysis functionality.

   </div>

.. code:: python

   class AnalysisBase:
       """
       Analysis Base Class
       
       All analysis modules inherit from this class, providing common analysis functionality.
       """
       
       def __init__(self, gro_file, xtc_file, residues, **kwargs):
           """
           Initialize analyzer
           
           Parameters:
           - gro_file (str): GRO topology file path
           - xtc_file (str): XTC trajectory file path
           - residues (dict): Residue group dictionary
           - **kwargs: Other parameters
           """
       
       def run(self, start_frame=0, stop_frame=-1, step_frame=1, **kwargs):
           """
           Run analysis
           
           Parameters:
           - start_frame (int): Start frame
           - stop_frame (int): Stop frame
           - step_frame (int): Frame step
           
           Returns:
           - dict: Analysis results
           """
       
       def save_results(self, output_file, results):
           """
           Save analysis results
           
           Parameters:
           - output_file (str): Output file path
           - results (dict): Analysis results
           """

PCA Analysis
~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 15px; border-radius: 8px; border-left: 4px solid #2196f3;">

**PCA Class**

Principal Component Analysis class for analyzing molecular conformational changes.

   </div>

.. code:: python

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

.. code:: python

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

.. code:: python

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

.. code:: python

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

.. code:: python

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


VMD控制模块API
--------------

VMD连接器
~~~~~~~~~

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**VMDTcp类**

VMD TCP连接器。

   </div>

.. code:: python

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

.. code:: python

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
-----------

图表生成器
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 15px; border-radius: 8px; border-left: 4px solid #2196f3;">

**FigurePage类**

图表页面类。

   </div>

.. code:: python

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
-----------

文件工具
~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**BtnGetPath类**

文件路径获取工具。

   </div>

.. code:: python

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

.. code:: python

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

.. code:: python

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

.. code:: python

   def lipidsSelect(ui):
       """
       脂质选择函数
       
       参数:
       - ui: 用户界面对象
       """

Error Handling API
------------------

Exception Classes
~~~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #ffebee; padding: 15px; border-radius: 8px; border-left: 4px solid #f44336;">

**LNBException Class**

LNB-MDT custom exception class.

   </div>

.. code:: python

   class LNBException(Exception):
       """
       LNB-MDT Custom Exception Class
       
       Used to handle LNB-MDT specific errors.
       """
       
       def __init__(self, message, error_code=None):
           """
           Initialize exception
           
           Parameters:
           - message (str): Error message
           - error_code (str): Error code
           """
           super().__init__(message)
           self.error_code = error_code

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**AnalysisError类**

分析错误异常类。

   </div>

.. code:: python

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

.. code:: python

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

Usage Examples
--------------

Basic Usage
~~~~~~~~~~~

.. code:: python

   # Basic analysis usage example
   from analysis.pca import PCA
   
   # Create PCA analyzer
   analyzer = PCA(
       gro_file="cases/lnb.gro",
       xtc_file="cases/md.xtc",
       residues={'DPPC': ['PO4']},
       n_components=3
   )
   
   # Run analysis
   results = analyzer.run(start_frame=0, stop_frame=100)
   
   # Save results
   analyzer.save_results("results/pca_results.csv", results)


VMD Integration Usage
~~~~~~~~~~~~~~~~~~~~

.. code:: python

   # VMD integration usage example
   from modules.vmd_control import VMDTcp, VMDCommands
   
   # Create VMD connection
   vmd = VMDTcp("./remote_ctl.tcl", "C:/Program Files/VMD/vmd.exe")
   
   # Start VMD
   vmd.start()
   
   # Send commands
   vmd.send_command(VMDCommands.loadTopology("cases/lnb.gro"))
   vmd.send_command(VMDCommands.loadTrajectory("cases/md.xtc"))
   vmd.send_command(VMDCommands.gotoFrame("100"))
   
   # Stop VMD
   vmd.stop()

Error Handling Usage
~~~~~~~~~~~~~~~~~~~

.. code:: python

   # Error handling usage example
   from _exception import AnalysisError, VMDError
   
   try:
       # Run analysis
       analyzer = PCA(gro_file="invalid.gro", xtc_file="invalid.xtc", residues={})
       results = analyzer.run()
   except AnalysisError as e:
       print(f"Analysis error: {e}")
       print(f"Analysis type: {e.analysis_type}")
   except VMDError as e:
       print(f"VMD error: {e}")
       print(f"VMD command: {e.vmd_command}")
   except Exception as e:
       print(f"Unknown error: {e}")

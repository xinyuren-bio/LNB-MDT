LNB-MDT 文档
============

.. image:: https://img.shields.io/badge/Version-v1.0-blue.svg
   :target: https://github.com/xinyuren-bio/LNB-MDT
   :alt: Version

.. image:: https://img.shields.io/badge/Python-3.11+-green.svg
   :target: https://python.org
   :alt: Python

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License

**LNB-MDT** (Lipid NanoBubble Molecular Dynamics Toolbox) 是一个专为脂质纳米泡分子动力学分析设计的综合性工具箱。

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #2c3e50; margin-top: 0;">🚀 快速开始</h3>
   <p>在5分钟内开始使用LNB-MDT：</p>
   <pre style="background-color: #2c3e50; color: #ecf0f1; padding: 15px; border-radius: 5px; overflow-x: auto;">
   <code># 安装
   git clone https://github.com/xinyuren-bio/LNB-MDT.git
   cd LNB-MDT
   ./install.sh
   
   # 启动图形界面
   conda activate LNB-MDT
   python main.py
   
   # 命令行分析示例
   python analysis/pca.py --gro-file cases/lnb.gro --xtc-file cases/md.xtc --residues "{'DPPC': ['PO4']}" --parallel</code>
   </pre>
   </div>

主要特性
--------

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 10px;">
   <h4 style="margin-top: 0;">🧬 分子动力学分析</h4>
   <ul style="margin-bottom: 0;">
   <li>PCA主成分分析</li>
   <li>Voronoi镶嵌面积计算</li>
   <li>曲率分析（平均/高斯）</li>
   <li>高度分布分析</li>
   <li>聚类行为分析</li>
   <li>各向异性计算</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); color: white; padding: 20px; border-radius: 10px;">
   <h4 style="margin-top: 0;">🤖 机器学习集成</h4>
   <ul style="margin-bottom: 0;">
   <li>贝叶斯参数优化</li>
   <li>异常模式检测</li>
   <li>分子属性预测</li>
   <li>高级特征工程</li>
   <li>模型性能评估</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); color: white; padding: 20px; border-radius: 10px;">
   <h4 style="margin-top: 0;">🖥️ 现代用户界面</h4>
   <ul style="margin-bottom: 0;">
   <li>Qt6图形界面</li>
   <li>直观数据可视化</li>
   <li>拖拽文件操作</li>
   <li>VMD集成支持</li>
   <li>并行处理支持</li>
   </ul>
   </div>

   </div>

系统要求
--------

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**操作系统支持：**
- Windows 10/11 (64-bit)
- macOS 10.15+
- Ubuntu 18.04+ 或其他主流Linux发行版

**软件依赖：**
- Python 3.11+
- Conda (Miniconda或Anaconda)
- VMD 1.9.4+ (可选，用于可视化)

**硬件要求：**
- 内存：最低8GB，推荐16GB+
- 存储：至少2GB可用空间
- 处理器：支持并行计算的多核处理器

   </div>

安装指南
--------

.. raw:: html

   <div style="background-color: #fff3cd; padding: 15px; border-radius: 8px; border-left: 4px solid #ffc107;">

**方法1：使用安装脚本（推荐）**

.. code-block:: bash

   # Linux/macOS
   git clone https://github.com/xinyuren-bio/LNB-MDT.git
   cd LNB-MDT
   ./install.sh

   # Windows
   git clone https://github.com/xinyuren-bio/LNB-MDT.git
   cd LNB-MDT
   install.bat

**方法2：手动安装**

.. code-block:: bash

   # 创建环境
   conda create -n LNB-MDT python=3.11
   conda activate LNB-MDT
   
   # 安装依赖
   pip install -r requirements.txt

   </div>

.. toctree::
   :maxdepth: 3
   :hidden:

   installation
   quickstart
   user_guide
   analysis_modules
   machine_learning
   command_line
   api_reference
   examples
   troubleshooting
   contributing


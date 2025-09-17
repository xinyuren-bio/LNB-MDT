安装指南
========

Conda安装
---------

安装LNB-MDT最简单的方法是通过Conda环境管理器：

.. code-block:: bash

   conda config --add channels conda-forge
   conda create -n LNB-MDT -c conda-forge python=3.11
   conda activate LNB-MDT
   git clone https://github.com/xinyuren-bio/LNB-MDT.git
   cd LNB-MDT
   pip install -r requirements.txt

这将在一个新的虚拟环境中安装LNB-MDT及其所有依赖项。

如果您还没有安装Conda，我们建议下载并安装Miniconda，这是Conda的轻量级版本。

PyPI安装
--------

也可以从Python包索引安装LNB-MDT：

.. code-block:: bash

   pip install LNB-MDT

或者安装开发版本：

.. code-block:: bash

   pip install https://github.com/xinyuren-bio/LNB-MDT/archive/main.zip

依赖项
------

LNB-MDT使用MDAnalysis进行所有分析计算，使用NumPy进行数值计算。

如上所述，安装这些包以及LNB-MDT的最简单方法是使用Conda。但是，也可以使用pip或从源码安装MDAnalysis和NumPy。

系统要求
--------

- **Windows**: Windows 10/11 (64位)
- **macOS**: macOS 10.15 (Catalina) 或更高版本
- **Linux**: Ubuntu 18.04+, CentOS 7+, 或其他主流发行版

安装完成
--------

恭喜！您已成功安装LNB-MDT。现在可以开始使用这个强大的分子动力学分析工具箱了。

下一步：查看 :doc:`quickstart` 学习基本用法。

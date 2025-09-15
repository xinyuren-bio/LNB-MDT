安装指南
==========

本指南将帮助您在不同操作系统上安装LNB-MDT。

系统要求
--------

操作系统支持
~~~~~~~~~~~~

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background-color: #e3f2fd; padding: 15px; border-radius: 8px; border-left: 4px solid #2196f3;">
   <h4 style="margin-top: 0; color: #1976d2;">🪟 Windows</h4>
   <ul style="margin-bottom: 0;">
   <li>Windows 10/11 (64-bit)</li>
   <li>PowerShell 5.0+</li>
   <li>Git for Windows</li>
   </ul>
   </div>

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">
   <h4 style="margin-top: 0; color: #7b1fa2;">🍎 macOS</h4>
   <ul style="margin-bottom: 0;">
   <li>macOS 10.15 (Catalina)+</li>
   <li>Xcode Command Line Tools</li>
   <li>Homebrew (推荐)</li>
   </ul>
   </div>

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">
   <h4 style="margin-top: 0; color: #388e3c;">🐧 Linux</h4>
   <ul style="margin-bottom: 0;">
   <li>Ubuntu 18.04+</li>
   <li>CentOS 7+</li>
   <li>其他主流发行版</li>
   </ul>
   </div>

   </div>

软件依赖
~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**必需软件：**

- **Python**: 3.11 或更高版本
- **Conda**: Miniconda 或 Anaconda
- **Git**: 用于克隆仓库

**可选软件：**

- **VMD**: 1.9.4+ (用于分子可视化)
- **Visual Studio Code**: 推荐的代码编辑器

   </div>

硬件要求
~~~~~~~~

.. raw:: html

   <div style="background-color: #fce4ec; padding: 15px; border-radius: 8px; border-left: 4px solid #e91e63;">

**最低配置：**
- CPU: 双核处理器
- 内存: 8GB RAM
- 存储: 2GB 可用空间

**推荐配置：**
- CPU: 四核或更多处理器
- 内存: 16GB+ RAM
- 存储: 5GB+ 可用空间
- GPU: 支持CUDA的显卡（可选，用于加速）

   </div>

安装方法
--------

方法1：使用安装脚本（推荐）
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

这是最简单的安装方法，安装脚本会自动处理所有依赖。

Windows安装
^^^^^^^^^^^^

1. **下载并安装Git** (如果尚未安装)
   - 访问 https://git-scm.com/download/win
   - 下载并运行安装程序

2. **克隆仓库并运行安装脚本**
   
   .. code-block:: cmd

      git clone https://github.com/xinyuren-bio/LNB-MDT.git
      cd LNB-MDT
      install.bat

3. **等待安装完成**
   - 脚本会自动创建conda环境
   - 安装所有必需的Python包
   - 验证安装是否成功

macOS/Linux安装
^^^^^^^^^^^^^^^^

1. **确保已安装Git**
   
   .. code-block:: bash

      # macOS (使用Homebrew)
      brew install git
      
      # Ubuntu/Debian
      sudo apt update && sudo apt install git

2. **克隆仓库并运行安装脚本**
   
   .. code-block:: bash

      git clone https://github.com/xinyuren-bio/LNB-MDT.git
      cd LNB-MDT
      chmod +x install.sh
      ./install.sh

3. **等待安装完成**

方法2：手动安装
~~~~~~~~~~~~~~~~

如果您需要更多控制或遇到安装脚本问题，可以手动安装。

步骤1：安装Conda
^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**选择Conda发行版：**

- **Miniconda**: 轻量级，只包含conda和Python
- **Anaconda**: 完整版，包含大量科学计算包

**下载链接：**
- Miniconda: https://docs.conda.io/en/latest/miniconda.html
- Anaconda: https://www.anaconda.com/products/distribution

   </div>

步骤2：创建虚拟环境
^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # 创建新的conda环境
   conda create -n LNB-MDT python=3.11 -y
   
   # 激活环境
   conda activate LNB-MDT

步骤3：克隆仓库
^^^^^^^^^^^^^^^^

.. code-block:: bash

   git clone https://github.com/xinyuren-bio/LNB-MDT.git
   cd LNB-MDT

步骤4：安装依赖
^^^^^^^^^^^^^^^^

.. code-block:: bash

   # 安装基础依赖
   pip install -r requirements.txt
   
   # 安装机器学习依赖（可选）
   pip install scikit-learn scipy matplotlib seaborn joblib

步骤5：验证安装
^^^^^^^^^^^^^^^^

.. code-block:: bash

   # 检查Python版本
   python --version
   
   # 检查关键依赖
   python -c "import MDAnalysis, numpy, pandas, PySide6; print('所有依赖安装成功！')"
   
   # 测试主程序
   python main.py --version

VMD集成安装（可选）
~~~~~~~~~~~~~~~~~~~~

VMD用于分子可视化，安装后可实现与LNB-MDT的无缝集成。

Windows VMD安装
^^^^^^^^^^^^^^^^

1. **下载VMD**
   - 访问 https://www.ks.uiuc.edu/Research/vmd/
   - 下载Windows版本

2. **安装VMD**
   - 运行安装程序
   - 记住安装路径（通常是 `C:\Program Files\VMD\`）

3. **配置LNB-MDT**
   - 在LNB-MDT界面中设置VMD路径
   - 或修改 `main.py` 中的 `vmd_path` 变量

macOS VMD安装
^^^^^^^^^^^^^

.. code-block:: bash

   # 使用Homebrew安装
   brew install --cask vmd
   
   # 或手动下载安装
   # 访问 https://www.ks.uiuc.edu/Research/vmd/

Linux VMD安装
^^^^^^^^^^^^^^

.. code-block:: bash

   # Ubuntu/Debian
   wget https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.4.bin.LINUXAMD64-CUDA8-OptiX4-OSPRay111p1.opengl.tar.gz
   tar -xzf vmd-1.9.4.bin.LINUXAMD64-CUDA8-OptiX4-OSPRay111p1.opengl.tar.gz
   cd vmd-1.9.4
   ./configure
   cd src
   make install

故障排除
--------

常见问题及解决方案
~~~~~~~~~~~~~~~~~~

问题1：conda命令未找到
^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #ffebee; padding: 15px; border-radius: 8px; border-left: 4px solid #f44336;">

**解决方案：**

1. 确保conda已正确安装
2. 重新启动终端
3. 手动添加到PATH环境变量：

   - Windows: 添加 `C:\Users\YourName\miniconda3\Scripts` 到PATH
   - macOS/Linux: 添加 `~/miniconda3/bin` 到PATH

   </div>

问题2：Python包安装失败
^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**解决方案：**

1. 更新pip: `pip install --upgrade pip`
2. 使用conda安装: `conda install package_name`
3. 使用国内镜像: `pip install -i https://pypi.tuna.tsinghua.edu.cn/simple package_name`

   </div>

问题3：VMD连接失败
^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**解决方案：**

1. 检查VMD是否正确安装
2. 确认VMD路径设置正确
3. 检查防火墙设置
4. 尝试手动启动VMD

   </div>

问题4：内存不足
^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**解决方案：**

1. 关闭其他应用程序
2. 使用较小的数据集进行测试
3. 调整分析参数（减少帧数）
4. 使用并行处理选项

   </div>

获取帮助
--------

如果您在安装过程中遇到问题：

1. **查看日志文件**: 检查安装脚本生成的日志
2. **检查系统要求**: 确保满足所有系统要求
3. **搜索已知问题**: 查看GitHub Issues
4. **联系支持**: 发送邮件至 zy2310205@buaa.edu.cn

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 20px; border-radius: 8px; margin: 20px 0; text-align: center;">
   <h3 style="color: #1976d2; margin-top: 0;">🎉 安装完成！</h3>
   <p>恭喜您成功安装LNB-MDT！现在可以开始使用这个强大的分子动力学分析工具箱了。</p>
   <p><strong>下一步：</strong> 查看 <a href="quickstart.html">快速开始指南</a> 学习基本使用方法。</p>
   </div>

故障排除
==========

本页面提供LNB-MDT常见问题的解决方案。

常见问题
--------

安装问题
~~~~~~~~

问题1：conda命令未找到
^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #ffebee; padding: 15px; border-radius: 8px; border-left: 4px solid #f44336;">

**症状：**
- 运行 `conda` 命令时提示 "command not found"
- 安装脚本无法执行

**解决方案：**

1. **检查conda安装**
   .. code:: bash

      # 检查conda是否安装
      which conda
      conda --version

2. **重新安装conda**
   - 访问 https://docs.conda.io/en/latest/miniconda.html
   - 下载适合您操作系统的版本
   - 重新安装conda

3. **添加到PATH环境变量**
   - Windows: 添加 `C:\Users\YourName\miniconda3\Scripts` 到PATH
   - macOS/Linux: 添加 `~/miniconda3/bin` 到PATH

4. **重新启动终端**
   - 关闭所有终端窗口
   - 重新打开终端
   - 测试conda命令

   </div>

问题2：Python包安装失败
^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**症状：**
- pip安装包时出现错误
- 依赖包安装不完整

**解决方案：**

1. **更新pip**
   .. code:: bash

      pip install --upgrade pip

2. **使用conda安装**
   .. code:: bash

      conda install package_name

3. **使用国内镜像**
   .. code:: bash

      pip install -i https://pypi.tuna.tsinghua.edu.cn/simple package_name

4. **清理缓存**
   .. code:: bash

      pip cache purge
      conda clean --all

5. **检查网络连接**
   - 确保网络连接正常
   - 检查防火墙设置
   - 尝试使用VPN

   </div>

问题3：VMD安装问题
^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**症状：**
- VMD无法启动
- VMD路径设置错误

**解决方案：**

1. **检查VMD安装**
   - 确认VMD已正确安装
   - 检查安装路径是否正确

2. **设置VMD路径**
   - Windows: `C:\Program Files\VMD\vmd.exe`
   - macOS: `/Applications/VMD.app/Contents/MacOS/VMD`
   - Linux: `/usr/local/bin/vmd`

3. **检查权限**
   - 确保有执行VMD的权限
   - 以管理员身份运行（Windows）

4. **重新安装VMD**
   - 卸载现有VMD
   - 重新下载安装VMD
   - 确保选择正确的版本

   </div>

运行问题
~~~~~~~~

问题1：程序启动失败
^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #ffebee; padding: 15px; border-radius: 8px; border-left: 4px solid #f44336;">

**症状：**
- 运行 `python main.py` 时出现错误
- 图形界面无法启动

**解决方案：**

1. **检查Python版本**
   .. code:: bash

      python --version
      # 确保版本 >= 3.11

2. **检查依赖包**
   .. code:: bash

      python -c "import MDAnalysis, numpy, pandas, PySide6; print('所有依赖安装成功！')"

3. **重新安装依赖**
   .. code:: bash

      pip install -r requirements.txt

4. **检查环境变量**
   - 确保Python路径正确
   - 检查PYTHONPATH设置

5. **查看详细错误信息**
   .. code:: bash

      python main.py --verbose

   </div>

问题2：文件加载失败
^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**症状：**
- GRO或XTC文件无法加载
- 文件格式错误

**解决方案：**

1. **检查文件格式**
   - 确认GRO文件格式正确
   - 确认XTC文件格式正确
   - 检查文件是否损坏

2. **检查文件路径**
   - 确保文件路径正确
   - 检查文件是否存在
   - 使用绝对路径

3. **检查文件权限**
   - 确保有读取文件的权限
   - 检查文件是否被其他程序占用

4. **验证文件内容**
   .. code:: bash

      # 检查GRO文件
      head -5 your_file.gro
      
      # 检查XTC文件
      file your_file.xtc

   </div>

问题3：内存不足
^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**症状：**
- 分析过程中出现内存错误
- 系统运行缓慢

**解决方案：**

1. **检查系统内存**
   .. code:: bash

      # Linux/macOS
      free -h
      
      # Windows
      wmic memorychip get size

2. **减少并行数**
   .. code:: bash

      python analysis/pca.py --n-jobs 2  # 减少并行数

3. **分段处理**
   .. code:: bash

      # 分段处理大轨迹
      python analysis/pca.py --start-frame 0 --stop-frame 1000
      python analysis/pca.py --start-frame 1000 --stop-frame 2000

4. **关闭其他程序**
   - 关闭不必要的应用程序
   - 释放系统内存

5. **使用交换文件**
   - 增加虚拟内存
   - 使用SSD作为交换空间

   </div>

分析问题
~~~~~~~~

问题1：分析结果异常
^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**症状：**
- 分析结果数值异常
- 结果不符合预期

**解决方案：**

1. **检查参数设置**
   - 验证残基组格式
   - 检查k值设置
   - 确认帧范围

2. **验证数据质量**
   - 检查轨迹质量
   - 确认拓扑文件正确
   - 验证时间步长

3. **使用示例数据测试**
   .. code:: bash

      # 使用示例数据测试
      python analysis/pca.py --gro-file cases/lnb.gro --xtc-file cases/md.xtc --residues "{'DPPC': ['PO4']}"

4. **调整参数**
   - 尝试不同的k值
   - 调整截止距离
   - 修改分析范围

5. **检查日志信息**
   - 启用verbose模式
   - 查看详细错误信息

   </div>

问题2：分析速度很慢
^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**症状：**
- 分析过程耗时很长
- 系统响应缓慢

**解决方案：**

1. **启用并行处理**
   .. code:: bash

      python analysis/pca.py --parallel --n-jobs 4

2. **优化参数**
   - 减少k值
   - 调整截止距离
   - 限制分析帧数

3. **使用SSD存储**
   - 将轨迹文件放在SSD上
   - 提高I/O性能

4. **增加系统资源**
   - 增加内存
   - 使用更快的CPU
   - 优化系统设置

5. **分段处理**
   - 将大轨迹分段处理
   - 减少单次处理的数据量

   </div>

问题3：参数格式错误
^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**症状：**
- 参数解析错误
- 残基组格式不正确

**解决方案：**

1. **检查残基组格式**
   .. code:: bash

      # 正确格式
      --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}"
      
      # 错误格式
      --residues {'DPPC': ['PO4']}  # 缺少引号

2. **验证参数类型**
   - 确保数值参数为数字
   - 检查字符串参数格式

3. **使用引号包围路径**
   .. code:: bash

      # 包含空格的路径
      --gro-file "/path with spaces/file.gro"

4. **检查特殊字符**
   - 避免使用特殊字符
   - 使用标准ASCII字符

   </div>

VMD集成问题
~~~~~~~~~~~

问题1：VMD连接失败
^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**症状：**
- VMD无法启动
- 连接VMD失败

**解决方案：**

1. **检查VMD安装**
   - 确认VMD已正确安装
   - 检查VMD版本

2. **设置VMD路径**
   - 在LNB-MDT中设置正确的VMD路径
   - 检查路径是否存在

3. **检查防火墙**
   - 确保防火墙允许VMD通信
   - 检查端口是否被占用

4. **手动启动VMD**
   - 先手动启动VMD
   - 再在LNB-MDT中连接

5. **检查权限**
   - 确保有启动VMD的权限
   - 以管理员身份运行

   </div>

问题2：VMD命令执行失败
^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #fce4ec; padding: 15px; border-radius: 8px; border-left: 4px solid #e91e63;">

**症状：**
- VMD命令无法执行
- 可视化效果异常

**解决方案：**

1. **检查命令格式**
   - 确保VMD命令格式正确
   - 检查命令语法

2. **验证文件路径**
   - 确保文件路径正确
   - 检查文件是否存在

3. **检查VMD状态**
   - 确认VMD正在运行
   - 检查连接状态

4. **重启VMD**
   - 停止VMD
   - 重新启动VMD
   - 重新连接

   </div>

机器学习问题
~~~~~~~~~~~~

问题1：ML模块导入失败
^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #ffebee; padding: 15px; border-radius: 8px; border-left: 4px solid #f44336;">

**症状：**
- 无法导入机器学习模块
- ML功能不可用

**解决方案：**

1. **安装ML依赖**
   .. code:: bash

      pip install scikit-learn scipy matplotlib seaborn joblib

2. **检查Python版本**
   - 确保Python版本 >= 3.8
   - 检查兼容性

3. **重新安装依赖**
   .. code:: bash

      pip uninstall scikit-learn
      pip install scikit-learn

4. **检查环境**
   - 确保在正确的conda环境中
   - 检查环境变量

   </div>

问题2：优化过程失败
^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**症状：**
- 参数优化失败
- 优化过程异常

**解决方案：**

1. **检查目标函数**
   - 确保目标函数正确
   - 检查返回值类型

2. **调整优化参数**
   - 减少迭代次数
   - 调整参数边界

3. **检查数据质量**
   - 确保输入数据正确
   - 检查数据格式

4. **使用简单测试**
   - 先用简单数据测试
   - 逐步增加复杂度

   </div>

性能问题
~~~~~~~~

问题1：系统资源不足
^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**症状：**
- 系统运行缓慢
- 资源使用率过高

**解决方案：**

1. **监控系统资源**
   .. code:: bash

      # Linux/macOS
      top
      htop
      
      # Windows
      taskmgr

2. **优化并行设置**
   - 减少并行数
   - 调整批处理大小

3. **清理系统**
   - 清理临时文件
   - 释放磁盘空间

4. **升级硬件**
   - 增加内存
   - 使用SSD
   - 升级CPU

   </div>

问题2：I/O性能问题
^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**症状：**
- 文件读写缓慢
- 磁盘使用率高

**解决方案：**

1. **使用SSD存储**
   - 将轨迹文件放在SSD上
   - 提高I/O性能

2. **优化文件系统**
   - 使用NTFS或ext4
   - 避免网络文件系统

3. **减少I/O操作**
   - 批量处理文件
   - 减少文件访问次数

4. **使用压缩**
   - 压缩轨迹文件
   - 减少存储空间

   </div>

调试技巧
--------

日志记录
~~~~~~~~

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #1976d2; margin-top: 0;">📝 启用详细日志</h3>
   <p>使用详细日志模式获取更多调试信息：</p>
   </div>

**启用verbose模式**

.. code:: bash

   # 命令行详细输出
   python analysis/pca.py --verbose

   # 图形界面调试模式
   python main.py --debug

**日志文件位置**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- **系统日志**: `~/.lnb-mdt/logs/`
- **分析日志**: `results/logs/`
- **错误日志**: `~/.lnb-mdt/errors/`

   </div>

**自定义日志**

.. code:: python

   import logging

   # 设置日志
   logging.basicConfig(
       level=logging.DEBUG,
       format='%(asctime)s - %(levelname)s - %(message)s',
       handlers=[
           logging.FileHandler('debug.log'),
           logging.StreamHandler()
       ]
   )

   logger = logging.getLogger(__name__)
   logger.debug("调试信息")

错误追踪
~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #f57c00; margin-top: 0;">🔍 错误追踪技巧</h3>
   <p>使用以下技巧追踪和解决错误：</p>
   </div>

**Python错误追踪**

.. code:: python

   import traceback
   import sys

   try:
       # 您的代码
       pass
   except Exception as e:
       print(f"错误: {e}")
       traceback.print_exc()
       sys.exit(1)

**系统错误检查**

.. code:: bash

   # 检查系统错误
   dmesg | tail -20  # Linux
   
   # 检查Python错误
   python -c "import sys; print(sys.version)"

**网络连接检查**

.. code:: bash

   # 检查网络连接
   ping google.com
   
   # 检查端口
   netstat -an | grep 8080

性能分析
~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #388e3c; margin-top: 0;">⚡ 性能分析工具</h3>
   <p>使用性能分析工具识别瓶颈：</p>
   </div>

**Python性能分析**

.. code:: python

   import cProfile
   import pstats

   # 性能分析
   profiler = cProfile.Profile()
   profiler.enable()
   
   # 您的代码
   
   profiler.disable()
   stats = pstats.Stats(profiler)
   stats.sort_stats('cumulative')
   stats.print_stats(10)

**系统性能监控**

.. code:: bash

   # 监控CPU使用
   top -p $(pgrep python)
   
   # 监控内存使用
   ps aux | grep python
   
   # 监控磁盘I/O
   iostat -x 1

**内存分析**

.. code:: python

   import psutil
   import os

   # 获取内存使用
   process = psutil.Process(os.getpid())
   memory_info = process.memory_info()
   print(f"内存使用: {memory_info.rss / 1024 / 1024:.2f} MB")

获取帮助
--------

在线资源
~~~~~~~~

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #1976d2; margin-top: 0;">🌐 在线资源</h3>
   <p>获取更多帮助和支持：</p>
   </div>

**官方资源**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- **GitHub仓库**: https://github.com/xinyuren-bio/LNB-MDT
- **文档网站**: https://lnb-mdt.readthedocs.io
- **问题报告**: https://github.com/xinyuren-bio/LNB-MDT/issues
- **讨论区**: https://github.com/xinyuren-bio/LNB-MDT/discussions

   </div>

**社区支持**

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

- **邮件支持**: zy2310205@buaa.edu.cn
- **学术交流**: 相关学术会议和研讨会
- **用户群组**: 分子动力学研究社区

   </div>

**相关文档**

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

- **MDAnalysis文档**: https://www.mdanalysis.org/
- **VMD文档**: https://www.ks.uiuc.edu/Research/vmd/
- **Python文档**: https://docs.python.org/
- **Conda文档**: https://docs.conda.io/

   </div>

联系支持
~~~~~~~~

.. raw:: html

   <div style="background-color: #fce4ec; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #c2185b; margin-top: 0;">📞 联系支持</h3>
   <p>如果问题仍然存在，请联系技术支持：</p>
   </div>

**报告问题**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

**提供以下信息：**

1. **系统信息**
   - 操作系统版本
   - Python版本
   - LNB-MDT版本

2. **错误信息**
   - 完整的错误消息
   - 错误发生时的操作步骤
   - 相关的日志文件

3. **环境信息**
   - 硬件配置
   - 软件环境
   - 网络环境

4. **复现步骤**
   - 详细的操作步骤
   - 输入数据信息
   - 预期结果

   </div>

**技术支持邮箱**

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**技术支持**: zy2310205@buaa.edu.cn

**邮件主题格式**: [LNB-MDT] 问题描述

**邮件内容模板**:
.. code-block:: text

   主题: [LNB-MDT] 分析模块运行错误

   系统信息:
   - 操作系统: Windows 10
   - Python版本: 3.11.0
   - LNB-MDT版本: v1.0

   问题描述:
   运行PCA分析时出现内存不足错误

   错误信息:
   [粘贴完整错误信息]

   操作步骤:
   1. 启动LNB-MDT
   2. 加载轨迹文件
   3. 运行PCA分析
   4. 出现错误

   期望结果:
   正常完成PCA分析

   </div>

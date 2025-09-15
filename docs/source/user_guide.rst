用户指南
========

本指南详细介绍LNB-MDT的图形界面使用方法。

界面概览
--------

主界面布局
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #2c3e50; margin-top: 0;">🖥️ LNB-MDT主界面</h3>
   <p>LNB-MDT采用现代化的Qt6界面设计，提供直观的用户体验：</p>
   </div>

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; margin: 20px 0;">

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 10px;">
   <h4 style="margin-top: 0;">📋 左侧菜单</h4>
   <ul style="margin-bottom: 0;">
   <li>主页</li>
   <li>生成模块</li>
   <li>分析模块</li>
   <li>图表模块</li>
   <li>数据处理</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); color: white; padding: 20px; border-radius: 10px;">
   <h4 style="margin-top: 0;">🎛️ 中央工作区</h4>
   <ul style="margin-bottom: 0;">
   <li>参数配置</li>
   <li>文件选择</li>
   <li>进度显示</li>
   <li>结果展示</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); color: white; padding: 20px; border-radius: 10px;">
   <h4 style="margin-top: 0;">⚙️ 右侧面板</h4>
   <ul style="margin-bottom: 0;">
   <li>高级选项</li>
   <li>参数优化</li>
   <li>设置面板</li>
   <li>帮助信息</li>
   </ul>
   </div>

   </div>

功能模块详解
------------

主页模块
~~~~~~~~

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #1976d2; margin-top: 0;">🏠 主页功能</h3>
   <p>主页提供项目概览和快速导航：</p>
   
   <ul>
   <li><strong>项目信息</strong>: 显示LNB-MDT版本和基本信息</li>
   <li><strong>快速开始</strong>: 提供快速入门指导</li>
   <li><strong>最近项目</strong>: 显示最近使用的项目</li>
   <li><strong>系统状态</strong>: 显示系统资源使用情况</li>
   </ul>
   </div>

**主要功能：**

- 项目概览和状态显示
- 快速导航到各功能模块
- 系统资源监控
- 帮助和文档链接

生成模块
~~~~~~~~

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #7b1fa2; margin-top: 0;">🧬 脂质纳米泡生成</h3>
   <p>生成模块用于创建脂质纳米泡结构：</p>
   </div>

**功能特点：**

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #388e3c;">📐 几何参数</h4>
   <ul style="margin-bottom: 0;">
   <li>盒子尺寸 (X, Y, Z)</li>
   <li>脂质密度</li>
   <li>气体密度</li>
   <li>溶剂浓度</li>
   </ul>
   </div>

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #f57c00;">🧪 分子组成</h4>
   <ul style="margin-bottom: 0;">
   <li>脂质类型选择</li>
   <li>胆固醇比例</li>
   <li>添加剂配置</li>
   <li>离子浓度</li>
   </ul>
   </div>

   <div style="background-color: #fce4ec; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #c2185b;">⚙️ 生成选项</h4>
   <ul style="margin-bottom: 0;">
   <li>随机种子</li>
   <li>生成算法</li>
   <li>质量控制</li>
   <li>输出格式</li>
   </ul>
   </div>

   </div>

**使用步骤：**

1. **设置几何参数**
   - 输入盒子尺寸
   - 设置脂质和气体密度
   - 配置溶剂参数

2. **选择分子组成**
   - 选择脂质类型
   - 设置胆固醇比例
   - 配置添加剂

3. **配置生成选项**
   - 设置随机种子
   - 选择生成算法
   - 启用质量控制

4. **运行生成**
   - 点击"生成"按钮
   - 等待生成完成
   - 保存生成的结构

分析模块
~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #388e3c; margin-top: 0;">📊 分子动力学分析</h3>
   <p>分析模块提供多种分子动力学分析方法：</p>
   </div>

**分析类型选择**

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0;">📐 PCA分析</h4>
   <p style="margin-bottom: 0;">主成分分析</p>
   </div>

   <div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); color: white; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0;">📏 面积分析</h4>
   <p style="margin-bottom: 0;">Voronoi镶嵌</p>
   </div>

   <div style="background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); color: white; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0;">🌊 曲率分析</h4>
   <p style="margin-bottom: 0;">膜曲率计算</p>
   </div>

   <div style="background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%); color: white; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0;">📊 高度分析</h4>
   <p style="margin-bottom: 0;">高度分布</p>
   </div>

   <div style="background: linear-gradient(135deg, #fa709a 0%, #fee140 100%); color: white; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0;">🔗 聚类分析</h4>
   <p style="margin-bottom: 0;">分子聚集</p>
   </div>

   <div style="background: linear-gradient(135deg, #a8edea 0%, #fed6e3 100%); color: #333; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0;">🎯 各向异性</h4>
   <p style="margin-bottom: 0;">取向分析</p>
   </div>

   </div>

**分析流程：**

1. **加载数据文件**
   - 选择GRO拓扑文件
   - 选择XTC轨迹文件
   - 设置输出路径

2. **配置分析参数**
   - 选择残基组
   - 设置帧范围
   - 配置计算参数

3. **运行分析**
   - 点击"下一步"按钮
   - 选择分析类型
   - 启动分析过程

4. **查看结果**
   - 实时进度显示
   - 结果预览
   - 保存分析结果

**参数配置详解**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

**通用参数：**

- **起始帧**: 分析开始的时间帧
- **结束帧**: 分析结束的时间帧（-1表示到最后）
- **步长**: 帧之间的间隔
- **k值**: 局部邻域大小
- **残基组**: 要分析的分子类型和原子

**特定参数：**
- **面积分析**: 最大法线角度
- **曲率分析**: 曲率类型（平均/高斯）
- **聚类分析**: 截止距离
- **高度分析**: 参考原子组

   </div>

图表模块
~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #f57c00; margin-top: 0;">📈 数据可视化</h3>
   <p>图表模块提供丰富的数据可视化功能：</p>
   </div>

**图表类型**

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background-color: #e3f2fd; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #1976d2;">📊 线图</h4>
   <ul style="margin-bottom: 0;">
   <li>时间序列图</li>
   <li>趋势分析</li>
   <li>多线对比</li>
   </ul>
   </div>

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #7b1fa2;">📊 柱状图</h4>
   <ul style="margin-bottom: 0;">
   <li>分布直方图</li>
   <li>对比柱状图</li>
   <li>堆叠柱状图</li>
   </ul>
   </div>

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #388e3c;">📊 散点图</h4>
   <ul style="margin-bottom: 0;">
   <li>相关性分析</li>
   <li>聚类可视化</li>
   <li>异常检测</li>
   </ul>
   </div>

   </div>

**图表功能**

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**主要功能：**

- **数据导入**: 支持CSV、Excel格式
- **图表定制**: 颜色、样式、标签自定义
- **交互操作**: 缩放、平移、选择
- **导出功能**: 高质量图片和PDF导出
- **统计分析**: 内置统计计算

   </div>

**使用步骤：**

1. **加载数据**
   - 选择数据文件
   - 预览数据内容
   - 选择数据列

2. **选择图表类型**
   - 线图：时间序列数据
   - 柱状图：分类数据
   - 散点图：相关性数据

3. **自定义样式**
   - 设置颜色方案
   - 调整图表大小
   - 添加标题和标签

4. **生成图表**
   - 预览图表效果
   - 保存图表文件
   - 导出高质量图片

数据处理模块
~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #fce4ec; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #c2185b; margin-top: 0;">🔧 VMD集成和数据处理</h3>
   <p>数据处理模块提供VMD集成和高级数据处理功能：</p>
   </div>

**VMD集成功能**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

**主要功能：**

- **VMD启动**: 一键启动VMD程序
- **数据传递**: 自动传递分析结果到VMD
- **可视化控制**: 控制VMD的显示效果
- **交互操作**: 在LNB-MDT中选择分子，VMD同步显示

   </div>

**VMD操作流程**

1. **启动VMD**
   - 点击"Start VMD"按钮
   - 等待VMD程序启动
   - 确认连接状态

2. **加载数据**
   - 拖拽CSV文件到VMD窗口
   - 或使用文件菜单加载
   - 查看数据表格

3. **可视化操作**
   - 选择要显示的帧
   - 选择要高亮的分子
   - VMD自动跳转和高亮

4. **停止VMD**
   - 点击"Stop VMD"按钮
   - 关闭VMD程序

**数据处理功能**

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #388e3c;">🧹 数据清洗</h4>
   <ul style="margin-bottom: 0;">
   <li>去除异常值</li>
   <li>填充缺失值</li>
   <li>数据验证</li>
   </ul>
   </div>

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #f57c00;">🔄 格式转换</h4>
   <ul style="margin-bottom: 0;">
   <li>CSV转Excel</li>
   <li>数据格式标准化</li>
   <li>编码转换</li>
   </ul>
   </div>

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #7b1fa2;">⚡ 批量处理</h4>
   <ul style="margin-bottom: 0;">
   <li>批量文件处理</li>
   <li>自动化分析</li>
   <li>结果汇总</li>
   </ul>
   </div>

   </div>

高级功能
--------

参数优化
~~~~~~~~

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #1976d2; margin-top: 0;">🤖 智能参数优化</h3>
   <p>使用机器学习技术自动优化分析参数：</p>
   </div>

**优化功能**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

**支持的优化：**

- **k值优化**: 自动寻找最佳k值
- **截止距离优化**: 优化聚类参数
- **帧范围优化**: 选择最佳分析范围
- **并行参数优化**: 优化并行处理参数

   </div>

**使用步骤：**

1. **选择优化类型**
   - 在右侧面板选择优化选项
   - 设置优化目标
   - 配置优化参数

2. **运行优化**
   - 点击"开始优化"按钮
   - 等待优化完成
   - 查看优化结果

3. **应用优化结果**
   - 自动应用最佳参数
   - 重新运行分析
   - 比较优化效果

并行处理
~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #f57c00; margin-top: 0;">⚡ 并行处理加速</h3>
   <p>利用多核CPU加速分析计算：</p>
   </div>

**并行选项**

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**并行设置：**

- **启用并行**: 勾选并行处理选项
- **核心数量**: 设置使用的CPU核心数
- **内存管理**: 自动管理内存使用
- **进度监控**: 实时显示并行进度

   </div>

**性能优化建议**

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**优化建议：**

- **CPU核心**: 使用所有可用核心（-1）
- **内存考虑**: 大系统减少并行数
- **I/O优化**: 使用SSD存储轨迹文件
- **网络优化**: 避免网络文件系统

   </div>

主题和个性化
~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #fce4ec; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #c2185b; margin-top: 0;">🎨 界面主题和个性化</h3>
   <p>LNB-MDT支持多种界面主题：</p>
   </div>

**主题选择**

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background-color: #2c3e50; color: white; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0;">🌙 深色主题</h4>
   <p style="margin-bottom: 0;">适合长时间使用</p>
   </div>

   <div style="background-color: #ecf0f1; color: #2c3e50; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0;">☀️ 浅色主题</h4>
   <p style="margin-bottom: 0;">适合明亮环境</p>
   </div>

   <div style="background-color: #8e44ad; color: white; padding: 15px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0;">🎨 自定义主题</h4>
   <p style="margin-bottom: 0;">用户自定义颜色</p>
   </div>

   </div>

**个性化设置**

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**可自定义项目：**

- **界面主题**: 深色/浅色/自定义
- **字体大小**: 调整界面字体
- **语言设置**: 中文/英文界面
- **快捷键**: 自定义快捷键
- **默认路径**: 设置默认文件路径

   </div>

快捷键参考
----------

常用快捷键
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

**全局快捷键：**

- **Ctrl+N**: 新建项目
- **Ctrl+O**: 打开项目
- **Ctrl+S**: 保存项目
- **Ctrl+Q**: 退出程序
- **F1**: 显示帮助

**模块快捷键：**

- **Ctrl+1**: 切换到主页
- **Ctrl+2**: 切换到生成模块
- **Ctrl+3**: 切换到分析模块
- **Ctrl+4**: 切换到图表模块
- **Ctrl+5**: 切换到数据处理模块

   </div>

故障排除
--------

常见问题解决
~~~~~~~~~~~~

问题1：界面无法启动
^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #ffebee; padding: 15px; border-radius: 8px; border-left: 4px solid #f44336;">

**可能原因和解决方案：**

1. **Python环境问题**
   - 检查Python版本 >= 3.11
   - 确认PySide6已安装
   - 重新安装依赖包

2. **显示问题**
   - 检查显卡驱动
   - 尝试软件渲染模式
   - 调整显示缩放

   </div>

问题2：文件加载失败
^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**解决方案：**

1. **文件格式检查**
   - 确认GRO/XTC文件格式正确
   - 检查文件是否损坏
   - 验证文件路径

2. **权限问题**
   - 检查文件读取权限
   - 以管理员身份运行
   - 检查防火墙设置

   </div>

问题3：分析结果异常
^^^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**检查步骤：**

1. **参数验证**
   - 检查残基组格式
   - 验证帧范围设置
   - 确认k值合理性

2. **数据质量**
   - 检查轨迹质量
   - 验证拓扑文件
   - 确认时间步长

   </div>

问题4：VMD连接失败
^^^^^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**解决方案：**

1. **VMD安装检查**
   - 确认VMD已正确安装
   - 检查VMD路径设置
   - 验证VMD版本

2. **网络连接**
   - 检查防火墙设置
   - 确认端口未被占用
   - 尝试重启VMD

   </div>

性能优化建议
~~~~~~~~~~~~

界面响应优化
^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**优化建议：**

- **硬件升级**: 使用SSD和充足内存
- **系统优化**: 关闭不必要的后台程序
- **界面设置**: 降低动画效果
- **数据管理**: 定期清理临时文件

   </div>

内存使用优化
^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #fce4ec; padding: 15px; border-radius: 8px; border-left: 4px solid #e91e63;">

**内存管理：**

- **分段处理**: 大文件分段加载
- **缓存清理**: 定期清理内存缓存
- **并行控制**: 合理设置并行数量
- **数据压缩**: 使用压缩格式存储

   </div>

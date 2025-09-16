贡献指南
==========

感谢您对LNB-MDT项目的关注！我们欢迎各种形式的贡献。

如何贡献
--------

贡献方式
~~~~~~~~

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; margin: 20px 0;">

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 10px;">
   <h4 style="margin-top: 0;">🐛 报告问题</h4>
   <p>发现bug或问题？请报告给我们</p>
   </div>

   <div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); color: white; padding: 20px; border-radius: 10px;">
   <h4 style="margin-top: 0;">💡 功能建议</h4>
   <p>有好的想法？欢迎提出新功能建议</p>
   </div>

   <div style="background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); color: white; padding: 20px; border-radius: 10px;">
   <h4 style="margin-top: 0;">📝 改进文档</h4>
   <p>帮助改进文档和教程</p>
   </div>

   <div style="background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%); color: white; padding: 20px; border-radius: 10px;">
   <h4 style="margin-top: 0;">🔧 代码贡献</h4>
   <p>提交代码改进和新功能</p>
   </div>

   </div>

报告问题
~~~~~~~~

.. raw:: html

   <div style="background-color: #ffebee; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #d32f2f; margin-top: 0;">🐛 如何报告问题</h3>
   <p>在报告问题之前，请先检查是否已有类似问题：</p>
   </div>

**检查现有问题**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

1. 访问 https://github.com/xinyuren-bio/LNB-MDT/issues
2. 搜索相关关键词
3. 查看是否有类似问题
4. 如果存在，请在现有问题下评论

   </div>

**创建新问题**

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**问题报告模板：**

.. code:: text

   **问题描述**
   简要描述遇到的问题

   **复现步骤**
   1. 执行的操作
   2. 期望的结果
   3. 实际的结果

   **环境信息**
   - 操作系统: 
   - Python版本: 
   - LNB-MDT版本: 
   - 其他相关信息

   **错误信息**
   [粘贴完整的错误信息]

   **附加信息**
   任何其他相关信息

   </div>

功能建议
~~~~~~~~

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #1976d2; margin-top: 0;">💡 如何提出功能建议</h3>
   <p>我们欢迎新功能建议，请使用以下格式：</p>
   </div>

**功能建议模板**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

.. code:: text

   **功能描述**
   简要描述建议的功能

   **使用场景**
   描述该功能的使用场景和好处

   **实现建议**
   如果有实现想法，请提供建议

   **优先级**
   高/中/低

   **附加信息**
   任何其他相关信息

   </div>

代码贡献
--------

开发环境设置
~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #388e3c; margin-top: 0;">🔧 设置开发环境</h3>
   <p>按照以下步骤设置开发环境：</p>
   </div>

**1. Fork仓库**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

1. 访问 https://github.com/xinyuren-bio/LNB-MDT
2. 点击右上角的 "Fork" 按钮
3. 克隆您的fork到本地

.. code:: bash

   git clone https://github.com/your-username/LNB-MDT.git
   cd LNB-MDT

   </div>

**2. 创建开发分支**

.. code:: bash

   # 创建并切换到新分支
   git checkout -b feature/your-feature-name
   
   # 或者修复bug
   git checkout -b bugfix/issue-number

**3. 安装开发依赖**

.. code:: bash

   # 创建开发环境
   conda create -n LNB-MDT-dev python=3.11
   conda activate LNB-MDT-dev
   
   # 安装依赖
   pip install -r requirements.txt
   pip install -r requirements-dev.txt  # 开发依赖

**4. 安装预提交钩子**

.. code:: bash

   # 安装pre-commit
   pip install pre-commit
   pre-commit install

代码规范
~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #f57c00; margin-top: 0;">📋 代码规范</h3>
   <p>请遵循以下代码规范：</p>
   </div>

**Python代码规范**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- **PEP 8**: 遵循Python官方代码规范
- **类型提示**: 使用类型提示提高代码可读性
- **文档字符串**: 为所有函数和类添加文档字符串
- **命名规范**: 使用清晰的变量和函数名

   </div>

**代码格式**

.. code:: bash

   # 使用black格式化代码
   black your_file.py
   
   # 使用isort排序导入
   isort your_file.py
   
   # 使用flake8检查代码
   flake8 your_file.py

**文档字符串格式**

.. code:: python

   def example_function(param1: str, param2: int) -> bool:
       """
       示例函数的文档字符串
       
       参数:
       - param1 (str): 第一个参数
       - param2 (int): 第二个参数
       
       返回:
       - bool: 返回值说明
       
       异常:
       - ValueError: 当参数无效时抛出
       """
       pass

测试规范
~~~~~~~~

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #7b1fa2; margin-top: 0;">🧪 测试规范</h3>
   <p>所有代码贡献都应包含相应的测试：</p>
   </div>

**测试结构**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

.. code:: text

   tests/
   ├── test_analysis/          # 分析模块测试
   │   ├── test_pca.py
   │   ├── test_area.py
   │   └── test_curvature.py
   ├── test_machine_learning/  # 机器学习测试
   │   ├── test_optimizer.py
   │   └── test_detector.py
   ├── test_modules/           # 模块测试
   │   └── test_vmd_control.py
   └── test_utils/             # 工具函数测试
       └── test_file_utils.py

   </div>

**测试示例**

.. code:: python

   import unittest
   import numpy as np
   from analysis.pca import PCA

   class TestPCA(unittest.TestCase):
       """PCA分析测试类"""
       
       def setUp(self):
           """测试前准备"""
           self.analyzer = PCA(
               gro_file="test_data/test.gro",
               xtc_file="test_data/test.xtc",
               residues={'DPPC': ['PO4']}
           )
       
       def test_initialization(self):
           """测试初始化"""
           self.assertIsNotNone(self.analyzer)
           self.assertEqual(self.analyzer.n_components, 3)
       
       def test_run_analysis(self):
           """测试分析运行"""
           results = self.analyzer.run(start_frame=0, stop_frame=10)
           self.assertIsInstance(results, dict)
           self.assertIn('frames', results)
           self.assertIn('values', results)

   if __name__ == '__main__':
       unittest.main()

**运行测试**

.. code:: bash

   # 运行所有测试
   python -m pytest tests/
   
   # 运行特定测试
   python -m pytest tests/test_analysis/test_pca.py
   
   # 运行测试并生成覆盖率报告
   python -m pytest --cov=analysis tests/

提交代码
~~~~~~~~

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #03a9f4; margin-top: 0;">📤 提交代码</h3>
   <p>按照以下步骤提交您的代码：</p>
   </div>

**1. 提交更改**

.. code:: bash

   # 添加更改的文件
   git add .
   
   # 提交更改
   git commit -m "feat: 添加新功能描述"
   
   # 推送到您的fork
   git push origin feature/your-feature-name

**2. 创建Pull Request**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

1. 访问您的GitHub fork页面
2. 点击 "New Pull Request" 按钮
3. 选择您的分支
4. 填写PR描述
5. 提交PR

   </div>

**PR描述模板**

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

.. code:: text

   **变更描述**
   简要描述此PR的变更内容

   **变更类型**
   - [ ] Bug修复
   - [ ] 新功能
   - [ ] 文档更新
   - [ ] 性能优化
   - [ ] 重构

   **测试**
   - [ ] 添加了测试
   - [ ] 所有测试通过
   - [ ] 手动测试完成

   **检查清单**
   - [ ] 代码遵循项目规范
   - [ ] 文档已更新
   - [ ] 没有破坏性变更

   **相关Issue**
   关联的Issue编号: #123

   </div>

文档贡献
--------

文档类型
~~~~~~~~

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background-color: #e3f2fd; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #1976d2;">📖 用户文档</h4>
   <ul style="margin-bottom: 0;">
   <li>安装指南</li>
   <li>使用教程</li>
   <li>示例代码</li>
   </ul>
   </div>

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #388e3c;">🔧 开发者文档</h4>
   <ul style="margin-bottom: 0;">
   <li>API参考</li>
   <li>架构说明</li>
   <li>贡献指南</li>
   </ul>
   </div>

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #f57c00;">📝 教程文档</h4>
   <ul style="margin-bottom: 0;">
   <li>快速开始</li>
   <li>最佳实践</li>
   <li>故障排除</li>
   </ul>
   </div>

   </div>

文档格式
~~~~~~~~

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #2c3e50; margin-top: 0;">📝 文档格式要求</h3>
   <p>文档使用reStructuredText格式，请遵循以下规范：</p>
   </div>

**reStructuredText语法**

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

.. code:: rst

   ===========
   标题
   ===========

   这是段落文本。

   **粗体文本**

   *斜体文本*

   .. code:: python

       # 代码示例
       def example():
           pass

   - 列表项1
   - 列表项2

   .. note::

       这是一个注释。

   </div>

**文档结构**

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

.. code:: text

   docs/source/
   ├── index.rst              # 主页
   ├── installation.rst       # 安装指南
   ├── quickstart.rst         # 快速开始
   ├── user_guide.rst         # 用户指南
   ├── analysis_modules.rst   # 分析模块
   ├── machine_learning.rst   # 机器学习
   ├── command_line.rst       # 命令行工具
   ├── examples.rst           # 使用示例
   ├── api_reference.rst      # API参考
   ├── troubleshooting.rst    # 故障排除
   └── contributing.rst       # 贡献指南

   </div>

**本地构建文档**

.. code:: bash

   # 安装文档依赖
   pip install sphinx sphinx-rtd-theme myst-parser

   # 构建文档
   cd docs
   make html

   # 查看文档
   open build/html/index.html

社区参与
--------

讨论和反馈
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #1976d2; margin-top: 0;">💬 参与讨论</h3>
   <p>我们欢迎社区成员参与讨论：</p>
   </div>

**讨论渠道**

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #2c3e50;">📧 邮件列表</h4>
   <p style="margin-bottom: 0;">zy2310205@buaa.edu.cn</p>
   </div>

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #2c3e50;">💬 GitHub讨论</h4>
   <p style="margin-bottom: 0;">项目讨论区</p>
   </div>

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #2c3e50;">🏛️ 学术会议</h4>
   <p style="margin-bottom: 0;">相关学术会议</p>
   </div>

   </div>

**反馈类型**

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

- **Bug报告**: 详细描述问题和复现步骤
- **功能请求**: 描述新功能的需求和使用场景
- **文档改进**: 指出文档中的错误或改进建议
- **使用经验**: 分享使用LNB-MDT的经验和技巧

   </div>

代码审查
~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #f57c00; margin-top: 0;">👀 代码审查流程</h3>
   <p>所有代码贡献都会经过审查：</p>
   </div>

**审查标准**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

1. **代码质量**: 代码是否清晰、可读
2. **功能正确性**: 功能是否按预期工作
3. **测试覆盖**: 是否有足够的测试
4. **文档更新**: 相关文档是否已更新
5. **性能影响**: 是否影响系统性能

   </div>

**审查流程**

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

1. **自动检查**: CI/CD系统自动运行测试
2. **人工审查**: 维护者审查代码
3. **反馈修改**: 根据反馈进行修改
4. **最终审查**: 维护者最终确认
5. **合并代码**: 代码合并到主分支

   </div>

发布流程
~~~~~~~~

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #03a9f4; margin-top: 0;">🚀 发布流程</h3>
   <p>了解LNB-MDT的发布流程：</p>
   </div>

**版本管理**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- **语义化版本**: 使用SemVer版本号格式
- **主版本**: 重大变更或不兼容变更
- **次版本**: 新功能添加
- **修订版本**: Bug修复

   </div>

**发布周期**

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

- **主要版本**: 每年1-2次
- **次要版本**: 每季度1次
- **修订版本**: 根据需要发布
- **预发布版本**: 用于测试

   </div>

**发布检查清单**

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

- [ ] 所有测试通过
- [ ] 文档已更新
- [ ] 版本号已更新
- [ ] 变更日志已更新
- [ ] 发布说明已准备

   </div>

贡献者认可
----------

贡献者名单
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #1976d2; margin-top: 0;">👥 贡献者认可</h3>
   <p>我们感谢所有贡献者的努力：</p>
   </div>

**贡献类型**

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #2c3e50;">💻 代码贡献</h4>
   <ul style="margin-bottom: 0;">
   <li>功能开发</li>
   <li>Bug修复</li>
   <li>性能优化</li>
   </ul>
   </div>

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #2c3e50;">📝 文档贡献</h4>
   <ul style="margin-bottom: 0;">
   <li>用户指南</li>
   <li>API文档</li>
   <li>教程编写</li>
   </ul>
   </div>

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #2c3e50;">🧪 测试贡献</h4>
   <ul style="margin-bottom: 0;">
   <li>单元测试</li>
   <li>集成测试</li>
   <li>性能测试</li>
   </ul>
   </div>

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #2c3e50;">💬 社区贡献</h4>
   <ul style="margin-bottom: 0;">
   <li>问题报告</li>
   <li>功能建议</li>
   <li>用户支持</li>
   </ul>
   </div>

   </div>

**认可方式**

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

- **贡献者名单**: 在README中列出贡献者
- **GitHub贡献**: 在GitHub上显示贡献记录
- **学术认可**: 在相关论文中致谢
- **社区认可**: 在社区中公开感谢

   </div>

行为准则
~~~~~~~~

.. raw:: html

   <div style="background-color: #fce4ec; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #c2185b; margin-top: 0;">🤝 行为准则</h3>
   <p>我们致力于为每个人提供友好、安全的环境：</p>
   </div>

**我们的承诺**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- **包容性**: 欢迎所有背景的人参与
- **尊重**: 尊重不同的观点和经验
- **合作**: 促进建设性的合作
- **安全**: 提供安全的环境

   </div>

**不可接受的行为**

.. raw:: html

   <div style="background-color: #ffebee; padding: 15px; border-radius: 8px; border-left: 4px solid #f44336;">

- 使用性暗示的语言或图像
- 恶意评论、侮辱或人身攻击
- 公开或私下骚扰
- 未经许可发布他人私人信息
- 其他不专业的行为

   </div>

**报告问题**

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

如果您遇到违反行为准则的情况，请通过以下方式报告：

- **邮件**: zy2310205@buaa.edu.cn
- **GitHub**: 通过GitHub的举报功能

   </div>

许可证
------

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #388e3c; margin-top: 0;">📄 开源许可证</h3>
   <p>LNB-MDT使用MIT许可证：</p>
   </div>

**MIT许可证**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

MIT许可证是一个宽松的开源许可证，允许：

- ✅ 商业使用
- ✅ 修改
- ✅ 分发
- ✅ 私人使用

**要求**:
- 📋 包含许可证和版权声明

   </div>

**许可证文本**

.. code:: text

   MIT License

   Copyright (c) 2025 LNB-MDT Contributors

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

联系我们
--------

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #1976d2; margin-top: 0;">📞 联系我们</h3>
   <p>如果您有任何问题或建议，请随时联系我们：</p>
   </div>

**联系方式**

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #2c3e50;">📧 邮件联系</h4>
   <p style="margin-bottom: 0;">zy2310205@buaa.edu.cn</p>
   </div>

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #2c3e50;">🐙 GitHub</h4>
   <p style="margin-bottom: 0;">@xinyuren-bio</p>
   </div>

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px;">
   <h4 style="margin-top: 0; color: #2c3e50;">🏛️ 学术机构</h4>
   <p style="margin-bottom: 0;">北京航空航天大学</p>
   </div>

   </div>

**项目信息**

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

- **项目名称**: LNB-MDT (Lipid NanoBubble Molecular Dynamics Toolbox)
- **版本**: v1.0
- **许可证**: MIT License
- **开发语言**: Python
- **主要作者**: XinyuRen

   </div>

感谢您的贡献！🎉

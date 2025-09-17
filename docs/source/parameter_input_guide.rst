参数输入指南
============

LNB-MDT提供了多种灵活的参数输入方式，让您可以根据需要选择最适合的方法。

概述
----

LNB-MDT支持以下参数输入方式：

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; margin: 20px 0;">

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">⚡ 短参数别名</h3>
   <p>所有参数都有简短的别名，让命令行更简洁</p>
   <ul style="margin-bottom: 0;">
   <li>-g 代替 --gro-file</li>
   <li>-r 代替 --residues</li>
   <li>-p 代替 --parallel</li>
   </ul>
   </div>

   <div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">📝 简化格式</h3>
   <p>支持更直观的参数格式</p>
   <ul style="margin-bottom: 0;">
   <li>DPPC:PO4,CHOL:ROH</li>
   <li>N2:N2,O2:O2</li>
   <li>DPPC:PO4+GLY</li>
   </ul>
   </div>


   <div style="background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">🔄 向后兼容</h3>
   <p>传统格式仍然完全支持</p>
   <ul style="margin-bottom: 0;">
   <li>字典字符串</li>
   <li>长参数名</li>
   <li>复杂格式</li>
   </ul>
   </div>

   </div>

短参数别名
----------

所有命令行参数都有对应的短别名，让命令行更加简洁：

.. list-table:: 参数别名对照表
   :header-rows: 1
   :widths: 10 20 30 40

   * - 短参数
     - 长参数
     - 说明
     - 示例
   * - ``-g``
     - ``--gro-file``
     - GRO文件路径
     - ``-g cases/lnb.gro``
   * - ``-x``
     - ``--xtc-file``
     - XTC文件路径
     - ``-x cases/md.xtc``
   * - ``-o``
     - ``--output-csv``
     - 输出CSV文件路径
     - ``-o results.csv``
   * - ``-r``
     - ``--residues``
     - 残基组定义
     - ``-r DPPC:PO4``
   * - ``-a``
     - ``--gas-group``
     - 气体组定义
     - ``-a N2:N2``
   * - ``-m``
     - ``--MW``
     - 分子量
     - ``-m 14``
   * - ``-R``
     - ``--radius``
     - 半径
     - ``-R 50``
   * - ``-p``
     - ``--parallel``
     - 启用并行处理
     - ``-p``
   * - ``-j``
     - ``--n-jobs``
     - 并行任务数
     - ``-j 4``
   * - ``-s``
     - ``--start-frame``
     - 起始帧
     - ``-s 0``
   * - ``-e``
     - ``--stop-frame``
     - 结束帧
     - ``-e 100``
   * - ``-t``
     - ``--step-frame``
     - 帧步长
     - ``-t 5``
   * - ``-v``
     - ``--verbose``
     - 详细输出
     - ``-v``
   * - ``-k``
     - ``--k-value``
     - k值
     - ``-k 20``
   * - ``-M``
     - ``--method``
     - 计算方法
     - ``-M mean``
   * - ``-T``
     - ``--threshold``
     - 阈值
     - ``-T 0.5``
   * - ``-P``
     - ``--plot-type``
     - 图表类型
     - ``-P line``
   * - ``-d``
     - ``--plot-dir``
     - 图表目录
     - ``-d plots/``

简化格式
--------

residues和gas-group参数现在支持更直观的输入格式：

基本格式
~~~~~~~~

**简单格式（推荐）:**
.. code:: bash

   # 基本格式: RESIDUE:ATOM
   -r DPPC:PO4,CHOL:ROH
   -a N2:N2
   
   # 多个残基/气体
   -r DPPC:PO4,DUPC:PO4,CHOL:ROH
   -a N2:N2,O2:O2

**多原子格式:**
.. code:: bash

   # 多原子: RESIDUE:ATOM1+ATOM2
   -r DPPC:PO4+GLY,CHOL:ROH
   -r DPPC:PO4+GLY+CH2,CHOL:ROH

**只有名称格式:**
.. code:: bash

   # 只有残基/气体名（原子名与名称相同）
   -r DPPC
   -a N2


传统格式
~~~~~~~~

**字典字符串格式（仍然支持）:**
.. code:: bash

   # 传统字典格式
   -r "{'DPPC': ['PO4'], 'CHOL': ['ROH']}"
   -a "{'N2': ['N2']}"

使用示例对比
------------

密度分析示例
~~~~~~~~~~~~

**传统方式:**
.. code:: bash

   python analysis/densitywithframe.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --residues "{'DPPC': ['PO4'], 'DUPC': ['PO4'], 'CHOL': ['ROH']}" \
     --gas-group "{'N2': ['N2']}" \
     --MW 14 \
     --radius 50 \
     --output-csv results.csv \
     --parallel \
     --n-jobs 4

**简化方式:**
.. code:: bash

   python analysis/densitywithframe.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -r DPPC:PO4,DUPC:PO4,CHOL:ROH \
     -a N2:N2 \
     -m 14 \
     -R 50 \
     -o results.csv \
     -p \
     -j 4


PCA分析示例
~~~~~~~~~~~

**传统方式:**
.. code:: bash

   python analysis/pca.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" \
     --start-frame 0 \
     --stop-frame 100 \
     --parallel \
     --verbose

**简化方式:**
.. code:: bash

   python analysis/pca.py \
     -g cases/lnb.gro \
     -x cases/md.xtc \
     -r DPPC:PO4,CHOL:ROH \
     -s 0 \
     -e 100 \
     -p \
     -v


Python API使用
--------------

在Python代码中使用简化参数解析：

.. code:: python

   from analysis.parameter_utils import parse_residues_simple, parse_gas_group_simple

   # 简单格式
   residues = parse_residues_simple('DPPC:PO4,CHOL:ROH')
   gas_group = parse_gas_group_simple('N2:N2')

   # 多原子格式
   residues = parse_residues_simple('DPPC:PO4+GLY,CHOL:ROH')


   # 传统格式（仍然支持）
   residues = parse_residues_simple("{'DPPC': ['PO4'], 'CHOL': ['ROH']}")

支持的模块
----------

所有analysis模块都支持简化的参数输入：

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0;">

   <div style="background-color: #f3e5f5; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #7b1fa2;">densitywithframe.py</h4>
   </div>

   <div style="background-color: #e8f5e8; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #388e3c;">densitywithradius.py</h4>
   </div>

   <div style="background-color: #fff3e0; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #f57c00;">area.py</h4>
   </div>

   <div style="background-color: #fce4ec; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #c2185b;">height.py</h4>
   </div>

   <div style="background-color: #e3f2fd; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #1976d2;">curvature.py</h4>
   </div>

   <div style="background-color: #f1f8e9; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #689f38;">pca.py</h4>
   </div>

   <div style="background-color: #fef7e0; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #f9a825;">cluster.py</h4>
   </div>

   <div style="background-color: #f3e5f5; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #7b1fa2;">anisotropy.py</h4>
   </div>

   <div style="background-color: #e8f5e8; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #388e3c;">gyration.py</h4>
   </div>

   <div style="background-color: #fff3e0; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #f57c00;">sz.py</h4>
   </div>

   <div style="background-color: #fce4ec; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #c2185b;">n_cluster.py</h4>
   </div>

   <div style="background-color: #e3f2fd; padding: 10px; border-radius: 8px; text-align: center;">
   <h4 style="margin-top: 0; color: #1976d2;">rad.py</h4>
   </div>

   </div>

优势总结
--------

使用简化的参数输入方式有以下优势：

.. raw:: html

   <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;">

   <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">⚡ 更快速</h3>
   <p>短参数别名让命令行更简洁，输入更快</p>
   </div>

   <div style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">📝 更直观</h3>
   <p>简单格式更接近自然语言，易于理解</p>
   </div>


   <div style="background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%); color: white; padding: 20px; border-radius: 10px;">
   <h3 style="margin-top: 0;">🔄 兼容性</h3>
   <p>向后兼容，传统格式仍然支持</p>
   </div>

   </div>

注意事项
--------

1. **空格处理**: 参数中的空格会被自动处理
2. **大小写敏感**: 残基名和原子名区分大小写
3. **错误处理**: 如果格式不正确，会显示详细的错误信息和格式说明
4. **向后兼容**: 所有传统格式仍然完全支持

现在您可以享受更简单、更直观的命令行体验了！

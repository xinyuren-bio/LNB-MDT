命令行工具
============

LNB-MDT提供了完整的命令行界面，支持批量处理和自动化分析。

命令行概览
----------

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #2c3e50; margin-top: 0;">💻 命令行工具优势</h3>
   <p>命令行工具提供以下优势：</p>
   
   <ul>
   <li><strong>批量处理</strong>: 处理大量文件和数据</li>
   <li><strong>自动化</strong>: 编写脚本实现自动化分析</li>
   <li><strong>高性能</strong>: 更高效的资源利用</li>
   <li><strong>远程执行</strong>: 支持远程服务器运行</li>
   <li><strong>集成</strong>: 易于集成到工作流程中</li>
   </ul>
   </div>

通用参数
--------

所有分析模块都支持以下通用参数：

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**必需参数：**

- `--gro-file`: GRO拓扑文件路径
- `--xtc-file`: XTC轨迹文件路径
- `--output-csv`: 输出CSV文件路径
- `--residues`: 残基组字典字符串

**可选参数：**

- `--parallel`: 启用并行处理
- `--n-jobs`: 并行作业数量（-1表示使用所有CPU核心）
- `--start-frame`: 分析起始帧（0索引）
- `--stop-frame`: 分析结束帧（独占）
- `--step-frame`: 帧步长
- `--verbose`: 启用详细输出

   </div>

参数详解
~~~~~~~~

残基组参数格式
^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**residues参数格式：**

.. code:: bash

   # 基本格式
   --residues "{'DPPC': ['PO4']}"
   
   # 多分子类型
   --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}"
   
   # 多原子组（高度分析）
   --residues "{'DPPC': (['PO4'], ['C4B', 'C4A'])}"
   
   # 复杂组合
   --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH'], 'DUPC': ['PO4']}"

   </div>

**注意事项：**
- 必须使用双引号包围整个字典
- 字典键使用单引号
- 列表使用方括号
- 元组使用圆括号

并行处理参数
^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 15px; border-radius: 8px; border-left: 4px solid #2196f3;">

**并行处理选项：**

.. code:: bash

   # 启用并行处理
   --parallel
   
   # 指定并行核数
   --n-jobs 4
   
   # 使用所有可用核心
   --n-jobs -1
   
   # 禁用并行处理（默认）
   # 不添加 --parallel 参数

   </div>

**性能建议：**
- 小系统：使用2-4个核心
- 中等系统：使用4-8个核心
- 大系统：使用8-16个核心
- 内存不足时减少并行数

帧范围参数
^^^^^^^^^^

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**帧范围设置：**

.. code:: bash

   # 分析所有帧
   # 不指定 start-frame 和 stop-frame
   
   # 分析前100帧
   --start-frame 0 --stop-frame 100
   
   # 分析100-200帧
   --start-frame 100 --stop-frame 200
   
   # 分析最后100帧
   --stop-frame -1 --start-frame -100
   
   # 每10帧分析一次
   --step-frame 10

   </div>

分析模块详解
------------

PCA分析 (pca.py)
~~~~~~~~~~~~~~~~

**功能描述**
主成分分析，用于研究分子构象变化和运动模式。

**特定参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- `--n-components`: 主成分数量（默认：3）

   </div>

**使用示例**

.. code:: bash

   # 基本PCA分析
   python analysis/pca.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/pca_results.csv \
     --residues "{'DPPC': ['PO4']}" \
     --verbose

   # 高级PCA分析
   python analysis/pca.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/pca_advanced.csv \
     --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" \
     --n-components 5 \
     --start-frame 100 \
     --stop-frame 1000 \
     --step-frame 5 \
     --parallel \
     --n-jobs 4 \
     --verbose

面积分析 (area.py)
~~~~~~~~~~~~~~~~~~

**功能描述**
使用Voronoi镶嵌方法计算脂质分子的面积分布。

**特定参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- `--k-value`: Voronoi镶嵌的k值（默认：20）
- `--max-normal-angle`: 最大法线角度（默认：140度）

   </div>

**使用示例**

.. code:: bash

   # 基本面积分析
   python analysis/area.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/area_results.csv \
     --residues "{'DPPC': ['PO4']}" \
     --verbose

   # 高级面积分析
   python analysis/area.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/area_advanced.csv \
     --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" \
     --k-value 25 \
     --max-normal-angle 135 \
     --start-frame 0 \
     --stop-frame 500 \
     --parallel \
     --verbose

曲率分析 (curvature.py)
~~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
计算脂质膜的平均曲率和高斯曲率。

**特定参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- `--method`: 曲率类型（'mean' 或 'gaussian'）
- `--k-value`: 曲率计算的k值（默认：20）

   </div>

**使用示例**

.. code:: bash

   # 平均曲率分析
   python analysis/curvature.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/mean_curvature.csv \
     --residues "{'DPPC': ['PO4']}" \
     --method mean \
     --k-value 20 \
     --verbose

   # 高斯曲率分析
   python analysis/curvature.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/gaussian_curvature.csv \
     --residues "{'DPPC': ['PO4']}" \
     --method gaussian \
     --k-value 20 \
     --parallel \
     --verbose

高度分析 (height.py)
~~~~~~~~~~~~~~~~~~~~

**功能描述**
分析脂质分子的高度分布和膜厚度。

**特定参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- `--k-value`: 高度计算的k值（默认：20）

   </div>

**使用示例**

.. code:: bash

   # 基本高度分析
   python analysis/height.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/height_results.csv \
     --residues "{'DPPC': (['PO4'], ['C4B', 'C4A'])}" \
     --verbose

   # 多分子高度分析
   python analysis/height.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/height_multi.csv \
     --residues "{'DPPC': (['PO4'], ['C4B', 'C4A']), 'CHOL': (['ROH'], ['R5'])}" \
     --k-value 25 \
     --parallel \
     --verbose

聚类分析 (cluster.py)
~~~~~~~~~~~~~~~~~~~~~

**功能描述**
分析脂质分子的聚集行为和聚类模式。

**特定参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- `--cutoff`: 聚类截止距离（默认：8.0埃）

   </div>

**使用示例**

.. code:: bash

   # 基本聚类分析
   python analysis/cluster.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/cluster_results.csv \
     --residues "{'DPPC': ['PO4']}" \
     --verbose

   # 高级聚类分析
   python analysis/cluster.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/cluster_advanced.csv \
     --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" \
     --cutoff 10.0 \
     --start-frame 0 \
     --stop-frame 1000 \
     --parallel \
     --verbose

各向异性分析 (anisotropy.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
计算分子取向的各向异性参数。

**使用示例**

.. code:: bash

   # 各向异性分析
   python analysis/anisotropy.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/anisotropy_results.csv \
     --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" \
     --parallel \
     --verbose

回转半径分析 (gyration.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
计算分子的回转半径，反映分子的紧凑程度。

**使用示例**

.. code:: bash

   # 回转半径分析
   python analysis/gyration.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/gyration_results.csv \
     --residues "{'DPPC': ['PO4'], 'CHOL': ['ROH']}" \
     --parallel \
     --verbose

Sz序参数分析 (sz.py)
~~~~~~~~~~~~~~~~~~~~

**功能描述**
计算脂质链的Sz序参数，反映链的有序程度。

**特定参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- `--chain`: 链类型（'sn1', 'sn2', 或 'both'）
- `--k-value`: Sz计算的k值（默认：15）

   </div>

**使用示例**

.. code:: bash

   # sn1链序参数分析
   python analysis/sz.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/sz_sn1.csv \
     --residues "{'DPPC': ['PO4'], 'DUPC': ['PO4']}" \
     --chain sn1 \
     --k-value 15 \
     --verbose

   # 双链序参数分析
   python analysis/sz.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/sz_both.csv \
     --residues "{'DPPC': ['PO4']}" \
     --chain both \
     --k-value 20 \
     --parallel \
     --verbose

N-聚类分析 (n_cluster.py)
~~~~~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
统计聚类数量，分析聚集模式。

**特定参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- `--cutoff`: 聚类截止距离（默认：12.0埃）
- `--n-cutoff`: 最小聚类大小阈值（默认：10）

   </div>

**使用示例**

.. code:: bash

   # N-聚类分析
   python analysis/n_cluster.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/ncluster_results.csv \
     --residues "{'DAPC': ['GL1', 'GL2'], 'DPPC': ['PO4']}" \
     --cutoff 12.0 \
     --n-cutoff 10 \
     --parallel \
     --verbose

径向分布分析 (rad.py)
~~~~~~~~~~~~~~~~~~~~~~

**功能描述**
计算径向分布函数，分析分子间的距离分布。

**特定参数**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

- `--output-excel`: 输出Excel文件路径
- `--n-circle`: 径向分析的同心圆数量（默认：50）

   </div>

**使用示例**

.. code:: bash

   # 径向分布分析
   python analysis/rad.py \
     --gro-file cases/lnb.gro \
     --output-excel results/radial_distribution.xlsx \
     --residues "{'DPPC': ['NC3'], 'CHOL': ['ROH']}" \
     --n-circle 50

批量处理
--------

脚本自动化
~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #388e3c; margin-top: 0;">📜 批量处理脚本</h3>
   <p>使用脚本实现自动化批量分析：</p>
   </div>

**Python脚本示例**

.. code:: python

   #!/usr/bin/env python3
   """
   批量分析脚本示例
   """
   import os
   import subprocess
   import glob
   from pathlib import Path

   def run_analysis(gro_file, xtc_file, output_dir):
       """运行完整的分析流程"""
       
       # 创建输出目录
       os.makedirs(output_dir, exist_ok=True)
       
       # 分析参数
       residues = "{'DPPC': ['PO4'], 'CHOL': ['ROH']}"
       
       # 分析列表
       analyses = [
           ('pca', 'analysis/pca.py'),
           ('area', 'analysis/area.py'),
           ('curvature', 'analysis/curvature.py'),
           ('cluster', 'analysis/cluster.py'),
       ]
       
       for analysis_name, script_path in analyses:
           output_file = os.path.join(output_dir, f"{analysis_name}_results.csv")
           
           cmd = [
               'python', script_path,
               '--gro-file', gro_file,
               '--xtc-file', xtc_file,
               '--output-csv', output_file,
               '--residues', residues,
               '--parallel',
               '--verbose'
           ]
           
           print(f"运行 {analysis_name} 分析...")
           try:
               subprocess.run(cmd, check=True)
               print(f"{analysis_name} 分析完成")
           except subprocess.CalledProcessError as e:
               print(f"{analysis_name} 分析失败: {e}")

   def main():
       """主函数"""
       # 数据文件路径
       gro_files = glob.glob("data/*.gro")
       xtc_files = glob.glob("data/*.xtc")
       
       for gro_file in gro_files:
           # 找到对应的xtc文件
           base_name = Path(gro_file).stem
           xtc_file = f"data/{base_name}.xtc"
           
           if os.path.exists(xtc_file):
               output_dir = f"results/{base_name}"
               print(f"分析 {base_name}...")
               run_analysis(gro_file, xtc_file, output_dir)
           else:
               print(f"未找到对应的xtc文件: {xtc_file}")

   if __name__ == "__main__":
       main()

**Shell脚本示例**

.. code:: bash

   #!/bin/bash
   # 批量分析Shell脚本

   # 设置参数
   GRO_DIR="data"
   XTC_DIR="data"
   OUTPUT_DIR="results"
   RESIDUES="{'DPPC': ['PO4'], 'CHOL': ['ROH']}"

   # 创建输出目录
   mkdir -p $OUTPUT_DIR

   # 遍历所有gro文件
   for gro_file in $GRO_DIR/*.gro; do
       if [ -f "$gro_file" ]; then
           # 获取文件名（不含扩展名）
           base_name=$(basename "$gro_file" .gro)
           xtc_file="$XTC_DIR/${base_name}.xtc"
           
           if [ -f "$xtc_file" ]; then
               echo "分析 $base_name..."
               
               # 创建子目录
               mkdir -p "$OUTPUT_DIR/$base_name"
               
               # 运行各种分析
               echo "  PCA分析..."
               python analysis/pca.py \
                   --gro-file "$gro_file" \
                   --xtc-file "$xtc_file" \
                   --output-csv "$OUTPUT_DIR/$base_name/pca_results.csv" \
                   --residues "$RESIDUES" \
                   --parallel --verbose
               
               echo "  面积分析..."
               python analysis/area.py \
                   --gro-file "$gro_file" \
                   --xtc-file "$xtc_file" \
                   --output-csv "$OUTPUT_DIR/$base_name/area_results.csv" \
                   --residues "$RESIDUES" \
                   --parallel --verbose
               
               echo "  曲率分析..."
               python analysis/curvature.py \
                   --gro-file "$gro_file" \
                   --xtc-file "$xtc_file" \
                   --output-csv "$OUTPUT_DIR/$base_name/curvature_results.csv" \
                   --residues "$RESIDUES" \
                   --method mean \
                   --parallel --verbose
               
               echo "  $base_name 分析完成"
           else
               echo "未找到对应的xtc文件: $xtc_file"
           fi
       fi
   done

   echo "所有分析完成！"

参数优化
--------

k值优化
~~~~~~~~

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #7b1fa2; margin-top: 0;">🎯 k值自动优化</h3>
   <p>使用机器学习技术自动寻找最佳k值：</p>
   </div>

**k值优化脚本**

.. code:: python

   #!/usr/bin/env python3
   """
   k值优化脚本
   """
   from machine_learning import KValueOptimizer
   import json

   def optimize_k_values():
       """优化不同分析类型的k值"""
       
       # 分析类型和参数
       analyses = {
           'area': {
               'gro_file': 'cases/lnb.gro',
               'xtc_file': 'cases/md.xtc',
               'residues': {'DPPC': ['PO4']}
           },
           'curvature': {
               'gro_file': 'cases/lnb.gro',
               'xtc_file': 'cases/md.xtc',
               'residues': {'DPPC': ['PO4']}
           },
           'height': {
               'gro_file': 'cases/lnb.gro',
               'xtc_file': 'cases/md.xtc',
               'residues': {'DPPC': (['PO4'], ['C4B', 'C4A'])}
           }
       }
       
       optimized_params = {}
       
       for analysis_type, params in analyses.items():
           print(f"优化 {analysis_type} 的k值...")
           
           # 创建优化器
           optimizer = KValueOptimizer(analysis_type)
           
           # 运行优化
           best_k = optimizer.optimize(**params)
           
           optimized_params[analysis_type] = best_k
           print(f"{analysis_type} 最佳k值: {best_k}")
       
       # 保存优化结果
       with open('optimized_k_values.json', 'w') as f:
           json.dump(optimized_params, f, indent=2)
       
       print("k值优化完成！结果已保存到 optimized_k_values.json")

   if __name__ == "__main__":
       optimize_k_values()

**使用优化后的参数**

.. code:: python

   #!/usr/bin/env python3
   """
   使用优化后的参数进行分析
   """
   import json
   import subprocess

   def load_optimized_params():
       """加载优化后的参数"""
       with open('optimized_k_values.json', 'r') as f:
           return json.load(f)

   def run_optimized_analysis():
       """使用优化参数运行分析"""
       
       params = load_optimized_params()
       
       # 使用优化后的k值进行分析
       for analysis_type, k_value in params.items():
           print(f"使用优化k值 {k_value} 运行 {analysis_type} 分析...")
           
           cmd = [
               'python', f'analysis/{analysis_type}.py',
               '--gro-file', 'cases/lnb.gro',
               '--xtc-file', 'cases/md.xtc',
               '--output-csv', f'results/{analysis_type}_optimized.csv',
               '--residues', "{'DPPC': ['PO4']}",
               '--k-value', str(k_value),
               '--parallel',
               '--verbose'
           ]
           
           subprocess.run(cmd, check=True)

   if __name__ == "__main__":
       run_optimized_analysis()

性能优化
--------

并行处理优化
~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #03a9f4; margin-top: 0;">⚡ 并行处理优化</h3>
   <p>优化并行处理性能：</p>
   </div>

**性能测试脚本**

.. code:: python

   #!/usr/bin/env python3
   """
   并行处理性能测试
   """
   import time
   import subprocess
   import multiprocessing
   import psutil

   def test_parallel_performance():
       """测试不同并行数的性能"""
       
       # 测试参数
       gro_file = "cases/lnb.gro"
       xtc_file = "cases/md.xtc"
       residues = "{'DPPC': ['PO4']}"
       
       # 获取CPU核心数
       cpu_count = multiprocessing.cpu_count()
       print(f"系统CPU核心数: {cpu_count}")
       
       # 测试不同并行数
       test_jobs = [1, 2, 4, 8, cpu_count]
       results = {}
       
       for n_jobs in test_jobs:
           if n_jobs > cpu_count:
               continue
               
           print(f"测试 {n_jobs} 个并行作业...")
           
           start_time = time.time()
           
           cmd = [
               'python', 'analysis/pca.py',
               '--gro-file', gro_file,
               '--xtc-file', xtc_file,
               '--output-csv', f'test_pca_{n_jobs}jobs.csv',
               '--residues', residues,
               '--parallel',
               '--n-jobs', str(n_jobs),
               '--start-frame', '0',
               '--stop-frame', '100'
           ]
           
           try:
               subprocess.run(cmd, check=True, capture_output=True)
               end_time = time.time()
               
               execution_time = end_time - start_time
               results[n_jobs] = execution_time
               
               print(f"  {n_jobs} 个作业: {execution_time:.2f} 秒")
               
           except subprocess.CalledProcessError as e:
               print(f"  {n_jobs} 个作业失败: {e}")
       
       # 分析结果
       print("\n性能分析结果:")
       best_jobs = min(results.keys(), key=lambda k: results[k])
       print(f"最佳并行数: {best_jobs}")
       print(f"最佳执行时间: {results[best_jobs]:.2f} 秒")
       
       # 计算加速比
       serial_time = results[1]
       for n_jobs, exec_time in results.items():
           speedup = serial_time / exec_time
           efficiency = speedup / n_jobs * 100
           print(f"{n_jobs} 个作业: 加速比 {speedup:.2f}, 效率 {efficiency:.1f}%")

   if __name__ == "__main__":
       test_parallel_performance()

**内存使用监控**

.. code:: python

   #!/usr/bin/env python3
   """
   内存使用监控脚本
   """
   import psutil
   import time
   import subprocess
   import threading

   class MemoryMonitor:
       def __init__(self):
           self.monitoring = False
           self.max_memory = 0
           self.memory_history = []
       
       def start_monitoring(self):
           """开始监控内存使用"""
           self.monitoring = True
           self.max_memory = 0
           self.memory_history = []
           
           monitor_thread = threading.Thread(target=self._monitor_loop)
           monitor_thread.daemon = True
           monitor_thread.start()
       
       def stop_monitoring(self):
           """停止监控"""
           self.monitoring = False
       
       def _monitor_loop(self):
           """监控循环"""
           while self.monitoring:
               memory_percent = psutil.virtual_memory().percent
               self.memory_history.append(memory_percent)
               self.max_memory = max(self.max_memory, memory_percent)
               time.sleep(1)
       
       def get_stats(self):
           """获取统计信息"""
           if not self.memory_history:
               return None
           
           return {
               'max_memory': self.max_memory,
               'avg_memory': sum(self.memory_history) / len(self.memory_history),
               'min_memory': min(self.memory_history)
           }

   def run_analysis_with_monitoring():
       """带内存监控的分析"""
       
       monitor = MemoryMonitor()
       monitor.start_monitoring()
       
       print("开始分析（监控内存使用）...")
       
       cmd = [
           'python', 'analysis/pca.py',
           '--gro-file', 'cases/lnb.gro',
           '--xtc-file', 'cases/md.xtc',
           '--output-csv', 'monitored_analysis.csv',
           '--residues', "{'DPPC': ['PO4']}",
           '--parallel',
           '--verbose'
       ]
       
       try:
           subprocess.run(cmd, check=True)
       finally:
           monitor.stop_monitoring()
       
       stats = monitor.get_stats()
       if stats:
           print(f"内存使用统计:")
           print(f"  最大内存使用: {stats['max_memory']:.1f}%")
           print(f"  平均内存使用: {stats['avg_memory']:.1f}%")
           print(f"  最小内存使用: {stats['min_memory']:.1f}%")

   if __name__ == "__main__":
       run_analysis_with_monitoring()

错误处理
--------

常见错误解决
~~~~~~~~~~~~

参数格式错误
^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #ffebee; padding: 15px; border-radius: 8px; border-left: 4px solid #f44336;">

**常见参数错误：**

.. code:: bash

   # 错误：缺少引号
   --residues {'DPPC': ['PO4']}
   
   # 正确：使用双引号
   --residues "{'DPPC': ['PO4']}"
   
   # 错误：文件路径包含空格
   --gro-file /path with spaces/file.gro
   
   # 正确：使用引号包围路径
   --gro-file "/path with spaces/file.gro"

   </div>

文件不存在错误
^^^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

**文件检查脚本：**

.. code:: python

   #!/usr/bin/env python3
   """
   文件检查脚本
   """
   import os
   import sys

   def check_files(gro_file, xtc_file):
       """检查文件是否存在"""
       
       errors = []
       
       if not os.path.exists(gro_file):
           errors.append(f"GRO文件不存在: {gro_file}")
       
       if not os.path.exists(xtc_file):
           errors.append(f"XTC文件不存在: {xtc_file}")
       
       if errors:
           print("文件检查失败:")
           for error in errors:
               print(f"  - {error}")
           return False
       
       print("所有文件检查通过")
       return True

   def main():
       if len(sys.argv) != 3:
           print("用法: python check_files.py <gro_file> <xtc_file>")
           sys.exit(1)
       
       gro_file = sys.argv[1]
       xtc_file = sys.argv[2]
       
       if not check_files(gro_file, xtc_file):
           sys.exit(1)

   if __name__ == "__main__":
       main()

   </div>

内存不足错误
^^^^^^^^^^^^

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0;">

**内存优化策略：**

.. code:: python

   #!/usr/bin/env python3
   """
   内存优化分析脚本
   """
   import psutil
   import subprocess

   def check_memory():
       """检查可用内存"""
       memory = psutil.virtual_memory()
       return memory.available / (1024**3)  # GB

   def run_memory_optimized_analysis():
       """运行内存优化的分析"""
       
       available_memory = check_memory()
       print(f"可用内存: {available_memory:.1f} GB")
       
       # 根据可用内存调整参数
       if available_memory < 4:
           # 低内存：减少并行数和帧数
           n_jobs = 1
           stop_frame = 100
           print("低内存模式：使用单线程，限制帧数")
       elif available_memory < 8:
           # 中等内存：适度并行
           n_jobs = 2
           stop_frame = 500
           print("中等内存模式：使用2个线程")
       else:
           # 高内存：完全并行
           n_jobs = -1
           stop_frame = -1
           print("高内存模式：使用所有可用核心")
       
       cmd = [
           'python', 'analysis/pca.py',
           '--gro-file', 'cases/lnb.gro',
           '--xtc-file', 'cases/md.xtc',
           '--output-csv', 'memory_optimized.csv',
           '--residues', "{'DPPC': ['PO4']}",
           '--parallel',
           '--n-jobs', str(n_jobs),
           '--stop-frame', str(stop_frame),
           '--verbose'
       ]
       
       subprocess.run(cmd, check=True)

   if __name__ == "__main__":
       run_memory_optimized_analysis()

   </div>

最佳实践
--------

脚本组织
~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #388e3c; margin-top: 0;">📁 脚本组织最佳实践</h3>
   <p>建议的脚本组织结构：</p>
   </div>

**目录结构**

.. code:: text

   scripts/
   ├── batch_analysis.py          # 批量分析脚本
   ├── parameter_optimization.py  # 参数优化脚本
   ├── performance_test.py       # 性能测试脚本
   ├── utils/
   │   ├── file_utils.py         # 文件工具函数
   │   ├── analysis_utils.py     # 分析工具函数
   │   └── plot_utils.py          # 绘图工具函数
   └── config/
       ├── analysis_config.json   # 分析配置
       └── system_config.json    # 系统配置

**配置管理**

.. code:: json

   {
     "analysis": {
       "default_residues": {
         "DPPC": ["PO4"],
         "CHOL": ["ROH"]
       },
       "default_params": {
         "k_value": 20,
         "cutoff": 8.0,
         "max_normal_angle": 140
       }
     },
     "system": {
       "max_memory_gb": 16,
       "default_n_jobs": -1,
       "output_dir": "results"
     }
   }

日志记录
~~~~~~~~

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 15px; border-radius: 8px; border-left: 4px solid #03a9f4;">

**日志记录最佳实践：**

.. code:: python

   #!/usr/bin/env python3
   """
   带日志记录的分析脚本
   """
   import logging
   import sys
   from datetime import datetime

   def setup_logging():
       """设置日志记录"""
       
       # 创建日志目录
       log_dir = "logs"
       os.makedirs(log_dir, exist_ok=True)
       
       # 生成日志文件名
       timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
       log_file = f"{log_dir}/analysis_{timestamp}.log"
       
       # 配置日志
       logging.basicConfig(
           level=logging.INFO,
           format='%(asctime)s - %(levelname)s - %(message)s',
           handlers=[
               logging.FileHandler(log_file),
               logging.StreamHandler(sys.stdout)
           ]
       )
       
       return logging.getLogger(__name__)

   def run_logged_analysis():
       """带日志记录的分析"""
       
       logger = setup_logging()
       
       logger.info("开始分析")
       logger.info(f"GRO文件: cases/lnb.gro")
       logger.info(f"XTC文件: cases/md.xtc")
       
       try:
           # 运行分析
           cmd = [
               'python', 'analysis/pca.py',
               '--gro-file', 'cases/lnb.gro',
               '--xtc-file', 'cases/md.xtc',
               '--output-csv', 'results/pca_logged.csv',
               '--residues', "{'DPPC': ['PO4']}",
               '--parallel',
               '--verbose'
           ]
           
           logger.info(f"执行命令: {' '.join(cmd)}")
           
           result = subprocess.run(cmd, capture_output=True, text=True)
           
           if result.returncode == 0:
               logger.info("分析成功完成")
           else:
               logger.error(f"分析失败: {result.stderr}")
               
       except Exception as e:
           logger.error(f"分析过程中发生错误: {e}")

   if __name__ == "__main__":
       run_logged_analysis()

   </div>

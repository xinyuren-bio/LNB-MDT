使用示例
========

本页面提供LNB-MDT的详细使用示例，帮助您快速上手。

基础示例
--------

示例1：PCA分析
~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #f3e5f5; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #7b1fa2; margin-top: 0;">🧬 PCA主成分分析示例</h3>
   <p>分析脂质分子的构象变化：</p>
   </div>

**图形界面操作**

1. **启动程序**
   .. code:: bash

      conda activate LNB-MDT
      python main.py

2. **加载数据**
   - 点击"分析"模块
   - 选择GRO文件：`cases/lnb.gro`
   - 选择XTC文件：`cases/md.xtc`
   - 设置输出路径：`results/pca_results.csv`

3. **配置参数**
   - 残基组：`{'DPPC': ['PO4'], 'CHOL': ['ROH']}`
   - 起始帧：0
   - 结束帧：1000
   - 主成分数：3

4. **运行分析**
   - 点击"下一步"
   - 选择"PCA分析"
   - 点击"运行"

**命令行操作**

.. code:: bash

   # 基本PCA分析
   python analysis/pca.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/pca_basic.csv \
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

**结果解读**

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

**PCA结果文件格式：**

.. code-block:: csv

   # Created by LNB-MDT v1.0
   # PCA Analysis
   # TYPE:Bubble
   # Parameters:{'DPPC': ['PO4'], 'CHOL': ['ROH']}
   Frames,PC1,PC2,PC3
   0,0.787,0.234,0.156
   1,0.801,0.241,0.162
   2,0.800,0.238,0.159
   ...

**结果解读：**
- **PC1, PC2, PC3**: 前三个主成分的值
- **数值变化**: 反映分子构象变化
- **趋势分析**: 可用于识别相变或稳定态

   </div>

示例2：面积分析
~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #388e3c; margin-top: 0;">📏 Voronoi面积分析示例</h3>
   <p>计算脂质分子的Voronoi镶嵌面积：</p>
   </div>

**图形界面操作**

1. **加载数据**
   - 选择GRO和XTC文件
   - 设置输出路径

2. **配置参数**
   - 残基组：`{'DPPC': ['PO4']}`
   - k值：20
   - 最大法线角度：140度

3. **运行分析**
   - 选择"面积分析"
   - 启动分析

**命令行操作**

.. code:: bash

   # 基本面积分析
   python analysis/area.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/area_basic.csv \
     --residues "{'DPPC': ['PO4']}" \
     --k-value 20 \
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

**结果可视化**

.. code:: python

   # 面积分析结果可视化
   import pandas as pd
   import matplotlib.pyplot as plt

   # 读取结果
   data = pd.read_csv('results/area_basic.csv')

   # 创建图表
   plt.figure(figsize=(10, 6))
   plt.plot(data['Frames'], data['Values'], 'b-', linewidth=1)
   plt.xlabel('Frame')
   plt.ylabel('Area (nm²)')
   plt.title('Voronoi Area Analysis')
   plt.grid(True, alpha=0.3)
   plt.tight_layout()
   plt.savefig('area_analysis.png', dpi=300)
   plt.show()

示例3：曲率分析
~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #f57c00; margin-top: 0;">🌊 膜曲率分析示例</h3>
   <p>计算脂质膜的平均曲率和高斯曲率：</p>
   </div>

**平均曲率分析**

.. code:: bash

   # 平均曲率分析
   python analysis/curvature.py \
     --gro-file cases/lnb.gro \
     --xtc-file cases/md.xtc \
     --output-csv results/mean_curvature.csv \
     --residues "{'DPPC': ['PO4']}" \
     --method mean \
     --k-value 20 \
     --parallel \
     --verbose

**高斯曲率分析**

.. code:: bash

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

**曲率对比分析**

.. code:: python

   # 曲率对比分析
   import pandas as pd
   import matplotlib.pyplot as plt
   import numpy as np

   # 读取两种曲率结果
   mean_data = pd.read_csv('results/mean_curvature.csv')
   gaussian_data = pd.read_csv('results/gaussian_curvature.csv')

   # 创建对比图
   fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

   # 平均曲率
   ax1.plot(mean_data['Frames'], mean_data['Values'], 'b-', linewidth=1)
   ax1.set_ylabel('Mean Curvature (nm⁻¹)')
   ax1.set_title('Mean Curvature Analysis')
   ax1.grid(True, alpha=0.3)

   # 高斯曲率
   ax2.plot(gaussian_data['Frames'], gaussian_data['Values'], 'r-', linewidth=1)
   ax2.set_xlabel('Frame')
   ax2.set_ylabel('Gaussian Curvature (nm⁻²)')
   ax2.set_title('Gaussian Curvature Analysis')
   ax2.grid(True, alpha=0.3)

   plt.tight_layout()
   plt.savefig('curvature_comparison.png', dpi=300)
   plt.show()

高级示例
--------

示例4：多分析组合
~~~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #1976d2; margin-top: 0;">🔬 多分析组合示例</h3>
   <p>同时运行多种分析并进行综合比较：</p>
   </div>

**批量分析脚本**

.. code:: python

   #!/usr/bin/env python3
   """
   多分析组合脚本
   """
   import subprocess
   import pandas as pd
   import matplotlib.pyplot as plt
   import numpy as np

   def run_multiple_analyses():
       """运行多种分析"""
       
       # 分析配置
       gro_file = "cases/lnb.gro"
       xtc_file = "cases/md.xtc"
       residues = "{'DPPC': ['PO4'], 'CHOL': ['ROH']}"
       
       analyses = [
           {
               'name': 'PCA',
               'script': 'analysis/pca.py',
               'params': ['--n-components', '3']
           },
           {
               'name': 'Area',
               'script': 'analysis/area.py',
               'params': ['--k-value', '20']
           },
           {
               'name': 'Curvature',
               'script': 'analysis/curvature.py',
               'params': ['--method', 'mean', '--k-value', '20']
           },
           {
               'name': 'Cluster',
               'script': 'analysis/cluster.py',
               'params': ['--cutoff', '8.0']
           }
       ]
       
       results = {}
       
       for analysis in analyses:
           print(f"运行 {analysis['name']} 分析...")
           
           output_file = f"results/{analysis['name'].lower()}_results.csv"
           
           cmd = [
               'python', analysis['script'],
               '--gro-file', gro_file,
               '--xtc-file', xtc_file,
               '--output-csv', output_file,
               '--residues', residues,
               '--parallel',
               '--verbose'
           ] + analysis['params']
           
           try:
               subprocess.run(cmd, check=True)
               results[analysis['name']] = output_file
               print(f"{analysis['name']} 分析完成")
           except subprocess.CalledProcessError as e:
               print(f"{analysis['name']} 分析失败: {e}")
       
       return results

   def create_comparison_plot(results):
       """创建对比图表"""
       
       fig, axes = plt.subplots(2, 2, figsize=(15, 10))
       axes = axes.flatten()
       
       for i, (name, file_path) in enumerate(results.items()):
           if i >= 4:
               break
               
           try:
               data = pd.read_csv(file_path)
               
               axes[i].plot(data['Frames'], data['Values'], linewidth=1)
               axes[i].set_title(f'{name} Analysis')
               axes[i].set_xlabel('Frame')
               axes[i].set_ylabel('Values')
               axes[i].grid(True, alpha=0.3)
               
           except Exception as e:
               print(f"读取 {name} 结果失败: {e}")
       
       plt.tight_layout()
       plt.savefig('results/multiple_analysis_comparison.png', dpi=300)
       plt.show()

   def main():
       """主函数"""
       print("开始多分析组合...")
       
       # 运行分析
       results = run_multiple_analyses()
       
       # 创建对比图
       create_comparison_plot(results)
       
       print("多分析组合完成！")

   if __name__ == "__main__":
       main()

示例5：机器学习优化
~~~~~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #fce4ec; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #c2185b; margin-top: 0;">🤖 机器学习参数优化示例</h3>
   <p>使用机器学习技术自动优化分析参数：</p>
   </div>

**参数优化脚本**

.. code:: python

   #!/usr/bin/env python3
   """
   机器学习参数优化示例
   """
   from machine_learning import AnalysisParameterOptimizer, KValueOptimizer
   import time
   import json

   def optimize_area_analysis():
       """优化面积分析参数"""
       
       print("开始面积分析参数优化...")
       
       # 创建优化器
       optimizer = AnalysisParameterOptimizer('area')
       
       # 定义目标函数
       def objective_function(params):
           try:
               from analysis.area import Area
               
               analyzer = Area(
                   gro_file="cases/lnb.gro",
                   xtc_file="cases/md.xtc",
                   residues={'DPPC': ['PO4']},
                   **params
               )
               
               start_time = time.time()
               results = analyzer.run()
               computation_time = time.time() - start_time
               
               # 目标函数：平衡计算时间和结果质量
               objective = computation_time + len(results) * 0.001
               
               print(f"参数 {params}: 目标值 {objective:.4f}")
               return objective
               
           except Exception as e:
               print(f"参数 {params} 评估失败: {e}")
               return float('inf')
       
       # 运行优化
       results = optimizer.optimize()
       
       print(f"优化完成！")
       print(f"最佳参数: {results['best_parameters']}")
       print(f"最佳得分: {results['best_score']:.4f}")
       
       return results

   def optimize_k_values():
       """优化k值"""
       
       print("开始k值优化...")
       
       # 优化不同分析类型的k值
       analyses = ['area', 'curvature', 'height']
       optimized_k_values = {}
       
       for analysis_type in analyses:
           print(f"优化 {analysis_type} 的k值...")
           
           optimizer = KValueOptimizer(analysis_type)
           
           best_k = optimizer.optimize(
               gro_file="cases/lnb.gro",
               xtc_file="cases/md.xtc",
               residues={'DPPC': ['PO4']}
           )
           
           optimized_k_values[analysis_type] = best_k
           print(f"{analysis_type} 最佳k值: {best_k}")
       
       return optimized_k_values

   def run_optimized_analysis(optimized_params):
       """使用优化参数运行分析"""
       
       print("使用优化参数运行分析...")
       
       # 使用优化后的k值
       for analysis_type, k_value in optimized_params.items():
           print(f"使用k值 {k_value} 运行 {analysis_type} 分析...")
           
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

   def main():
       """主函数"""
       print("开始机器学习参数优化...")
       
       # 优化面积分析参数
       area_results = optimize_area_analysis()
       
       # 优化k值
       k_values = optimize_k_values()
       
       # 保存优化结果
       optimization_results = {
           'area_optimization': area_results,
           'k_values': k_values
       }
       
       with open('results/optimization_results.json', 'w') as f:
           json.dump(optimization_results, f, indent=2)
       
       # 使用优化参数运行分析
       run_optimized_analysis(k_values)
       
       print("机器学习参数优化完成！")

   if __name__ == "__main__":
       main()

示例6：异常检测
~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e1f5fe; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #03a9f4; margin-top: 0;">🔍 异常检测示例</h3>
   <p>使用机器学习技术检测分子动力学轨迹中的异常：</p>
   </div>

**异常检测脚本**

.. code:: python

   #!/usr/bin/env python3
   """
   异常检测示例
   """
   from machine_learning import MDAnomalyDetector
   import matplotlib.pyplot as plt
   import numpy as np

   def detect_anomalies():
       """检测轨迹异常"""
       
       print("开始异常检测...")
       
       # 创建异常检测器
       detector = MDAnomalyDetector(
           method='isolation_forest',
           contamination=0.1
       )
       
       # 分析轨迹
       results = detector.analyze_trajectory(
           gro_file="cases/lnb.gro",
           xtc_file="cases/md.xtc",
           residues={'DPPC': ['PO4'], 'CHOL': ['ROH']},
           start_frame=0,
           stop_frame=1000,
           step_frame=5
       )
       
       # 打印结果
       print(f"总帧数: {len(results['predictions'])}")
       print(f"异常帧数: {results['n_anomalies']}")
       print(f"异常比例: {results['anomaly_ratio']:.2%}")
       
       # 显示异常帧
       anomaly_frames = results['anomaly_indices']
       print(f"异常帧: {anomaly_frames}")
       
       return results

   def visualize_anomalies(results):
       """可视化异常检测结果"""
       
       print("创建异常检测可视化...")
       
       # 创建图表
       fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
       
       # 原始数据
       frames = range(len(results['predictions']))
       predictions = results['predictions']
       probabilities = results['probabilities']
       
       # 预测结果
       ax1.scatter(frames, predictions, c=predictions, cmap='RdYlBu', alpha=0.6)
       ax1.set_ylabel('Prediction (-1: Anomaly, 1: Normal)')
       ax1.set_title('Anomaly Detection Results')
       ax1.grid(True, alpha=0.3)
       
       # 异常概率
       ax2.plot(frames, probabilities, 'b-', linewidth=1)
       ax2.axhline(y=0.5, color='r', linestyle='--', alpha=0.7)
       ax2.set_xlabel('Frame')
       ax2.set_ylabel('Anomaly Probability')
       ax2.set_title('Anomaly Probabilities')
       ax2.grid(True, alpha=0.3)
       
       plt.tight_layout()
       plt.savefig('results/anomaly_detection.png', dpi=300)
       plt.show()

   def analyze_anomaly_patterns(results):
       """分析异常模式"""
       
       print("分析异常模式...")
       
       anomaly_frames = results['anomaly_indices']
       
       if len(anomaly_frames) > 0:
           print("异常模式分析:")
           
           # 计算异常间隔
           intervals = np.diff(anomaly_frames)
           if len(intervals) > 0:
               print(f"异常间隔统计:")
               print(f"  平均间隔: {np.mean(intervals):.1f} 帧")
               print(f"  最小间隔: {np.min(intervals)} 帧")
               print(f"  最大间隔: {np.max(intervals)} 帧")
           
           # 分析异常分布
           total_frames = len(results['predictions'])
           anomaly_ratio = len(anomaly_frames) / total_frames
           
           print(f"异常分布:")
           print(f"  异常比例: {anomaly_ratio:.2%}")
           print(f"  异常集中度: {'高' if anomaly_ratio > 0.2 else '低'}")

   def main():
       """主函数"""
       print("开始异常检测分析...")
       
       # 检测异常
       results = detect_anomalies()
       
       # 可视化结果
       visualize_anomalies(results)
       
       # 分析异常模式
       analyze_anomaly_patterns(results)
       
       print("异常检测分析完成！")

   if __name__ == "__main__":
       main()

VMD集成示例
----------

示例7：VMD可视化
~~~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #388e3c; margin-top: 0;">🎬 VMD可视化集成示例</h3>
   <p>使用VMD进行分子可视化：</p>
   </div>

**VMD集成脚本**

.. code:: python

   #!/usr/bin/env python3
   """
   VMD集成可视化示例
   """
   from modules.vmd_control import VMDTcp, VMDCommands
   import pandas as pd
   import time

   def setup_vmd():
       """设置VMD连接"""
       
       print("启动VMD...")
       
       # VMD配置
       rctl_path = "./remote_ctl.tcl"
       vmd_path = "C:/Program Files/VMD/vmd.exe"  # Windows路径
       
       # 创建VMD连接
       vmd = VMDTcp(rctl_path, vmd_path)
       
       # 启动VMD
       response = vmd.start()
       if response == -1:
           print("VMD启动失败！")
           return None
       
       print("VMD启动成功！")
       return vmd

   def load_trajectory(vmd, gro_file, xtc_file):
       """加载轨迹文件"""
       
       print(f"加载轨迹: {gro_file}, {xtc_file}")
       
       # 加载拓扑文件
       vmd.send_command(VMDCommands.loadTopology(gro_file))
       
       # 加载轨迹文件
       vmd.send_command(VMDCommands.loadTrajectory(xtc_file))
       
       # 设置显示样式
       vmd.send_command(VMDCommands.setRepresentation("CPK"))
       vmd.send_command(VMDCommands.setColoringMethod("Name"))

   def visualize_analysis_results(vmd, csv_file):
       """可视化分析结果"""
       
       print(f"加载分析结果: {csv_file}")
       
       # 读取分析结果
       data = pd.read_csv(csv_file)
       
       # 获取帧数和值
       frames = data['Frames'].tolist()
       values = data['Values'].tolist()
       
       # 找到极值帧
       max_frame = frames[values.index(max(values))]
       min_frame = frames[values.index(min(values))]
       
       print(f"最大值帧: {max_frame}, 值: {max(values):.4f}")
       print(f"最小值帧: {min_frame}, 值: {min(values):.4f}")
       
       # 跳转到极值帧
       print("跳转到最大值帧...")
       vmd.send_command(VMDCommands.gotoFrame(str(max_frame)))
       time.sleep(1)
       
       print("跳转到最小值帧...")
       vmd.send_command(VMDCommands.gotoFrame(str(min_frame)))
       time.sleep(1)
       
       return max_frame, min_frame

   def highlight_molecules(vmd, residues):
       """高亮特定分子"""
       
       print(f"高亮分子: {residues}")
       
       # 高亮DPPC分子
       if 'DPPC' in residues:
           vmd.send_command(VMDCommands.highlightResname("DPPC"))
       
       # 高亮胆固醇分子
       if 'CHOL' in residues:
           vmd.send_command(VMDCommands.highlightResname("CHOL"))

   def create_animation(vmd, start_frame, end_frame, step=10):
       """创建动画"""
       
       print(f"创建动画: 帧 {start_frame} 到 {end_frame}")
       
       for frame in range(start_frame, end_frame, step):
           vmd.send_command(VMDCommands.gotoFrame(str(frame)))
           time.sleep(0.1)

   def main():
       """主函数"""
       print("开始VMD集成可视化...")
       
       # 设置VMD
       vmd = setup_vmd()
       if vmd is None:
           return
       
       try:
           # 加载轨迹
           load_trajectory(vmd, "cases/lnb.gro", "cases/md.xtc")
           
           # 可视化分析结果
           max_frame, min_frame = visualize_analysis_results(vmd, "results/pca_results.csv")
           
           # 高亮分子
           highlight_molecules(vmd, {'DPPC': ['PO4'], 'CHOL': ['ROH']})
           
           # 创建动画
           create_animation(vmd, 0, 100, step=5)
           
           print("VMD可视化完成！")
           
       except Exception as e:
           print(f"VMD可视化过程中发生错误: {e}")
       
       finally:
           # 停止VMD
           vmd.stop()
           print("VMD已停止")

   if __name__ == "__main__":
       main()

性能优化示例
-----------

示例8：性能测试
~~~~~~~~~~~~~~~

.. raw:: html

   <div style="background-color: #fff3e0; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #f57c00; margin-top: 0;">⚡ 性能测试和优化示例</h3>
   <p>测试和优化分析性能：</p>
   </div>

**性能测试脚本**

.. code:: python

   #!/usr/bin/env python3
   """
   性能测试脚本
   """
   import time
   import subprocess
   import multiprocessing
   import psutil
   import matplotlib.pyplot as plt

   class PerformanceTester:
       def __init__(self):
           self.results = {}
       
       def test_parallel_performance(self):
           """测试并行性能"""
           
           print("测试并行性能...")
           
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
           
           self.results['parallel'] = results
           return results
       
       def test_memory_usage(self):
           """测试内存使用"""
           
           print("测试内存使用...")
           
           # 监控内存使用
           memory_before = psutil.virtual_memory().used / (1024**3)
           
           cmd = [
               'python', 'analysis/pca.py',
               '--gro-file', 'cases/lnb.gro',
               '--xtc-file', 'cases/md.xtc',
               '--output-csv', 'test_memory.csv',
               '--residues', "{'DPPC': ['PO4']}",
               '--parallel',
               '--verbose'
           ]
           
           try:
               subprocess.run(cmd, check=True)
           except subprocess.CalledProcessError as e:
               print(f"内存测试失败: {e}")
           
           memory_after = psutil.virtual_memory().used / (1024**3)
           memory_used = memory_after - memory_before
           
           print(f"内存使用: {memory_used:.2f} GB")
           
           self.results['memory'] = memory_used
           return memory_used
       
       def test_different_parameters(self):
           """测试不同参数的性能影响"""
           
           print("测试不同参数的性能影响...")
           
           # 测试不同k值
           k_values = [10, 15, 20, 25, 30]
           k_results = {}
           
           for k_value in k_values:
               print(f"测试k值: {k_value}")
               
               start_time = time.time()
               
               cmd = [
                   'python', 'analysis/area.py',
                   '--gro-file', 'cases/lnb.gro',
                   '--xtc-file', 'cases/md.xtc',
                   '--output-csv', f'test_k{k_value}.csv',
                   '--residues', "{'DPPC': ['PO4']}",
                   '--k-value', str(k_value),
                   '--verbose'
               ]
               
               try:
                   subprocess.run(cmd, check=True, capture_output=True)
                   end_time = time.time()
                   
                   execution_time = end_time - start_time
                   k_results[k_value] = execution_time
                   
                   print(f"  k值 {k_value}: {execution_time:.2f} 秒")
                   
               except subprocess.CalledProcessError as e:
                   print(f"  k值 {k_value} 失败: {e}")
           
           self.results['k_values'] = k_results
           return k_results
       
       def create_performance_plots(self):
           """创建性能图表"""
           
           print("创建性能图表...")
           
           fig, axes = plt.subplots(2, 2, figsize=(15, 10))
           
           # 并行性能图
           if 'parallel' in self.results:
               parallel_data = self.results['parallel']
               jobs = list(parallel_data.keys())
               times = list(parallel_data.values())
               
               axes[0, 0].plot(jobs, times, 'bo-')
               axes[0, 0].set_xlabel('并行作业数')
               axes[0, 0].set_ylabel('执行时间 (秒)')
               axes[0, 0].set_title('并行性能测试')
               axes[0, 0].grid(True, alpha=0.3)
           
           # k值性能图
           if 'k_values' in self.results:
               k_data = self.results['k_values']
               k_values = list(k_data.keys())
               k_times = list(k_data.values())
               
               axes[0, 1].plot(k_values, k_times, 'ro-')
               axes[0, 1].set_xlabel('k值')
               axes[0, 1].set_ylabel('执行时间 (秒)')
               axes[0, 1].set_title('k值性能影响')
               axes[0, 1].grid(True, alpha=0.3)
           
           # 内存使用图
           if 'memory' in self.results:
               memory_used = self.results['memory']
               axes[1, 0].bar(['内存使用'], [memory_used], color='green', alpha=0.7)
               axes[1, 0].set_ylabel('内存使用 (GB)')
               axes[1, 0].set_title('内存使用测试')
           
           # 性能总结
           axes[1, 1].text(0.1, 0.5, '性能测试总结', fontsize=16, fontweight='bold')
           axes[1, 1].text(0.1, 0.3, f'最佳并行数: {min(self.results.get("parallel", {}), key=self.results.get("parallel", {}).get) if "parallel" in self.results else "N/A"}', fontsize=12)
           axes[1, 1].text(0.1, 0.2, f'最佳k值: {min(self.results.get("k_values", {}), key=self.results.get("k_values", {}).get) if "k_values" in self.results else "N/A"}', fontsize=12)
           axes[1, 1].text(0.1, 0.1, f'内存使用: {self.results.get("memory", 0):.2f} GB', fontsize=12)
           axes[1, 1].set_xlim(0, 1)
           axes[1, 1].set_ylim(0, 1)
           axes[1, 1].axis('off')
           
           plt.tight_layout()
           plt.savefig('results/performance_test.png', dpi=300)
           plt.show()

   def main():
       """主函数"""
       print("开始性能测试...")
       
       tester = PerformanceTester()
       
       # 运行各种测试
       tester.test_parallel_performance()
       tester.test_memory_usage()
       tester.test_different_parameters()
       
       # 创建性能图表
       tester.create_performance_plots()
       
       print("性能测试完成！")

   if __name__ == "__main__":
       main()

最佳实践总结
-----------

.. raw:: html

   <div style="background-color: #e3f2fd; padding: 20px; border-radius: 8px; margin: 20px 0;">
   <h3 style="color: #1976d2; margin-top: 0;">📋 最佳实践总结</h3>
   <p>基于以上示例的最佳实践建议：</p>
   </div>

**分析流程建议**

.. raw:: html

   <div style="background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;">

1. **数据预处理**: 检查轨迹质量，去除异常帧
2. **参数优化**: 使用ML模块优化关键参数
3. **批量分析**: 使用脚本自动化分析流程
4. **结果验证**: 检查结果的合理性
5. **可视化**: 使用图表和VMD可视化结果
6. **性能监控**: 定期进行性能测试和优化

   </div>

**代码组织建议**

.. raw:: html

   <div style="background-color: #e8f5e8; padding: 15px; border-radius: 8px; border-left: 4px solid #4caf50;">

- **模块化**: 将功能分解为独立的模块
- **配置管理**: 使用配置文件管理参数
- **错误处理**: 添加适当的错误处理机制
- **日志记录**: 记录分析过程和结果
- **文档化**: 为脚本添加详细的注释

   </div>

**性能优化建议**

.. raw:: html

   <div style="background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800;">

- **并行处理**: 合理使用并行处理提高效率
- **内存管理**: 监控内存使用，避免内存溢出
- **参数调优**: 根据系统特点调整分析参数
- **硬件优化**: 使用SSD和充足内存

   </div>

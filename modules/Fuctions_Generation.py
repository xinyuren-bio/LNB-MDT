import os
import subprocess
import sys
import shutil

from main import *
from functools import partial
from generation.lipidsInfo_martini3 import *
from .Tools import *
global window


class InfoGeneration:
    BOX_D: str = '1'
    TEST_PATH: str = 'E:/excel/1.gro'

    def __init__(self, ui):
        self.ui = ui
        # box
        self.box_x = None
        self.box_y = None
        self.box_z = None
        self.box_d = self.BOX_D
        # path
        self.pathGene = None
        # lnb
        self.gas_density = None
        self.gas_type = None
        self.area = None
        self.r = None
        # solvent
        self.solvent_type = None
        self.salt = None
        self.lipids = None

        # 函数执行
        self._get_info()

    def _get_info(self):
        # box
        self.box_x = str(self.ui.spin_box_x.value())
        self.box_y = str(self.ui.spin_box_y.value())
        self.box_z = str(self.ui.spin_box_z.value())
        # path - 如果用户输入的是文件名（不包含.gro后缀），则添加.gro后缀
        user_input = str(self.ui.edit_gene_path.text() if self.ui.edit_gene_path.text() else self.TEST_PATH)
        if not user_input.endswith('.gro'):
            # 如果用户只输入了文件名，添加.gro后缀
            user_input = user_input + '.gro'
        self.pathGene = user_input
        # lnb
        self.gas_density = str(self.ui.spin_gas_density.value())
        self.gas_type = str(self.ui.como_gas.currentText())
        self.area = str(self.ui.spin_area_5.value())
        self.r = str(self.ui.spin_r.value())
        # solvent
        self.solvent_type = str(self.ui.como_solvent.currentText())
        self.salt = str(self.ui.spin_salt.value())
        self.lipids = {}
        for i in ALL_P:
            value_p = getattr(self.ui, f'Gene{i}Spin').value()
            if value_p != 0:
                self.lipids[i] = value_p


class BtnGeneClick:
    LNB_PATH: str = 'generation//lnb_gener_martini3.py'

    def __init__(self, ui):
        self.ui = ui
        try:
            self.Info = InfoGeneration(self.ui)
        except:
            create_warn_dialog(text="Please check the enter infomation", title="Warning")
            return
        self.runClick()

    def runClick(self):
        # 使用当前解释器的 Python 路径
        command = [sys.executable, self.LNB_PATH]
        arguments = self.getArguments()
        full_command = command + arguments

        result = subprocess.run(full_command, capture_output=True, text=True)

        if result.stderr:
            folder_path = self.writeTop(result.stderr)
            success_message = (
                'system infomation:\n'
                f"{result.stderr}\n"
                "Files were organized and saved in folder:\n"
                f"{folder_path if folder_path else 'Failed to create folder structure'}\n\n"
                "Folder structure:\n"
                "  - system.gro\n"
                "  - system.top\n"
                "  - parameter.txt (generation parameters)\n"
                "  - README\n"
                "  - *.mdp files (simulation parameters)\n"
                "  - toppar/\n"
                "    - *.itp files (topology parameters)"
            )
            create_warn_dialog(text=success_message, title="Generation")


    def getArguments(self):
        arguments = [
            '-d', self.Info.box_d,
            '-r', self.Info.r,
            '-x', self.Info.box_x,
            '-y', self.Info.box_y,
            '-z', self.Info.box_z,
            '-sol', self.Info.solvent_type,
            '-salt', self.Info.salt,
            '-a', self.Info.area,
            '-gas', self.Info.gas_type,
            '-gden', self.Info.gas_density,
            '-o', self.Info.pathGene
        ]
        for p in self.Info.lipids:
            arguments.append('-u')
            arguments.append(f'{p}:{self.Info.lipids[p]}')
        return arguments

    def writeTop(self, additional_text):
        itps = '\n'.join(f'#include  "toppar/{i}.itp"' for i in self.Info.lipids)

        content = """;
;  
; Example topology file for MARTINI 3 
;  

; First include the file containing all particle definitions,  
; the interaction matrix, plus the topology for water.  


; Then include the file(s) containing the topologies of other  
; molecules present in your system.  

#include "toppar/martini_v3.0.0.itp"
#include "toppar/martini_v3.0.0_ions.itp"  
#include "toppar/martini_v3.0.0_phospholipids.itp"
#include "toppar/martini_v3.0.0_gas.itp" 
#include "toppar/martini_v3.0.0_solvents.itp"


; Define a name for your system  

[ system ]  
Lipid-Nanobubble  

; Define the composition of your system  
; The molecule names should correspond to those defined in the itp file(s).  

[ molecules ]
{}  

""".format(additional_text)
        with open(f'{os.path.dirname(self.Info.pathGene)}/topol.top', "w") as file:
            file.write(content)
        
        # 创建文件夹结构并复制文件
        self.create_folder_structure()

    def create_parameter_file(self, main_folder_path):
        """创建包含所有参数的parameter.txt文件"""
        parameter_content = f"""LNB-MDT Generation Parameters
===============================

Generation Time: {subprocess.run(['date'], capture_output=True, text=True).stdout.strip()}

System Parameters:
-----------------
Box Dimensions:
  - X: {self.Info.box_x}
  - Y: {self.Info.box_y}
  - Z: {self.Info.box_z}
  - D: {self.Info.box_d}

LNB Parameters:
---------------
Gas Type: {self.Info.gas_type}
Gas Density: {self.Info.gas_density}
Area: {self.Info.area}
Radius: {self.Info.r}

Solvent Parameters:
------------------
Solvent Type: {self.Info.solvent_type}
Salt Concentration: {self.Info.salt}

Lipid Composition:
-----------------"""
        
        # 添加脂质信息
        if self.Info.lipids:
            for lipid, count in self.Info.lipids.items():
                parameter_content += f"\n  - {lipid}: {count}"
        else:
            parameter_content += "\n  - No lipids selected"
        
        parameter_content += f"""

Output Information:
------------------
Output Path: {self.Info.pathGene}
Base Filename: {os.path.splitext(os.path.basename(self.Info.pathGene))[0]}

Generation Command:
------------------"""
        
        # 添加生成命令参数
        arguments = self.getArguments()
        command_line = "python lnb_gener_martini3.py " + " ".join(arguments)
        parameter_content += f"\n{command_line}\n"
        
        parameter_content += """
Notes:
------
- This file contains all parameters used for system generation
- Generated by LNB-MDT v1.0
- For questions or support, please refer to the documentation
"""
        
        # 写入parameter.txt文件
        parameter_file_path = os.path.join(main_folder_path, "parameter.txt")
        try:
            with open(parameter_file_path, 'w', encoding='utf-8') as f:
                f.write(parameter_content)
            return parameter_file_path
        except Exception as e:
            print(f"创建parameter.txt文件时出错: {e}")
            return None

    def create_folder_structure(self):
        """创建文件夹结构并复制所有文件"""
        # 获取输出目录和文件名
        output_dir = os.path.dirname(self.Info.pathGene)
        gro_filename = os.path.basename(self.Info.pathGene)
        base_filename = os.path.splitext(gro_filename)[0]
        
        # 创建主文件夹
        main_folder_name = f"LNB_MDT_{base_filename}"
        main_folder_path = os.path.join(output_dir, main_folder_name)
        
        # 创建主文件夹
        os.makedirs(main_folder_path, exist_ok=True)
        
        # 创建toppar子文件夹
        toppar_folder_path = os.path.join(main_folder_path, "toppar")
        os.makedirs(toppar_folder_path, exist_ok=True)
        
        try:
            # 重命名并移动gro文件为system.gro
            gro_path = os.path.join(output_dir, gro_filename)
            system_gro_path = os.path.join(main_folder_path, "system.gro")
            if os.path.exists(gro_path):
                shutil.copy2(gro_path, system_gro_path)
                # 删除原始gro文件
                os.remove(gro_path)

            # 重命名并移动top文件为system.top
            top_path = os.path.join(output_dir, "topol.top")
            system_top_path = os.path.join(main_folder_path, "system.top")
            if os.path.exists(top_path):
                shutil.copy2(top_path, system_top_path)
                # 删除原始top文件
                os.remove(top_path)

            # 移动ndx文件为system.ndx
            ndx_path = os.path.join(output_dir, f"{base_filename}.ndx")
            system_ndx_path = os.path.join(main_folder_path, "system.ndx")
            if os.path.exists(ndx_path):
                shutil.copy2(ndx_path, system_ndx_path)
                os.remove(ndx_path)

            # 获取files文件夹路径
            files_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'generation', 'files')
            
            # 复制文件到相应位置
            if os.path.exists(files_dir):
                for file in os.listdir(files_dir):
                    source_path = os.path.join(files_dir, file)
                    if file.endswith('.itp'):
                        # .itp文件复制到toppar文件夹
                        dest_path = os.path.join(toppar_folder_path, file)
                        shutil.copy2(source_path, dest_path)
                    else:
                        # 其他文件（README, .mdp等）复制到主文件夹
                        dest_path = os.path.join(main_folder_path, file)
                        shutil.copy2(source_path, dest_path)
            
            # 创建参数文件
            self.create_parameter_file(main_folder_path)
            
            return main_folder_path
            
        except Exception as e:
            print(f"创建文件夹结构时出错: {e}")
            return None

def lipidsSelect(ui):
    ui.GeneExtraLayout = QVBoxLayout(ui.Gene_extraBox)
    ui.GeneExtraLayout.addWidget(GeneWidgetLipidsSel(ui))


class GeneWidgetLipidsSel(QWidget):
    def __init__(self, ui):
        super().__init__()
        self.ui = ui
        self.ui.GenemainLayout = QVBoxLayout()
        self.LipidsLayout()

    def LipidsLayout(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)

        scrollArea = QScrollArea()
        scrollArea.setWidgetResizable(True)
        container = QWidget()
        containerLayout = QVBoxLayout()  # 创建布局

        for type in lipids:
            groupBox = UIItemsMake.make_widget()
            setattr(self.ui,f'Gene{type}Widget', groupBox)
            btn = UIItemsMake.make_btn(btnName='▼', background_color='#6272a4')  # 浅蓝色
            setattr(self.ui,f'Genebtn{type}', btn)
            groupBoxLayout = QGridLayout(groupBox)
            label = UIItemsMake.make_label(type)

            frame = QFrame()
            hlayout = QHBoxLayout(frame)
            hlayout.addWidget(label)
            hlayout.addWidget(btn)
            btn.clicked.connect(partial(self.toggle_container, groupBox, btn))

            for lipid in lipids[type]:
                label = UIItemsMake.make_label(lipid)
                spin = UIItemsMake.make_spin_box(value=0, max=100)
                setattr(self.ui, f'Gene{lipid}Spin', spin)
                groupBoxLayout.addWidget(label, lipids[type].index(lipid), 0)
                groupBoxLayout.addWidget(spin, lipids[type].index(lipid), 1)
            groupBox.hide()
            containerLayout.addWidget(frame)
            containerLayout.addWidget(groupBox)
        container.setLayout(containerLayout)  # 为container设置布局
        scrollArea.setWidget(container)  # 将container设置为scrollArea的子部件

        layout.addWidget(scrollArea)
        self.ui.GenemainLayout.addWidget(widget)
        self.setLayout(self.ui.GenemainLayout)

    def toggle_container(self, container, btn):
        if container.isVisible():
            container.hide()
            btn.setText('▼')
        else:
            container.show()
            btn.setText('▲')
import os
import subprocess
import sys

from main import *
from functools import partial
from generation.lipidsInfo import *
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
        # path
        self.pathGene = str(self.ui.edit_gene_path.text() if self.ui.edit_gene_path.text() else self.TEST_PATH)
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
    LNB_PATH: str = 'generation//lnb_gener.py'

    def __init__(self, ui):
        self.ui = ui
        self.Info = InfoGeneration(self.ui)
        self.runClick()

    def runClick(self):
        # 使用当前解释器的 Python 路径
        command = [sys.executable, self.LNB_PATH]
        arguments = self.getArguments()
        full_command = command + arguments

        result = subprocess.run(full_command, capture_output=True, text=True)
        if result.stderr:
            print(result.stderr)
            self.writeTop(result.stderr)
        print(result.stdout)

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
; Example topology file for MARTINI 2.1  
;  

; First include the file containing all particle definitions,  
; the interaction matrix, plus the topology for water.  


; Then include the file(s) containing the topologies of other  
; molecules present in your system.  

#include "toppar/martini_v2.1-gas.itp" 
{} 
#include "toppar/martini_v2.0_ions.itp"  


; Define a name for your system  

[ system ]  
Lipid-Nanobubble  

; Define the composition of your system  
; The molecule names should correspond to those defined in the itp file(s).  

{}  

""".format(itps, additional_text)
        with open(f'{os.path.dirname(self.Info.pathGene)}/topol.top', "w") as file:
            file.write(content)

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
import time
from collections import defaultdict
from threading import Thread
from dataclasses import dataclass, field

import numpy as np
from PySide6.QtCore import QObject, Signal, QPropertyAnimation, QEasingCurve, QTimer, QThread
from PySide6.QtGui import QFont, QColor
import MDAnalysis as mda


from figure import *
from .Tools import *
from analysis import *

__all__ = ['NextClick']


# Next Click
# ///////////////////////////////////////////////////////////////

def NextClick(ui):
    if not hasattr(ui, 'VLayoutRightMain'):
        AnalysisBtnClick(ui)
    else:
        if ui.extraRightBox.width() == 0:
            AnalysisBtnClick(ui)
        else:
            ui.extraRightBox.setLayout(None)
            while ui.VLayoutRightMain.count():
                item = ui.VLayoutRightMain.takeAt(0)
                if item.widget():
                    item.widget().deleteLater()
                elif item.layout():
                    item.layout().deleteLater()


class AnalysisUtils:
    """静态工具类，用于处理残基和原子相关逻辑。"""

    @staticmethod
    def get_residue(residues_list, residue_click):
        residue_click.clear()
        for box in residues_list:
            if box.isChecked():
                residue_click.append(box.text())

    @staticmethod
    def get_atom(Group, btnType, storeList):
        storeList.clear()
        for group in Group:
            atoms = []
            for box in group.findChildren(btnType):
                if box.isChecked():
                    atoms.append(box.text())
            storeList[group.title()] = atoms
        print('dict', storeList)

    @staticmethod
    def get_key_value(dict1, dict2):
        dict3 = {}
        for key in dict1.keys() & dict2.keys():
            atoms_1 = dict1[key]
            atoms_2 = dict2[key]
            atomsAdd = (atoms_1, atoms_2)
            dict3[key] = atomsAdd
        return dict3

    @staticmethod
    def get_spin_value(spin_box):
        return spin_box.value()


class AnalysisBtnClick:
    """
    创建出右侧的窗口，并且显示残基，后续不同的runmethod对应
    的不同的显示，由AnalysisLayout完成
    """
    FONT_SIZE = '15pt'  # 字体大小
    TEST_GRO_PATH = 'E:/ach.gro'
    TEST_XTC_PATH = 'E:/ach.xtc'

    def __init__(self, ui):
        self.ui = ui
        self.readFile()

    def readFile(self):
        self.Info = InfoAnalysis(self.ui)
        self.method = self.Info.method_analysis
        self.Box = SelBox()
        # try:
        #     self.Box.u = mda.Universe(self.Info.path_structure
        #                               , self.Info.path_trajectory
        #                               , all_coordinates=False)
        #     self.Box.residues = np.unique(self.Box.u.atoms.resnames)
        #     for sp in self.Box.residues:
        #         self.Box.residues_atoms[sp] = np.unique(self.Box.u.select_atoms('resname %s' % sp).names)
        #     self.stackLayout()
        # except:
        #     create_warn_dialog(text='Failed to import files')

        self.Box.u = mda.Universe(self.Info.path_structure
                                  , self.Info.path_trajectory
                                  , all_coordinates=False)
        self.Box.residues = np.unique(self.Box.u.atoms.resnames)
        for sp in self.Box.residues:
            self.Box.residues_atoms[sp] = np.unique(self.Box.u.select_atoms('resname %s' % sp).names)
        self.stackLayout()

            # self.Box.u = mda.Universe(self.TEST_GRO_PATH, self.TEST_XTC_PATH)

    def stackLayout(self):
        # 创建叠加窗口
        self.ui.stackedWidget_Analysis = QStackedWidget()
        self.ui.stackedWidget_Analysis.setStyleSheet("font-size:%s;" % self.FONT_SIZE)
        # 右侧窗口添加垂直布局
        if not hasattr(self.ui, 'VLayoutRightMain'):
            self.ui.VLayoutRightMain = QVBoxLayout(self.ui.extraRightBox)

        # 创建顶部窗口，储存Next和Back按钮
        self.ui.widgetUpNextBack = QWidget()
        # 创建顶部布局用来刷新事件
        self.ui.btnBack = UIItemsMake.make_btn('Back'
                                                , background_color='white'
                                                , font_color='black'
                                                )
        self.ui.btnBack.clicked.connect(self.switchWidgetBack)
        self.ui.btnRefresh = UIItemsMake.make_btn('Refresh'
                                              , background_color='white'
                                              , font_color='black'
                                              )
        self.ui.btnRefresh.clicked.connect(self.refreshWidget)
        # 添加水平布局到顶部窗口
        HLayoutNextBack = QHBoxLayout(self.ui.widgetUpNextBack)
        HLayoutNextBack.addWidget(self.ui.btnBack)
        HLayoutNextBack.addWidget(self.ui.btnRefresh)

        # 创建显示残基的窗口，并添加到StackedWidget以及在该窗口添加垂直布局
        self.ui.widgetResidue = QWidget()
        self.ui.VLayoutResidue = QVBoxLayout(self.ui.widgetResidue)
        self.ui.stackedWidget_Analysis.addWidget(self.ui.widgetResidue)
        # 在垂直布局中添加标签
        self.ui.Label = UIItemsMake.make_label('Select Residues', font_size='16pt')
        self.ui.VLayoutResidue.addWidget(self.ui.Label)
        for sp in self.Box.residues:
            checkBox = QCheckBox(sp)
            self.ui.VLayoutResidue.addWidget(checkBox)
            self.Box.residues_list.append(checkBox)

        layout_strategies = {
            'Height': HeightLayout
            , 'SZ': SZLayout
            , 'Area': AreaLayout
            , 'MeanCurvature': MeanCurvatureLayout
            , 'Anisotropy': AnisotropyLayout
            , 'RadialDistribution': RDLayout
            , 'Cluster': ClusterLayout
            , 'NCluster': NClusterLayout
            , 'Gyration': GyrationLayout
            , 'PCA': PCALayout
        }

        layout_analysis = layout_strategies[self.method]
        self._create_layout(layout_analysis)

        self.ui.VLayoutRightMain.addWidget(self.ui.widgetUpNextBack)
        self.ui.VLayoutRightMain.addWidget(self.ui.stackedWidget_Analysis)

    def _create_layout(self, layout):
        init_layout = layout(self.ui, self.Box, self.Info)
        btn = UIItemsMake.make_btn('Next'
                                   , callback=lambda: AnalysisUtils.get_residue(self.Box.residues_list,
                                                                                self.Box.residue_click)
                                   , layout=self.ui.VLayoutResidue)
        btn.clicked.connect(lambda: init_layout.step_1())

    def switchWidgetBack(self):
        currentIndex = self.ui.stackedWidget_Analysis.currentIndex()
        numWidgets = self.ui.stackedWidget_Analysis.count()
        if currentIndex == 0:
            pass
        else:
            try:
                self.ui.progressBar.deleteLater()
                self.ui.btnMakeFigure.deleteLater()
            except:
                pass
            self.ui.stackedWidget_Analysis.setCurrentIndex(currentIndex - 1)
            for i in range(currentIndex, numWidgets):
                widget = self.ui.stackedWidget_Analysis.widget(i)
                self.ui.stackedWidget_Analysis.removeWidget(widget)
                if widget:
                    widget.deleteLater()
                else:
                    pass

    def refreshWidget(self):
        def clearLayout(layout):
            """
            清空布局中的所有控件并删除它们
            """
            if layout is not None:
                while layout.count():
                    child = layout.takeAt(0)
                    if child.widget():
                        child.widget().deleteLater()  # 删除控件
                    elif child.layout():
                        clearLayout(child.layout())  # 如果是布局，则递归清理子布局
        clearLayout(self.ui.VLayoutRightMain)
        self.readFile()


class Worker(QObject):
    progressValueChanged = Signal(int)

    def __init__(self, cls_instance, *args, **kwargs):
        super().__init__()
        self.cls_instance = cls_instance
        self.args = args
        self.kwargs = kwargs

    def run(self):
        # 在此调用具体的分析方法

        self.cls_instance.run(*self.args, **self.kwargs, callBack=self.update_progress)
        self.update_progress(100)
        # 防止线程出现问题
        time.sleep(1)


    def update_progress(self, value):
        self.progressValueChanged.emit(value)


class AnalysisLayout:
    def __init__(self, ui, Box, Info):
        self.ui = ui
        self.Box = Box
        self.Info = Info
        self.start_time = None  # 初始化开始时间

    @classmethod
    def _addProgressBar(cls, func):
        def wrapper(self, *args, **kwargs):
            if not hasattr(self, 'progressBar'):
                self.ui.progressBar = QProgressBar()
                self.ui.progressBar.setStyleSheet("""
                                   QProgressBar {
                                       border: 2px solid grey;
                                       border-radius: 5px;
                                   }
                                   QProgressBar::chunk {
                                       background-color: #6272a4;
                                       width: 20px;
                                   }
                               """)
                self.ui.VLayoutRightMain.addWidget(self.ui.progressBar)
            self.start_time = time.time()
            result = func(self, *args, **kwargs)
            worker = Worker(self.cls, self.Info.frame_first, self.Info.frame_last, self.Info.step)
            worker.progressValueChanged.connect(self.updateProgressBar)
            self.thread = Thread(target=worker.run)
            self.thread.start()
            return result

        return wrapper

    def updateProgressBar(self, value):
        self.ui.progressBar.setValue(value)
        if value == 100:
            if not hasattr(self, 'btnMakeFigure'):
                self.ui.btnAnalysisRun.deleteLater()
                self.ui.btnMakeFigure = UIItemsMake.make_btn('Make Figure', background_color='#BD93F9')
                self.ui.btnMakeFigure.clicked.connect(self.makeFigure)
                self.ui.VLayoutRightMain.addWidget(self.ui.btnMakeFigure)
            end_time = time.time()  # 记录结束时间
            elapsed_time = end_time - self.start_time  # 计算耗时
            formatted_time = time.strftime('%M:%S', time.gmtime(elapsed_time))  # 格式化耗时
            success_message = (
                'Analysis Completed\n'
                f"Time taken: {formatted_time}\n"
                "The gro file and topol file were saved at:\n"
                f"{self.Info.path_result}"
            )
            create_warn_dialog(success_message, 'Analysis')

    def AtomsLayout(self
                    , widgetName: str
                    , widgetLayout: str
                    , labelText: str
                    , groups
                    , radioOrCheck
                    , btnName: str
                    , stackID
                    , func
                    , connect=False
                    ):
        """
        :param widgetName:创新新的窗口名称
        :param widgetLayout: 创建新的布局的名称
        :param groups: 储存信息的Group
        :param radioOrCheck:
        :param stackID:
        :param func:
        connect:用来决定是否将变量全局化
        str_btn:用来输入btn的名称
        """
        widget = UIItemsMake.make_widget()
        layout = QVBoxLayout(widget)
        label = self.ui.Label = UIItemsMake.make_label(labelText, font_size='16pt')
        setattr(self.ui, widgetName, widget)
        setattr(self.ui, widgetLayout, layout)
        setattr(self.ui, widgetName + 'Label', label)
        layout.addWidget(label)
        scrollArea = QScrollArea()
        scrollArea.setWidgetResizable(True)
        container = QWidget()
        containerLayout = QVBoxLayout()  # 创建布局

        groups.clear()

        for sp in self.Box.residue_click:
            groupBox = UIItemsMake.make_group_box(sp)
            groupBoxLayout = QVBoxLayout(groupBox)
            for atom in self.Box.residues_atoms[sp]:
                btn = UIItemsMake.make_radio_check(radioOrCheck, atom)
                groupBoxLayout.addWidget(btn)
            groups.append(groupBox)
            containerLayout.addWidget(groupBox)
        container.setLayout(containerLayout)  # 为container设置布局
        scrollArea.setWidget(container)  # 将container设置为scrollArea的子部件

        layout.addWidget(scrollArea)
        btnNext = UIItemsMake.make_btn(btnName)
        layout.addWidget(btnNext)
        self.ui.stackedWidget_Analysis.addWidget(widget)
        self.ui.stackedWidget_Analysis.setCurrentIndex(stackID)
        btnNext.clicked.connect(func)
        if connect: setattr(self.ui, 'btnAnalysisRun', btnNext)

    def ListLayout(self, widgetName: str, widgetLayout: str, labelText: str,
                   groups, listStr, radioOrCheck, btnName: str, stackID, func, connect=False):
        widget = UIItemsMake.make_widget()
        layout = QVBoxLayout(widget)
        label = UIItemsMake.make_label(labelText)
        setattr(self.ui, widgetName, widget)
        setattr(self.ui, widgetLayout, layout)
        setattr(self.ui, widgetName + 'Label', label)
        layout.addWidget(label)
        groups.clear()
        for part in listStr:
            btn = UIItemsMake.make_radio_check(radioOrCheck, part)
            layout.addWidget(btn)
            groups.append(btn)
        btnNext = UIItemsMake.make_btn(btnName)
        layout.addWidget(btnNext)
        self.ui.stackedWidget_Analysis.addWidget(widget)
        self.ui.stackedWidget_Analysis.setCurrentIndex(stackID)
        btnNext.clicked.connect(func)
        if connect: setattr(self.ui, 'btnAnalysisRun', btnNext)

    def SpinLayout(self, widgetName: str
                   , widgetLayout: str
                   , labelText
                   , spin_value
                   , spin_min
                   , spin_max
                   , btnName: str
                   , stackID
                   , func
                   , connect=False
                   , num=1):
        widget = UIItemsMake.make_widget()
        layout = QVBoxLayout(widget)
        setattr(self.ui, widgetName, widget)
        setattr(self.ui, widgetLayout, layout)
        if num != 1:
            for i in range(num):
                label = UIItemsMake.make_label(labelText[i])
                spin_box = UIItemsMake.make_spin_box(spin_value[i], spin_min[i], spin_max[i])
                setattr(self.ui, widgetName + 'Label' + str(i), label)
                setattr(self.ui, widgetName + 'SpinBox' + str(i), spin_box)
                layout.addWidget(label)
                layout.addWidget(spin_box)
        else:
            label = UIItemsMake.make_label(labelText)
            spin_box = UIItemsMake.make_spin_box(spin_value, spin_min, spin_max)
            setattr(self.ui, widgetName + 'Label', label)
            setattr(self.ui, widgetName + 'SpinBox', spin_box)
            layout.addWidget(label)
            layout.addWidget(spin_box)
        btnNext = UIItemsMake.make_btn(btnName)
        layout.addWidget(btnNext)
        self.ui.stackedWidget_Analysis.addWidget(widget)
        self.ui.stackedWidget_Analysis.setCurrentIndex(stackID)
        btnNext.clicked.connect(func)
        if connect: setattr(self.ui, 'btnAnalysisRun', btnNext)

    def makeFigure(self):
        # window = AnalysisFigureWidget(residues=self.Box.residueClick, runMethod=self.Info.runMethod,
        #                               resultPath=self.Info.path_result)
        self.window = AnalysisFigureWidget(self.ui
                                           , residues=self.Box.residue_click
                                           , runMethod=self.Info.method_analysis
                                           , resultPath=self.Info.path_result)
        self.window.show()


@dataclass
class InfoAnalysis:
    ui: object
    method_analysis: str = field(init=False)
    path_structure: str = field(init=False)
    path_trajectory: str = field(init=False)
    frame_first: int = field(init=False, default=0)
    frame_last: int = field(init=False, default=-1)
    step: int = field(init=False, default=1)
    K: int = field(init=False, default=21)
    path_result: str = field(init=False, default=None)

    TEST_PATH_RESULT: str = 'E:/excel/temp.xlsx'

    def __post_init__(self):
        self.method_analysis = self.ui.comboBoxMethod.currentText()
        self.path_structure = self.ui.editStructure.text()
        self.path_trajectory = self.ui.editTrajectory.text()

    def get_text(self):
        # Frame
        self.frame_first = self.ui.editFirstFrame.value()
        self.frame_last = self.ui.editLastFrame.value()
        self.step = self.ui.editStep.value()
        # Path
        self.path_result = self.ui.editResult.text() or self.TEST_PATH_RESULT
        # k
        self.K = self.ui.editK.value()


class SelBox:
    """储存选择的残基及原子信息"""
    __slots__ = ('u'
                 , 'residue_click'
                 , 'residues'
                 , 'residues_atoms'
                 , 'residues_list'
                 , '_configs')

    def __init__(self):
        self.u = None
        self.residues = None  # 记录结构文件中包含的全部残基名称
        self.residues_atoms = {}  # 储存全部残基名称及对应的原子名称
        self.residues_list = []  # 储存残基选择的按钮，用于后面筛选出残基
        self.residue_click = []  # 储存选择了的残基名称
        self._configs = defaultdict(lambda: None)

    def get_config(self, config_name):
        if config_name not in self._configs:
            self._configs[config_name] = ConfigFactory.create(config_name)
        return self._configs[config_name]

    def __getattr__(self, item):
        return self.get_config(item)


class HeightLayout(AnalysisLayout):

    def step_1(self):
        """得到选取的头部原子"""
        self.AtomsLayout('widgetAtomsHeight'
                         , 'VLayoutAtomsHeight'
                         , 'select head atom'
                         , self.Box.get_config('Height').HeightHeadGroups
                         , QRadioButton
                         , 'Next!'
                         , 1
                         , lambda: self.step_2())

    def step_2(self):
        """得到选取的尾部原子"""
        AnalysisUtils.get_atom(self.Box.get_config('Height').HeightHeadGroups
                               , QRadioButton
                               , self.Box.get_config('Height').HeightHeadAtoms)
        self.AtomsLayout(
                         'widgetTailAtomsHeight'
                         , 'VLayoutTailAtomsHeight'
                         , 'select tail atom(s)'
                         , self.Box.get_config('Height').HeightTailGroups
                         , QCheckBox
                         , 'Run!'
                         , 2
                         , func=self.run
                         , connect=True
                         )

    @AnalysisLayout._addProgressBar
    def run(self):

        self.Info.get_text()  # 得到选择的参数，例如帧数和K
        AnalysisUtils.get_atom(self.Box.get_config('Height').HeightTailGroups
                               , QCheckBox
                               , self.Box.get_config('Height').HeightTailAtoms)
        dictHeadTail = AnalysisUtils.get_key_value(self.Box.get_config('Height').HeightHeadAtoms
                                                   , self.Box.get_config('Height').HeightTailAtoms)

        print('head'
              , self.Box.get_config('Height').HeightHeadAtoms
              , 'tail'
              , self.Box.HeightTailAtoms)

        self.cls = Height(self.Box.u
                          , dictHeadTail
                          , file_path=self.Info.path_result
                          , k=self.Info.K)


class SZLayout(AnalysisLayout):
    FF_TYPE = None
    CHAIN = None

    @staticmethod
    def getInfo(boxes):
        for box in boxes:
            if box.isChecked():
                return box.text()

    def setValueFFtype(self, boxes):
        self.FF_TYPE = self.getInfo(boxes)

    def setValueChain(self, boxes):
        self.CHAIN = self.getInfo(boxes)

    def step_1(self):
        self.ListLayout('widgetFFSZ'
                        , 'VLayoutFFSZ'
                        , 'select force field '
                        , self.Box.get_config('SZ').SZFFGroups
                        , ['All-Atom', 'United-Atom', 'Coarse-Atom']
                        , QRadioButton
                        , 'Next!'
                        , 1
                        , lambda: self.step_2())

    def step_2(self):
        self.setValueFFtype(self.Box.get_config('SZ').SZFFGroups)
        self.AtomsLayout('widgetHeadSZ'
                         , 'VLayoutHeadSZ'
                         , 'select head atom to fit'
                         , self.Box.get_config('SZ').SZHeadGroups
                         , QRadioButton
                         , 'Next!'
                         , 2
                         , lambda: self.step_3())

    def step_3(self):
        AnalysisUtils.get_atom(self.Box.get_config('SZ').SZHeadGroups
                               , QRadioButton
                               , self.Box.get_config('SZ').SZHeadAtoms)
        self.ListLayout('widgetChainSZ'
                        , 'VLayoutChainSZ'
                        , 'select chain to analyse'
                        , self.Box.get_config('SZ').SZChainGroups
                        , ['sn1', 'sn2', 'sn1 and sn2']
                        , QRadioButton
                        , 'Run!'
                        , 3
                        , self.run
                        , connect=True
                        )

    @AnalysisLayout._addProgressBar
    def run(self):
        self.Info.get_text()
        self.setValueChain(self.Box.get_config('SZ').SZChainGroups)
        if self.FF_TYPE == 'Coarse-Atom':
            self.cls = SZ(self.Box.u
                          , self.Box.get_config('SZ').SZHeadAtoms
                          , self.CHAIN
                          , path=self.Info.path_result
                          , k=self.Info.K)
        else:
            print('暂时不支持%s力场' % self.FF_TYPE)


class AreaLayout(AnalysisLayout):
    def step_1(self):
        self.AtomsLayout('widgetAtomsArea'
                         , 'VLayoutAtomsArea'
                         , 'select head atom'
                         , self.Box.get_config('Area').AreaHeadGroups
                         , QRadioButton
                         , 'Run!'
                         , 1
                         , self.run
                         , connect=True
                         )

    @AnalysisLayout._addProgressBar
    def run(self):
        self.Info.get_text()
        AnalysisUtils.get_atom(self.Box.get_config('Area').AreaHeadGroups
                               , QRadioButton
                               , self.Box.get_config('Area').AreaHeadAtoms
                               )
        self.cls = Area(universe=self.Box.u
                        , residueGroup=self.Box.get_config('Area').AreaHeadAtoms
                        , file_path=self.Info.path_result
                        , k=self.Info.K)


class MeanCurvatureLayout(AnalysisLayout):
    def step_1(self):
        self.AtomsLayout('widgetAtomsMeanCurvature'
                         , 'VLayoutAtomsMeanCurvature'
                         , 'select head atom'
                         , self.Box.get_config('MeanCurvature').MCHeadGroups
                         , QRadioButton
                         , 'Run!'
                         , 1
                         , self.run
                         , connect=True)

    @AnalysisLayout._addProgressBar
    def run(self):
        self.Info.get_text()
        AnalysisUtils.get_atom(self.Box.get_config('MeanCurvature').MCHeadGroups
                               , QRadioButton
                               , self.Box.get_config('MeanCurvature').MeanCurvatureHeadAtoms
                               )
        self.cls = Curvature(self.Box.u
                             , self.Box.get_config('MeanCurvature').MeanCurvatureHeadAtoms
                             , path=self.Info.path_result
                             , k=self.Info.K
                             , method='mean')


class AnisotropyLayout(AnalysisLayout):
    def step_1(self):
        self.AtomsLayout('widgetAtomsAnisotropy'
                         , 'VLayoutAtomsAnisotropy'
                         , 'select head atom'
                         , self.Box.get_config('Anisotropy').AnsHeadGroups
                         , QRadioButton
                         , 'Run!'
                         , 1
                         , self.run
                         , connect=True)

    @AnalysisLayout._addProgressBar
    def run(self):
        self.Info.get_text()
        AnalysisUtils.get_atom(self.Box.get_config('Anisotropy').AnsHeadGroups
                               , QRadioButton
                               , self.Box.get_config('Anisotropy').AnisotropyHeadAtoms)
        self.cls = Anisotropy(self.Box.u
                              , self.Box.get_config('Anisotropy').AnisotropyHeadAtoms
                              , filePath=self.Info.path_result)


class PCALayout(AnalysisLayout):
    def step_1(self):
        self.AtomsLayout('widgetAtomsPCA'
                         , 'VLayoutAtomsPCA'
                         , 'select head atom'
                         , self.Box.get_config('PCA').PCAHeadGroups
                         , QRadioButton
                         , 'Run!'
                         , 1
                         , self.run
                         , connect=True)

    @AnalysisLayout._addProgressBar
    def run(self):
        self.Info.get_text()
        AnalysisUtils.get_atom(self.Box.get_config('PCA').PCAHeadGroups
                               , QRadioButton
                               , self.Box.get_config('PCA').PCAHeadAtoms)
        self.cls = PCA(self.Box.u
                              , self.Box.get_config('PCA').PCAHeadAtoms
                              , filePath=self.Info.path_result)


class RDLayout(AnalysisLayout):
    def step_1(self):
        self.AtomsLayout('widgetAtomsRD'
                         , 'VLayoutAtomsRD'
                         , 'select head atom'
                         , self.Box.get_config('RadialDistribution').RDHeadGroups
                         , QRadioButton
                         , 'Run!'
                         , 1
                         , self.run
                         , connect=True)

    @AnalysisLayout._addProgressBar
    def run(self):
        self.Info.get_text()
        AnalysisUtils.get_atom(self.Box.get_config('RadialDistribution').RDHeadGroups
                               , QRadioButton
                               , self.Box.get_config('RadialDistribution').RDHeadAtoms)
        self.cls = CalRad(self.Box.u
                          , self.Box.get_config('RadialDistribution').RDHeadAtoms
                          , filePath=self.Info.path_result)


class PressureLayout(AnalysisLayout):
    def step_1(self):
        self.AtomsLayout('widgetBallResiduesPressure'
                         , 'VLayoutBallResiduesPressure'
                         , 'select head atom'
                         , self.Box.get_config('Pressure').PressureHeadGroups
                         , QRadioButton
                         , 'Next!'
                         , 1
                         , func=self.step_2
                         , connect=True)

    def step_2(self):
        AnalysisUtils.get_atom(self.Box.get_config('Pressure').PressureHeadGroups
                               , QRadioButton
                               , self.Box.get_config('Pressure').PressureHeadAtoms)
        self.AtomsLayout('widgetGasResidues'
                         , 'VLayoutGasResidues'
                         , 'select gas residue'
                         , self.Box.get_config('Pressure').GasGroups
                         , QRadioButton
                         , 'Next!'
                         , 2
                         , func=self.step_3
                         , connect=True)

    def step_3(self):

        self.AtomsLayout('widgetGasAtoms'
                         , 'VLayoutGasAtoms'
                         , 'select gas atom'
                         , self.Box.get_config('Pressure').GasHeadAtoms
                         , QRadioButton
                         , 'Run!'
                         , 3
                         , self.step_4
                         , connect=True)

    def step_4(self):
        AnalysisUtils.get_atom(self.Box.get_config('Pressure').GasGroups
                               , QRadioButton
                               , self.Box.get_config('Pressure').GasHeadAtoms)
        self.SpinLayout('widgetNcircles'
                        , 'VLayoutNcircles'
                        , 'select num of circle'
                        , 10
                        , 0
                        , 1000000
                        , 'RUN!'
                        , 4
                        , self.run
                        , connect=True)

    @AnalysisLayout._addProgressBar
    def run(self):
        self.Info.get_text()
        self.cls = Pressure(self.Box.u
                            , self.Box.get_config('Pressure').PressureHeadAtoms
                            , self.Box.get_config('Pressure').GasHeadAtoms
                            , n_circles=AnalysisUtils.get_spin_value(self.ui.widgetNcirclesSpinBox)
                            , filePath=self.Info.path_result)


class GyrationLayout(AnalysisLayout):
    def step_1(self):
        self.AtomsLayout('widgetAtomsGyration'
                         , 'VLayoutAtomsGyration'
                         , 'select head atoms'
                         , self.Box.get_config('Gyration').GyrHeadGroups
                         , QCheckBox
                         , 'Run!'
                         , 1
                         , self.run
                         , connect=True)

    @AnalysisLayout._addProgressBar
    def run(self):
        self.Info.get_text()
        AnalysisUtils.get_atom(self.Box.get_config('Gyration').GyrHeadGroups
                               , QCheckBox
                               , self.Box.get_config('Gyration').GyrationHeadAtoms)
        self.cls = Gyration_py(self.Box.u
                               , self.Box.get_config('Gyration').GyrationHeadAtoms
                               , filePath=self.Info.path_result)


class ClusterLayout(AnalysisLayout):
    def step_1(self):
        self.AtomsLayout('widgetAtomsCluster'
                         , 'VLayoutAtomsCluster'
                         , 'select atoms'
                         , self.Box.get_config('Cluster').CLGroups
                         , QCheckBox
                         , 'Next'
                         , 1
                         , self.step_2)

    def step_2(self):
        self.SpinLayout('widgetCutoff'
                        , 'VLayoutCutoff'
                        , 'select cutoff value(A)'
                        , 12
                        , 0
                        , 1000000
                        , 'RUN!'
                        , 2
                        , self.run
                        , connect=True)

    @AnalysisLayout._addProgressBar
    def run(self):
        self.Info.get_text()
        AnalysisUtils.get_atom(self.Box.get_config('Cluster').CLGroups
                               , QCheckBox
                               , self.Box.get_config('Cluster').CLHeadAtoms)
        self.Box.get_config('Cluster').CLCutoff = AnalysisUtils.get_spin_value(self.ui.widgetCutoffSpinBox)
        self.cls = Cluster(self.Box.u
                           , self.Box.get_config('Cluster').CLHeadAtoms
                           , file_path=self.Info.path_result
                           , cutoff=self.Box.get_config('Cluster').CLCutoff)


class NClusterLayout(AnalysisLayout):
    def step_1(self):
        self.AtomsLayout('widgetAtomsNCluster'
                         , 'VLayoutAtomsNCluster'
                         , 'select atoms'
                         , self.Box.get_config('NCluster').NCLGroups
                         , QCheckBox
                         , 'Next'
                         , 1
                         , self.step_2)

    def step_2(self):
        self.SpinLayout('widgetNCutoff'
                        , 'VLayoutNCutoff'
                        , ['select cutoff value(A)', 'select cutoff number']
                        , [12, 10]
                        , [0, 0]
                        , [1000000, 1000000]
                        , 'RUN!'
                        , 2
                        , self.run
                        , connect=True
                        , num=2)

    @AnalysisLayout._addProgressBar
    def run(self):
        self.Info.get_text()
        AnalysisUtils.get_atom(self.Box.get_config('NCluster').NCLGroups
                               , QCheckBox
                               , self.Box.get_config('NCluster').NCLHeadAtoms)
        self.Box.get_config('NCluster').NCLCutoff = AnalysisUtils.get_spin_value(self.ui.widgetNCutoffSpinBox0)
        self.Box.get_config('NCluster').NCutoff = AnalysisUtils.get_spin_value(self.ui.widgetNCutoffSpinBox1)
        self.cls = NCluster(self.Box.u
                           , self.Box.get_config('NCluster').NCLHeadAtoms
                           , file_path=self.Info.path_result
                           , cutoff=self.Box.get_config('NCluster').NCLCutoff
                           , N_cutoff=self.Box.get_config('NCluster').NCutoff)


class BaseConfig:
    __slots__ = ()

    def __repr__(self):
        return f'<{self.__class__.__name__}>'


class HeightConfig(BaseConfig):
    __slots__ = ('HeightHeadAtoms', 'HeightTailAtoms', 'HeightHeadGroups', 'HeightTailGroups')

    def __init__(self):
        self.HeightHeadAtoms = {}
        self.HeightTailAtoms = {}
        self.HeightHeadGroups = []
        self.HeightTailGroups = []


class SZConfig(BaseConfig):
    __slots__ = ('SZHeadAtoms', 'SZFFGroups', 'SZHeadGroups', 'SZChainGroups')

    def __init__(self):
        self.SZHeadAtoms = {}
        self.SZFFGroups = []
        self.SZHeadGroups = []
        self.SZChainGroups = []


class AreaConfig(BaseConfig):
    __slots__ = ('AreaHeadAtoms', 'AreaHeadGroups')

    def __init__(self):
        self.AreaHeadAtoms = {}
        self.AreaHeadGroups = []


class MCConfig(BaseConfig):
    __slots__ = ('MeanCurvatureHeadAtoms', 'MCHeadGroups')

    def __init__(self):
        self.MeanCurvatureHeadAtoms = {}
        self.MCHeadGroups = []


class ASPConfig(BaseConfig):
    __slots__ = ('AnisotropyHeadAtoms', 'AnsHeadGroups')

    def __init__(self):
        self.AnisotropyHeadAtoms = {}
        self.AnsHeadGroups = []


class PCAConfig(BaseConfig):
    __slots__ = ('PCAHeadAtoms', 'PCAHeadGroups')

    def __init__(self):
        self.PCAHeadAtoms = {}
        self.PCAHeadGroups = []


class RDConfig(BaseConfig):
    __slots__ = ('RDHeadAtoms', 'RDHeadGroups')

    def __init__(self):
        self.RDHeadAtoms = {}
        self.RDHeadGroups = []


class PressureConfig(BaseConfig):
    __slots__ = ('n_circles', 'PressureHeadAtoms', 'GasGroups', 'PressureHeadGroups', 'GasHeadAtoms')

    def __init__(self):
        self.PressureHeadAtoms = {}
        self.PressureHeadGroups = []
        self.GasGroups = []
        self.GasHeadAtoms = {}
        self.n_circles = 10


class GRConfig(BaseConfig):
    __slots__ = ('GyrationHeadAtoms', 'GyrHeadGroups')

    def __init__(self):
        self.GyrationHeadAtoms = {}
        self.GyrHeadGroups = []


class CLConfig(BaseConfig):
    __slots__ = ('CLHeadAtoms', 'CLGroups', 'CLCutoff')

    def __init__(self):
        self.CLHeadAtoms = {}
        self.CLGroups = []
        self.CLCutoff = 12

class NCLConfig(BaseConfig):
    __slots__ = ('NCLHeadAtoms', 'NCLGroups', 'NCLCutoff', 'NCutoff')

    def __init__(self):
        self.NCLHeadAtoms = {}
        self.NCLGroups = []
        self.NCLCutoff = 12
        self.NCutoff = 10


class ConfigFactory:
    __slots__ = ()
    CONFIG_MAP = {
        'Height': HeightConfig
        , 'SZ': SZConfig
        , 'Area': AreaConfig
        , 'MeanCurvature': MCConfig
        , 'Anisotropy': ASPConfig
        , 'RadialDistribution': RDConfig
        , 'Pressure': PressureConfig
        , 'Gyration': GRConfig
        , 'Cluster': CLConfig
        , 'NCluster': NCLConfig
        , 'PCA': PCAConfig
    }

    @classmethod
    def create(cls, config_name):
        config_class = cls.CONFIG_MAP.get(config_name, None)
        if config_class:
            return config_class()


from .Analysis_Figure import *
from functools import partial


class AnalysisFigureWidget(QWidget):
    TYPE_ID = {'Height': 0, 'SZ': 0, 'MeanCurvature': 0, 'Area': 0, 'Anisotropy': 1, 'RadialDistribution': 2,
               'Gyration': 1, 'Pressure': 2}
    FIGURE_TYPE = {0: ['Bar', 'Line', 'Scatter', 'Map'], 1: ['Bar', 'Line'], 2: ['Line']}

    def __init__(self
                 , ui
                 , residues
                 , runMethod
                 , resultPath):

        super().__init__()
        self.ui_analysis = Ui_Form()
        self.ui_analysis.setupUi(self)
        self.ui = ui
        self.residues = residues  # 获取得到所有的残基名称

        self.runMethod = runMethod  # 获取得到当前使用的分析方法

        self.resultPath = resultPath  # 获取得到结果的保存路径

        # 参数设置
        self.ui_analysis.FigureColorLayout = None
        self.ui_analysis.FigureShapeWidget = None
        self.btnLanguageClick()

        # 储存信息
        self.LineInfo = defaultdict(lambda: None)
        self.BarInfo = defaultdict(lambda: None)
        self.ScaInfo = defaultdict(lambda: None)

        self.ColorInfo = defaultdict(lambda: None)
        self.ShapeInfo = []

        self.information, self.results = read_excel(self.resultPath)
        self.description = self.information
        self.lipids_type = self.results['Resname'].unique() if 'Resname' in self.results.columns else None
        # 函数绑定
        # Line
        self.ui_analysis.figure_line_btn_color_2.clicked.connect(self.analysisBtnColor)
        self.ui_analysis.figure_line_btn_color_2.clicked.connect(self.openCloseAnalysisColorBox)
        # Bar
        self.ui_analysis.figure_bar_btn_color_2.clicked.connect(self.analysisBtnColor)
        self.ui_analysis.figure_bar_btn_color_2.clicked.connect(self.openCloseAnalysisColorBox)
        self.ui_analysis.figure_bar_btn_trend_2.clicked.connect(partial(self.btnColorClicked,
                                                                        self.ui_analysis.figure_bar_btn_trend_2,
                                                                        self.BarInfo['trend_color']))
        # Scatter
        self.ui_analysis.figure_scatter_btn_shape_2.clicked.connect(self.analysisBtnShape)
        self.ui_analysis.figure_scatter_btn_shape_2.clicked.connect(self.openCloseAnalysisShapeBox)
        # Tab
        self.ui_analysis.tabWidget.currentChanged.connect(self.openCloseAnalysisExtra)
        # Language
        self.ui.btn_language.clicked.connect(self.btnLanguageClick)
        # Make Figure
        self.ui_analysis.btn_figure_run.clicked.connect(self.analysisBtnMakeFigure)

    def analysisBtnColor(self):
        if self.TYPE_ID[self.runMethod] == 0 or self.TYPE_ID[self.runMethod] == 2:  # 如果是Lipids或者径向函数分布
            if not getattr(self.ui_analysis, 'FigureColorLayout'):
                self.ui_analysis.figure_color_extra_box.setStyleSheet(u"font: 15pt \"\u534e\u6587\u7ec6\u9ed1\";")
                self.ui_analysis.FigureColorLayout = QVBoxLayout()
                for id, residue in enumerate(self.residues):
                    btn = QPushButton(residue)
                    btn.clicked.connect(partial(self.cellResidueClicked, id, btn))
                    self.ui_analysis.FigureColorLayout.addWidget(btn)
                self.ui_analysis.figure_color_extra_box.setLayout(self.ui_analysis.FigureColorLayout)

        elif self.TYPE_ID[self.runMethod] == 1:  # 如果是Bubble
            if self.ui_analysis.tabWidget.currentIndex() == 0:
                self.btnColorClicked(self.ui_analysis.figure_line_btn_color_2, self.LineInfo['bubble_color'])
            elif self.ui_analysis.tabWidget.currentIndex() == 1:
                self.btnColorClicked(self.ui_analysis.figure_bar_btn_color_2, self.BarInfo['bubble_color'])

    def analysisBtnShape(self):
        if not getattr(self.ui_analysis, 'FigureShapeWidget'):
            self.ui_analysis.FigureShapeLayout = QVBoxLayout(self.ui_analysis.figure_shape_extra_box)
            label_shape = UIItemsMake.make_label('Shape', color='rgb(33, 37, 43)')
            self.ui_analysis.FigureShapeLayout.addWidget(label_shape)
            scrollArea = QScrollArea()
            scrollArea.setWidgetResizable(True)
            self.ui_analysis.FigureShapeWidget = QWidget()
            containerLayout = QVBoxLayout()
            for sp in self.residues:
                groupBox = UIItemsMake.make_group_box(sp, title_color='rgb(33, 37, 43)')
                groupLayout = QVBoxLayout(groupBox)
                for sh in ['o', 'p', 's', '^', '*', 'x', '+']:
                    radio = UIItemsMake.make_radio_check(QRadioButton, sh, color='rgb(33, 37, 43)')
                    groupLayout.addWidget(radio)
                self.shapeInfo.append(groupBox)
                containerLayout.addWidget(groupBox)
            self.ui_analysis.FigureShapeWidget.setLayout(containerLayout)
            scrollArea.setWidget(self.ui_analysis.FigureShapeWidget)
            self.ui_analysis.FigureShapeLayout.addWidget(scrollArea)

    def analysisBtnMakeFigure(self):
        method = self.ui_analysis.tabWidget.currentIndex()
        if method == 0:  # Line
            FigureLine(self.description
                       , self.results
                       , self.getLine()).plot()
        elif method == 1:  # Bar
            FigureBar(self.description
                      , self.results
                      , self.getBar()).plot()
        elif method == 2:  # Scatter
            FigureScatter(self.description
                          , self.results
                          , self.getScatter()).plot()

    def btnColorClicked(self, btn, info):
        color = QColor()
        new_color = QColorDialog.getColor(color)
        if new_color.isValid():
            btn.setStyleSheet(f'background-color: rgb({new_color.red()}, {new_color.green()}, {new_color.blue()});')
            info.append((new_color.red() / 255, new_color.green() / 255, new_color.blue() / 255))

    def cellResidueClicked(self, id, btn):
        color = QColor()
        # 打开颜色选择对话框并获取新颜色
        new_color = QColorDialog.getColor(color)
        if new_color.isValid():
            # 将新颜色的RGB值存储到字典中
            self.ColorInfo[self.residues[id]] = (new_color.red() / 255
                                                 , new_color.green() / 255
                                                 , new_color.blue() / 255
                                                 )
            # 设置单元格的背景颜色为新颜色的RGB值，忽略Alpha通道
            btn.setStyleSheet(f'background-color: rgb({new_color.red()}, {new_color.green()}, {new_color.blue()});')

    def getLine(self):
        self.LineInfo.update({
            'axis_scale': self.ui_analysis.figure_line_spin_axis_scale_2.value()
            , 'axis_text': self.ui_analysis.figure_line_spin_axis_text_size_2.value()
            , 'grid_size': self.ui_analysis.figure_line_spin_grid_size_2.value()
            , 'x_title': self.ui_analysis.figure_line_edit_x_2.text() or "Frames"
            , 'y_title': self.ui_analysis.figure_line_edit_y_2.text() or self.description
            , 'x_min': self.ui_analysis.figure_line_spin_x_min_2.value()
            , 'x_max': self.ui_analysis.figure_line_spin_x_max_2.value()
            , 'y_min': self.ui_analysis.figure_line_spin_y_min_2.value()
            , 'y_max': self.ui_analysis.figure_line_spin_y_max_2.value()
            , 'marker_size': self.ui_analysis.figure_line_spin_marker_size_2.value()
            , 'marker_shape': self.ui_analysis.figure_line_como_marker.currentText()
            , 'color': self.ColorInfo
        })

        return self.LineInfo

    def getBar(self):

        self.BarInfo.update({
            'axis_scale': self.ui_analysis.figure_bar_spin_axis_scale_2.value()
            , 'axis_text': self.ui_analysis.figure_bar_spin_axis_text_size_2.value()
            , 'x_title': self.ui_analysis.figure_bar_edit_x_2.text()
            , 'y_title': self.ui_analysis.figure_bar_edit_y_2.text() or self.description
            , 'y_min': self.ui_analysis.figure_bar_spin_y_min_2.value()
            , 'y_max': self.ui_analysis.figure_bar_spin_y_max_2.value()
            , 'trend_size': self.ui_analysis.figure_bar_spin_trend_2.value()
            , 'up_bar_value': self.ui_analysis.figure_bar_spin_bar_2.value()
            , 'error_deci': self.ui_analysis.figure_bar_radio_error_2.isChecked()
            , 'color': self.ColorInfo
        })

        return self.BarInfo

    def getScatter(self):
        self.ScaInfo.update({
            'grid_size': self.ui_analysis.figure_scatter_spin_grid_size_2.value()
            , 'bar_min': self.ui_analysis.figure_scatter_color_min_2.value()
            , 'bar_max': self.ui_analysis.figure_scatter_color_max_2.value()
            , 'bar_color': self.ui_analysis.figure_scatter_como_color_2.currentText()
            , 'shape_size': self.ui_analysis.figure_scatter_spin_shape_size_2.value()
            , 'shape': {}
        })
        for group in self.ShapeInfo:
            for radio in group.findChildren(QRadioButton):
                if radio.isChecked():
                    self.ScaInfo['shape'][group.title()] = radio.text()
        return self.ScaInfo

    def btnLanguageClick(self):
            if self.ui.btn_language.text() == '中' and self.ui_analysis.figure_line_label_axis_title == 'Axis Tick Size':
                font_style = "font: 16pt '华文细黑';"
                # Figure
                self.ui_analysis.btn_figure_run.setText('开始绘图！')
                ## Liui_analysis.ne
                self.ui_analysis.figure_line_label_axis_tick.setText('坐标刻度大小')
                self.ui_analysis.figure_line_label_axis_tick.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_axis_title.setText('标题大小')
                self.ui_analysis.figure_line_label_axis_title.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_legend.setText('图例大小')
                self.ui_analysis.figure_line_label_legend.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_x.setText('X轴标题')
                self.ui_analysis.figure_line_label_x.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_y.setText('Y轴标题')
                self.ui_analysis.figure_line_label_y.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_x_range.setText('X轴显示范围')
                self.ui_analysis.figure_line_label_x_range.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_y_range.setText('Y轴显示范围')
                self.ui_analysis.figure_line_label_y_range.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_marker.setText('标记点大小')
                self.ui_analysis.figure_line_label_marker.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_color.setText('颜色设置')
                self.ui_analysis.figure_line_label_color.setStyleSheet(font_style)
                self.ui_analysis.figure_line_btn_color_2.setText('选择颜色')
                ## Baui_analysis.r
                self.ui_analysis.figure_bar_label_axis_tick.setText('坐标刻度大小')
                self.ui_analysis.figure_bar_label_axis_tick.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_axis_title.setText('标题大小')
                self.ui_analysis.figure_bar_label_axis_title.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_x.setText('X轴标题')
                self.ui_analysis.figure_bar_label_x.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_y.setText('Y轴标题')
                self.ui_analysis.figure_bar_label_y.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_y_range.setText('Y轴显示范围')
                self.ui_analysis.figure_bar_label_y_range.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_trend.setText('趋势线设置')
                self.ui_analysis.figure_bar_label_trend.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_bar.setText('柱图上数值设置')
                self.ui_analysis.figure_bar_label_bar.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_color.setText('柱图颜色设置')
                self.ui_analysis.figure_bar_label_color.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_error.setText('误差棒设置')
                self.ui_analysis.figure_bar_label_error.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_btn_trend_2.setText('选择颜色')
                self.ui_analysis.figure_bar_btn_color_2.setText('选择颜色')
                ## Scui_analysis.atter
                self.ui_analysis.figure_scatter_label_legend.setText('图例大小')
                self.ui_analysis.figure_scatter_label_legend.setStyleSheet(font_style)
                self.ui_analysis.figure_scatter_label_range.setText('数值显示范围')
                self.ui_analysis.figure_scatter_label_range.setStyleSheet(font_style)
                self.ui_analysis.figure_scatter_label_color.setText('颜色类型设置')
                self.ui_analysis.figure_scatter_label_color.setStyleSheet(font_style)
                self.ui_analysis.figure_scatter_label_shape.setText('形状设置')
                self.ui_analysis.figure_scatter_label_shape.setStyleSheet(font_style)
                self.ui_analysis.figure_scatter_label_shape_size.setText('形状大小')
                self.ui_analysis.figure_scatter_label_shape_size.setStyleSheet(font_style)

            elif self.ui.btn_language.text() == "ABC" and self.ui_analysis.figure_line_label_axis_title == '坐标刻度大小':
                font_style = "font: 16pt '华文细黑';"
                # Figure
                self.ui_analysis.btn_figure_run.setText('Run！')
                ## Line
                self.ui_analysis.figure_line_label_axis_tick.setText('Axis Tick Size')
                self.ui_analysis.figure_line_label_axis_tick.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_axis_title.setText('Axis Title Size')
                self.ui_analysis.figure_line_label_axis_title.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_legend.setText('Legend Size')
                self.ui_analysis.figure_line_label_legend.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_x.setText('X-Title')
                self.ui_analysis.figure_line_label_x.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_y.setText('Y-Title')
                self.ui_analysis.figure_line_label_y.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_x_range.setText('X-Range')
                self.ui_analysis.figure_line_label_x_range.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_y_range.setText('Y-Range')
                self.ui_analysis.figure_line_label_y_range.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_marker.setText('Marker Size')
                self.ui_analysis.figure_line_label_marker.setStyleSheet(font_style)
                self.ui_analysis.figure_line_label_color.setText('Color')
                self.ui_analysis.figure_line_label_color.setStyleSheet(font_style)
                self.ui_analysis.figure_line_btn_color_2.setText('Select Color')
                self.ui_analysis.figure_line_btn_color_2.setText('Select Color')
                ## Bar
                self.ui_analysis.figure_bar_label_axis_tick.setText('Axis Tick Size')
                self.ui_analysis.figure_bar_label_axis_tick.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_axis_title.setText('Axis Title Size')
                self.ui_analysis.figure_bar_label_axis_title.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_x.setText('X-Title')
                self.ui_analysis.figure_bar_label_x.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_y.setText('Y-Title')
                self.ui_analysis.figure_bar_label_y.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_y_range.setText('Y-Range')
                self.ui_analysis.figure_bar_label_y_range.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_trend.setText('Trend Line')
                self.ui_analysis.figure_bar_label_trend.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_bar.setText('Bar Value')
                self.ui_analysis.figure_bar_label_bar.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_color.setText('Color')
                self.ui_analysis.figure_bar_label_color.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_label_error.setText('Error Bar')
                self.ui_analysis.figure_bar_label_error.setStyleSheet(font_style)
                self.ui_analysis.figure_bar_btn_trend_2.setText('Select Color')
                self.ui_analysis.figure_bar_btn_color_2.setText('Select Color')
                ## Scatter
                self.ui_analysis.figure_scatter_label_legend.setText('Legend Size')
                self.ui_analysis.figure_scatter_label_legend.setStyleSheet(font_style)
                self.ui_analysis.figure_scatter_label_range.setText('Value Range')
                self.ui_analysis.figure_scatter_label_range.setStyleSheet(font_style)
                self.ui_analysis.figure_scatter_label_color.setText('Color Type')
                self.ui_analysis.figure_scatter_label_color.setStyleSheet(font_style)
                self.ui_analysis.figure_scatter_label_shape.setText('Shape')
                self.ui_analysis.figure_scatter_label_shape.setStyleSheet(font_style)
                self.ui_analysis.figure_scatter_label_shape_size.setText('Shape Size')
                self.ui_analysis.figure_scatter_label_shape_size.setStyleSheet(font_style)
                self.ui_analysis.figure_scatter_btn_shape_2.setText('Select Shape')

        # 底层的一些用法
    def openCloseAnalysisColorBox(self):
        if self.TYPE_ID[self.runMethod] == 0:
            self.toggleAnalysisColorBox(True)

    def openCloseAnalysisShapeBox(self):
        self.toggleAnalysisShapeBox(True)

    def openCloseAnalysisExtra(self, index):
        if index == 0 or index == 1:
            self.toggleAnalysisShapeBox(False)
        if index == 2:
            self.toggleAnalysisColorBox(False)

    def toggleAnalysisColorBox(self, enable):
        if enable:
            # GET WIDTH
            width = self.ui_analysis.figure_color_extra_box.width()
            maxExtend = 240
            standard = 0

            # SET MAX WIDTH
            if width == 0:
                widthExtended = maxExtend
            else:
                widthExtended = standard

            # ANIMATION
            self.Figure_colorBox = QPropertyAnimation(self.ui_analysis.figure_color_extra_box, b"minimumWidth")
            self.Figure_colorBox.setDuration(500)
            self.Figure_colorBox.setStartValue(width)
            self.Figure_colorBox.setEndValue(widthExtended)
            self.Figure_colorBox.setEasingCurve(QEasingCurve.InOutQuart)
            self.Figure_colorBox.start()
        else:
            if self.ui_analysis.figure_color_extra_box.width() != 0:
                self.Figure_colorBox = QPropertyAnimation(self.ui_analysis.figure_color_extra_box, b"minimumWidth")
                self.Figure_colorBox.setDuration(500)
                self.Figure_colorBox.setStartValue(self.ui_analysis.figure_color_extra_box.width())
                self.Figure_colorBox.setEndValue(0)
                self.Figure_colorBox.setEasingCurve(QEasingCurve.InOutQuart)
                self.Figure_colorBox.start()

    def toggleAnalysisShapeBox(self, enable):
        if enable:
            # GET WIDTH
            width = self.ui_analysis.figure_shape_extra_box.width()
            maxExtend = 240
            standard = 0

            # SET MAX WIDTH
            if width == 0:
                widthExtended = maxExtend
            else:
                widthExtended = standard

            # ANIMATION
            self.Figure_shapeBox = QPropertyAnimation(self.ui_analysis.figure_shape_extra_box, b"minimumWidth")
            self.Figure_shapeBox.setDuration(500)
            self.Figure_shapeBox.setStartValue(width)
            self.Figure_shapeBox.setEndValue(widthExtended)
            self.Figure_shapeBox.setEasingCurve(QEasingCurve.InOutQuart)
            self.Figure_shapeBox.start()
        else:
            if self.ui_analysis.figure_shape_extra_box.width() != 0:
                self.Figure_shapeBox = QPropertyAnimation(self.ui_analysis.figure_shape_extra_box, b"minimumWidth")
                self.Figure_shapeBox.setDuration(500)
                self.Figure_shapeBox.setStartValue(self.ui_analysis.figure_shape_extra_box.width())
                self.Figure_shapeBox.setEndValue(0)
                self.Figure_shapeBox.setEasingCurve(QEasingCurve.InOutQuart)
                self.Figure_shapeBox.start()

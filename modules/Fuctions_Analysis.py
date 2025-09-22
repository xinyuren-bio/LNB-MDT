import time
from collections import defaultdict
from threading import Thread
from dataclasses import dataclass, field

import numpy as np
from PySide6.QtCore import QObject, Signal, QPropertyAnimation, QEasingCurve, QTimer, QThread, QSize, Qt
from PySide6.QtGui import QFont, QColor, QIcon, QCursor
from PySide6.QtWidgets import QButtonGroup, QPushButton, QGroupBox, QCheckBox, QScrollArea, QWidget, QVBoxLayout, QHBoxLayout
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
        # 创建带图标的Back按钮
        print(f"*** Creating Back button ***")
        self.ui.btnBack = QPushButton()
        self.update_button_icons()  # 根据当前主题设置图标
        self.ui.btnBack.setIconSize(QSize(24, 24))
        self.ui.btnBack.setToolTip('Back')
        self.ui.btnBack.setCursor(QCursor(Qt.PointingHandCursor))
        
        # 添加测试连接
        def test_back_click():
            print(f"*** TEST: Back button clicked! ***")
            print(f"*** Current method: {getattr(self, 'method', 'No method')} ***")
        
        self.ui.btnBack.clicked.connect(test_back_click)
        self.ui.btnBack.clicked.connect(self.switchWidgetBack)
        print(f"*** Back button created and connected ***")
        
        # 创建带图标的Refresh按钮
        print(f"*** Creating Refresh button ***")
        self.ui.btnRefresh = QPushButton()
        self.update_button_icons()  # 根据当前主题设置图标
        self.ui.btnRefresh.setIconSize(QSize(24, 24))
        self.ui.btnRefresh.setToolTip('Refresh')
        self.ui.btnRefresh.setCursor(QCursor(Qt.PointingHandCursor))
        
        # 添加测试连接
        def test_refresh_click():
            print(f"*** TEST: Refresh button clicked! ***")
            print(f"*** Current method: {getattr(self, 'method', 'No method')} ***")
        
        self.ui.btnRefresh.clicked.connect(test_refresh_click)
        self.ui.btnRefresh.clicked.connect(self.refreshWidget)
        print(f"*** Refresh button created and connected ***")
        
        # 存储Handler实例，用于back/refresh功能
        self.density_handler = None
        self.density_multi_handler = None
        self.area_handler = None
        self.gyration_handler = None
        self.anisotropy_handler = None
        # 添加水平布局到顶部窗口
        HLayoutNextBack = QHBoxLayout(self.ui.widgetUpNextBack)
        HLayoutNextBack.addWidget(self.ui.btnBack)
        HLayoutNextBack.addWidget(self.ui.btnRefresh)

        # 定义使用统一界面的方法（不需要传统的残基选择界面）
        unified_methods = ['DensityRadius', 'DensityMultiRadius', 'Area', 'Gyration', 'Anisotropy']
        
        if self.method not in unified_methods:
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
        elif self.method == 'DensityRadius':
            print('DensityRadius')
            # 为DensityRadius直接创建SpinBox界面
            widget = UIItemsMake.make_widget()
            layout = QVBoxLayout(widget)
            self.ui.widgetspinbox = widget
            self.ui.VLayoutspinbox = layout
            
            # 创建MW SpinBox
            label_mw = UIItemsMake.make_label('MW(g/mol)')
            spin_box_mw = UIItemsMake.make_spin_box(14, 0, 1000000)
            self.ui.widgetspinboxLabel0 = label_mw
            self.ui.widgetspinboxSpinBox0 = spin_box_mw
            layout.addWidget(label_mw)
            layout.addWidget(spin_box_mw)
            
            # 创建Radius SpinBox
            label_radius = UIItemsMake.make_label('Radius(A)')
            spin_box_radius = UIItemsMake.make_spin_box(50, 0, 1000000)
            self.ui.widgetspinboxLabel1 = label_radius
            self.ui.widgetspinboxSpinBox1 = spin_box_radius
            layout.addWidget(label_radius)
            layout.addWidget(spin_box_radius)
            
            # 创建Next按钮
            btnNext = UIItemsMake.make_btn('Next')
            layout.addWidget(btnNext)
            self.ui.stackedWidget_Analysis.addWidget(widget)
            self.ui.stackedWidget_Analysis.setCurrentIndex(0)
            
            # 创建DensityRadiusHandler实例并连接按钮事件
            self.density_handler = DensityRadiusHandler(self.ui, self.Box, self.Info)
            btnNext.clicked.connect(lambda: self.density_handler.step_1())
            print("DensityRadius button connected successfully")
        elif self.method == 'DensityMultiRadius':
            print('DensityMultiRadius')
            # 为DensityMultiRadius直接创建SpinBox界面
            widget = UIItemsMake.make_widget()
            layout = QVBoxLayout(widget)
            self.ui.widgetspinbox = widget
            self.ui.VLayoutspinbox = layout
            
            # 创建MW SpinBox
            label_mw = UIItemsMake.make_label('MW(g/mol)')
            spin_box_mw = UIItemsMake.make_spin_box(14, 0, 1000000)
            self.ui.widgetspinboxLabel0 = label_mw
            self.ui.widgetspinboxSpinBox0 = spin_box_mw
            layout.addWidget(label_mw)
            layout.addWidget(spin_box_mw)
            
            # 创建MaxRadius SpinBox
            label_max_radius = UIItemsMake.make_label('MaxRadius(A)')
            spin_box_max_radius = UIItemsMake.make_spin_box(50, 0, 1000000)
            self.ui.widgetspinboxLabel1 = label_max_radius
            self.ui.widgetspinboxSpinBox1 = spin_box_max_radius
            layout.addWidget(label_max_radius)
            layout.addWidget(spin_box_max_radius)
            
            # 创建NumberSegments SpinBox
            label_number_segments = UIItemsMake.make_label('NumberSegments')
            spin_box_number_segments = UIItemsMake.make_spin_box(5, 1, 1000)
            self.ui.widgetspinboxLabel2 = label_number_segments
            self.ui.widgetspinboxSpinBox2 = spin_box_number_segments
            layout.addWidget(label_number_segments)
            layout.addWidget(spin_box_number_segments)
            
            # 创建Next按钮
            btnNext = UIItemsMake.make_btn('Next')
            layout.addWidget(btnNext)
            self.ui.stackedWidget_Analysis.addWidget(widget)
            self.ui.stackedWidget_Analysis.setCurrentIndex(0)
            
            # 创建DensityMultiRadiusHandler实例并连接按钮事件
            self.density_multi_handler = DensityMultiRadiusHandler(self.ui, self.Box, self.Info)
            btnNext.clicked.connect(lambda: self.density_multi_handler.step_1())
            print("DensityMultiRadius button connected successfully")
        elif self.method in ['Area', 'Gyration', 'Anisotropy']:
            print(f'{self.method} - Using unified handler')
            # 对于使用UnifiedAnalysisHandler的方法，创建对应的Handler实例
            handler_strategies = {
                'Area': AreaHandler,
                'Gyration': GyrationHandler,
                'Anisotropy': AnisotropyHandler
            }
            
            handler_class = handler_strategies[self.method]
            handler_instance = handler_class(self.ui, self.Box, self.Info)
            
            # 保存Handler实例到self中，防止被垃圾回收
            print(f"*** Saving Handler instance for {self.method} ***")
            if self.method == 'Area':
                self.area_handler = handler_instance
                print(f"*** Saved area_handler: {self.area_handler} ***")
            elif self.method == 'Gyration':
                self.gyration_handler = handler_instance
                print(f"*** Saved gyration_handler: {self.gyration_handler} ***")
            elif self.method == 'Anisotropy':
                self.anisotropy_handler = handler_instance
                print(f"*** Saved anisotropy_handler: {self.anisotropy_handler} ***")
            
            # 调用start方法开始分析流程
            handler_instance.start()
            
            # 确保按钮被正确设置为全局运行按钮
            if hasattr(self.ui, 'btnAnalysisRun'):
                print(f"Global run button set for {self.method}")
            else:
                print(f"Warning: Global run button not set for {self.method}")
            
            # 为Area、Gyration、Anisotropy添加back/refresh按钮和stackedWidget_Analysis
            print(f"*** Adding buttons to main layout for {self.method} ***")
            self.ui.VLayoutRightMain.addWidget(self.ui.widgetUpNextBack)
            self.ui.VLayoutRightMain.addWidget(self.ui.stackedWidget_Analysis)

        # 如果方法不是DensityRadius、DensityMultiRadius、Area、Gyration、Anisotropy，则使用layout_strategies处理
        if self.method not in ['DensityRadius', 'DensityMultiRadius', 'Area', 'Gyration', 'Anisotropy']:
            layout_strategies = {
                'Height': HeightLayout
                , 'SZ': SZLayout
                # , 'MeanCurvature': MeanCurvatureLayout
                # , 'RadialDistribution': RDLayout
                , 'Cluster': ClusterLayout
                # , 'NCluster': NClusterLayout
                # , 'PCA': PCALayout
            }

            layout_analysis = layout_strategies[self.method]
            self._create_layout(layout_analysis)

            self.ui.VLayoutRightMain.addWidget(self.ui.widgetUpNextBack)
            self.ui.VLayoutRightMain.addWidget(self.ui.stackedWidget_Analysis)
        elif self.method in ['DensityRadius', 'DensityMultiRadius']:
            # DensityRadius和DensityMultiRadius方法需要添加back/refresh按钮和stackedWidget_Analysis
            print(f"*** Adding buttons to main layout for {self.method} ***")
            self.ui.VLayoutRightMain.addWidget(self.ui.widgetUpNextBack)
            self.ui.VLayoutRightMain.addWidget(self.ui.stackedWidget_Analysis)

    def _create_layout(self, layout):
        init_layout = layout(self.ui, self.Box, self.Info)
        
        # 定义使用统一界面的方法（不需要传统的残基选择按钮）
        unified_methods = ['Area', 'Gyration', 'Anisotropy']
        
        if self.method not in unified_methods:
            # 传统方法：创建残基选择按钮
            btn = UIItemsMake.make_btn('Next'
                                       , callback=lambda: AnalysisUtils.get_residue(self.Box.residues_list,
                                                                                    self.Box.residue_click)
                                       , layout=self.ui.VLayoutResidue)
            btn.clicked.connect(lambda: init_layout.step_1())
        # 对于使用AtomsLayoutUnified的方法，不需要额外的按钮，因为按钮已经在AtomsLayoutUnified中创建了

    def switchWidgetBack(self):
        if self.method == 'DensityRadius' and self.density_handler:
            # DensityRadius的back功能
            self.density_handler.go_back()
        elif self.method == 'DensityMultiRadius' and self.density_multi_handler:
            # DensityMultiRadius的back功能
            self.density_multi_handler.go_back()
        elif self.method in ['Area', 'Gyration', 'Anisotropy']:
            # Area、Gyration、Anisotropy的back功能 - 由于使用单步界面，back应该返回到方法选择
            print(f"Back button clicked for {self.method} - refreshing to method selection")
            # 调用对应Handler的refresh方法
            if self.method == 'Area' and self.area_handler:
                self.area_handler.refresh()
            elif self.method == 'Gyration' and self.gyration_handler:
                self.gyration_handler.refresh()
            elif self.method == 'Anisotropy' and self.anisotropy_handler:
                self.anisotropy_handler.refresh()
            else:
                self.refreshWidget()
        else:
            # 传统分析方法的back功能
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

    def update_button_icons(self):
        """
        根据当前主题更新按钮图标
        """
        try:
            # 获取主题管理器
            theme_manager = None
            
            # 尝试从全局变量获取主题管理器
            import main
            if hasattr(main, 'window') and main.window and hasattr(main.window, 'theme_manager'):
                theme_manager = main.window.theme_manager
                print("通过全局变量获取到主题管理器")
            else:
                # 尝试从UI获取主题管理器
                if hasattr(self.ui, 'theme_manager') and self.ui.theme_manager:
                    theme_manager = self.ui.theme_manager
                    print("通过UI获取到主题管理器")
                else:
                    # 默认使用白天主题
                    is_dark = False
                    print("未找到主题管理器，使用默认白天主题")
            
            # 确定当前是否为暗色主题
            if theme_manager:
                is_dark = theme_manager.is_dark_theme()
                print(f"当前主题: {'暗色' if is_dark else '亮色'}")
            else:
                is_dark = False
            
            # 根据主题设置图标路径
            if is_dark:
                back_icon_path = 'images/icons/houtui_black.png'
                refresh_icon_path = 'images/icons/shuaxin_black.png'
            else:
                back_icon_path = 'images/icons/houtui.png'
                refresh_icon_path = 'images/icons/shuaxin.png'
            
            # 更新按钮图标
            if hasattr(self.ui, 'btnBack'):
                self.ui.btnBack.setIcon(QIcon(back_icon_path))
                print(f"Back按钮图标已更新: {back_icon_path}")
            
            if hasattr(self.ui, 'btnRefresh'):
                self.ui.btnRefresh.setIcon(QIcon(refresh_icon_path))
                print(f"Refresh按钮图标已更新: {refresh_icon_path}")
                
        except Exception as e:
            print(f"更新按钮图标时出错: {e}")
            # 如果出错，使用默认的白天主题图标
            if hasattr(self.ui, 'btnBack'):
                self.ui.btnBack.setIcon(QIcon('images/icons/houtui.png'))
            if hasattr(self.ui, 'btnRefresh'):
                self.ui.btnRefresh.setIcon(QIcon('images/icons/shuaxin.png'))

    def refreshWidget(self):
        print(f"*** refreshWidget() called for method: {self.method} ***")
        print(f"*** Current time: {time.time()} ***")
        
        # 检查Handler实例是否存在
        print(f"*** Checking Handler instances ***")
        print(f"hasattr(self, 'area_handler'): {hasattr(self, 'area_handler')}")
        print(f"hasattr(self, 'gyration_handler'): {hasattr(self, 'gyration_handler')}")
        print(f"hasattr(self, 'anisotropy_handler'): {hasattr(self, 'anisotropy_handler')}")
        
        if hasattr(self, 'area_handler'):
            print(f"self.area_handler: {self.area_handler}")
        if hasattr(self, 'gyration_handler'):
            print(f"self.gyration_handler: {self.gyration_handler}")
        if hasattr(self, 'anisotropy_handler'):
            print(f"self.anisotropy_handler: {self.anisotropy_handler}")
        
        # 所有分析方法都使用相同的refresh功能：清空布局并重新读取文件
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
        
        # 对于Area、Gyration、Anisotropy，使用传统的refresh方式重新读取文件
        if self.method in ['Area', 'Gyration', 'Anisotropy']:
            print(f"*** Using traditional refresh for {self.method} to re-read file ***")
            clearLayout(self.ui.VLayoutRightMain)
            self.readFile()
        else:
            # 其他分析方法使用传统的refresh方式
            print(f"*** Using traditional refresh for {self.method} ***")
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
            if not hasattr(self, 'figureTypeWidget'):
                self.ui.btnAnalysisRun.deleteLater()
                self._create_figure_type_selection()
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
            
            # 删除进度条
            if hasattr(self.ui, 'progressBar'):
                self.ui.progressBar.deleteLater()

    def _create_figure_type_selection(self):
        """创建图表类型选择界面"""
        # 创建图表类型选择widget
        self.ui.figureTypeWidget = QWidget()
        layout = QVBoxLayout(self.ui.figureTypeWidget)
        
        # 添加标题
        title_label = UIItemsMake.make_label('Select Figure Type')
        title_label.setStyleSheet(f"font-size: {AnalysisBtnClick.FONT_SIZE}; font-weight: bold;")
        layout.addWidget(title_label)
        
        # 获取分析类支持的图表类型
        if hasattr(self, 'cls') and self.cls and hasattr(self.cls, 'supported_figure_types'):
            figure_types = self.cls.supported_figure_types
            print(f"Available figure types for {type(self.cls).__name__}: {figure_types}")
        else:
            # 默认图表类型
            figure_types = ['Line Chart', 'Bar Chart']
            print(f"No supported_figure_types found, using default: {figure_types}")
        
        # 创建图表类型按钮
        for figure_type in figure_types:
            btn = UIItemsMake.make_btn(f'Create {figure_type}')
            btn.clicked.connect(lambda checked, ft=figure_type: self._create_figure(ft))
            layout.addWidget(btn)
        
        # 添加到主布局
        self.ui.VLayoutRightMain.addWidget(self.ui.figureTypeWidget)

    def _create_figure(self, figure_type):
        """创建指定类型的图表"""
        print(f"Creating {figure_type} for {self.Info.method_analysis} analysis")
        print(f"Result path: {self.Info.path_result}")
        print(f"Analysis completed in: {time.strftime('%M:%S', time.gmtime(time.time() - self.start_time))}")
        
        # 调用分析类的绘图方法
        if hasattr(self, 'cls') and self.cls:
            print(f"Analysis class: {type(self.cls).__name__}")
            
            # 根据图表类型调用相应的绘图方法
            if figure_type == 'Line Chart':
                if hasattr(self.cls, 'plot_line'):
                    print("Calling plot_line method...")
                    self.cls.plot_line()
                else:
                    print(f"{type(self.cls).__name__} does not support plot_line method")
            elif figure_type == 'Bar Chart':
                if hasattr(self.cls, 'plot_bar'):
                    print("Calling plot_bar method...")
                    self.cls.plot_bar()
                else:
                    print(f"{type(self.cls).__name__} does not support plot_bar method")
            elif figure_type == 'Scatter Plot':
                if hasattr(self.cls, 'plot_scatter'):
                    print("Calling plot_scatter method...")
                    self.cls.plot_scatter()
                else:
                    print(f"{type(self.cls).__name__} does not support plot_scatter method")
            elif figure_type == 'Heatmap':
                if hasattr(self.cls, 'plot_heatmap'):
                    print("Calling plot_heatmap method...")
                    self.cls.plot_heatmap()
                else:
                    print(f"{type(self.cls).__name__} does not support plot_heatmap method")
            elif figure_type == '3D Surface':
                if hasattr(self.cls, 'plot_3d_surface'):
                    print("Calling plot_3d_surface method...")
                    self.cls.plot_3d_surface()
                else:
                    print(f"{type(self.cls).__name__} does not support plot_3d_surface method")
            else:
                print(f"Unknown figure type: {figure_type}")
        else:
            print("No analysis class available for plotting")

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
                    , clickresidues=True
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

        for sp in self.Box.residue_click if clickresidues else self.Box.residues:
            groupBox = UIItemsMake.make_group_box(sp)
            groupBoxLayout = QVBoxLayout(groupBox)
            for atom in self.Box.residues_atoms[sp] if clickresidues else self.Box.residues_atoms[sp]:
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

class UnifiedAnalysisHandler:
    """统一分析处理器基类，用于Area、Gyration、Anisotropy等使用AtomsLayoutUnified的方法"""
    
    def __init__(self, ui, Box, Info, analysis_class, config_name, atom_type='checkbox'):
        self.ui = ui
        self.Box = Box
        self.Info = Info
        self.analysis_class = analysis_class
        self.config_name = config_name
        self.atom_type = atom_type
        self.current_step = 0
        self.step_widgets = []
        self.start_time = None
    
    def start(self):
        """开始分析流程"""
        print(f"{self.config_name} - Starting unified analysis")
        self.current_step = 0
        self._create_unified_selection_layout()
    
    def _create_unified_selection_layout(self):
        """创建统一的选择界面"""
        print(f"*** _create_unified_selection_layout() called for {self.config_name} ***")
        print(f"*** Current time: {time.time()} ***")
        
        widget = UIItemsMake.make_widget()
        layout = QVBoxLayout(widget)
        
        # 设置widget属性
        setattr(self.ui, f'widgetAtoms{self.config_name}', widget)
        setattr(self.ui, f'VLayoutAtoms{self.config_name}', layout)
        
        # 添加标题
        title_label = UIItemsMake.make_label(f'Select Residues and Atoms for {self.config_name} Analysis')
        title_label.setStyleSheet(f"font-size: {AnalysisBtnClick.FONT_SIZE}; font-weight: bold;")
        layout.addWidget(title_label)
        
        # 创建滚动区域
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # 清空groups列表
        groups = getattr(self.Box.get_config(self.config_name), f'{self.config_name}HeadGroups')
        groups.clear()
        
        # 为每个残基创建选择界面
        for residue_name in self.Box.residues:
            # 创建残基组
            residue_group = QGroupBox(f"{residue_name}")
            residue_group.setCheckable(True)
            residue_group.setChecked(False)
            residue_layout = QVBoxLayout(residue_group)
            
            # 获取该残基的所有原子
            atoms = self.Box.residues_atoms[residue_name]
            atom_checkboxes = []
            
            for atom_name in atoms:
                if self.atom_type == 'radio':
                    atom_widget = QRadioButton(atom_name)
                else:  # 默认使用checkbox
                    atom_widget = QCheckBox(atom_name)
                atom_widget.setEnabled(False)  # 默认禁用，只有选中残基时才启用
                atom_layout = QHBoxLayout()
                atom_layout.addWidget(atom_widget)
                atom_layout.addStretch()
                residue_layout.addLayout(atom_layout)
                atom_checkboxes.append(atom_widget)
            
            # 当残基被选中时，启用原子选择
            def make_residue_handler(residue_group, atom_checkboxes):
                def handler():
                    enabled = residue_group.isChecked()
                    for checkbox in atom_checkboxes:
                        checkbox.setEnabled(enabled)
                return handler
            
            residue_group.toggled.connect(make_residue_handler(residue_group, atom_checkboxes))
            
            # 存储残基组信息
            residue_group.atom_checkboxes = atom_checkboxes
            residue_group.residue_name = residue_name
            groups.append(residue_group)
            
            scroll_layout.addWidget(residue_group)
        
        scroll_area.setWidget(scroll_widget)
        scroll_area.setWidgetResizable(True)
        layout.addWidget(scroll_area)
        
        # 添加Run按钮
        btnRun = UIItemsMake.make_btn('Run!')
        layout.addWidget(btnRun)
        
        # 添加到堆叠窗口
        self.ui.stackedWidget_Analysis.addWidget(widget)
        self.ui.stackedWidget_Analysis.setCurrentIndex(0)
        
        # 连接按钮事件
        print(f"Connecting button for {self.config_name}")
        print(f"Button object: {btnRun}")
        print(f"Button text: {btnRun.text()}")
        
        # 添加一个简单的测试连接
        def test_click():
            print(f"TEST: Button clicked for {self.config_name}")
        
        btnRun.clicked.connect(test_click)
        btnRun.clicked.connect(self.run)
        
        setattr(self.ui, 'btnAnalysisRun', btnRun)
        print(f"Button connected to run and set as global run button for {self.config_name}")
        print(f"Global button set: {hasattr(self.ui, 'btnAnalysisRun')}")
        print(f"Global button object: {getattr(self.ui, 'btnAnalysisRun', 'Not found')}")
    
    def _collect_selections(self):
        """收集选择信息"""
        print(f"UnifiedAnalysisHandler._collect_selections() called for {self.config_name}")
        
        selected_atoms = {}
        # 使用正确的属性名称 - 统一使用全称+HeadGroups
        groups_attr_name = f'{self.config_name}HeadGroups'
        
        print(f"Getting groups from config: {groups_attr_name}")
        groups = getattr(self.Box.get_config(self.config_name), groups_attr_name)
        print(f"Groups: {groups}")
        print(f"Groups type: {type(groups)}")
        print(f"Groups length: {len(groups) if hasattr(groups, '__len__') else 'No length'}")
        
        print("Iterating through groups:")
        for i, residue_group in enumerate(groups):
            print(f"  Group {i}: {residue_group}")
            print(f"    Checked: {residue_group.isChecked()}")
            print(f"    Residue name: {getattr(residue_group, 'residue_name', 'No residue_name')}")
            print(f"    Atom checkboxes: {getattr(residue_group, 'atom_checkboxes', 'No atom_checkboxes')}")
            
            if residue_group.isChecked():
                residue_name = residue_group.residue_name
                selected_atoms[residue_name] = []
                print(f"    Processing checked residue: {residue_name}")
                
                for j, checkbox in enumerate(residue_group.atom_checkboxes):
                    print(f"      Checkbox {j}: {checkbox.text()}, checked: {checkbox.isChecked()}")
                    if checkbox.isChecked():
                        selected_atoms[residue_name].append(checkbox.text())
                        print(f"        Added atom: {checkbox.text()}")
        
        print(f"Final selected_atoms: {selected_atoms}")
        
        # 正确设置到配置中
        print(f"Setting {self.config_name}HeadAtoms in config")
        setattr(self.Box.get_config(self.config_name), f'{self.config_name}HeadAtoms', selected_atoms)
        
        # 验证设置是否成功
        verify_atoms = getattr(self.Box.get_config(self.config_name), f'{self.config_name}HeadAtoms')
        print(f"Verification - {self.config_name}HeadAtoms: {verify_atoms}")
        
        print(f"Collected selections for {self.config_name}: {len(selected_atoms)} residue groups")
    
    def run(self):
        """执行分析"""
        print(f"*** {self.config_name}Layout.run() called ***")
        print(f"*** Method: {self.config_name} ***")
        print(f"*** Time: {time.time()} ***")
        self.start_time = time.time()
        
        # 收集选择
        print(f"Step 1: Collecting selections for {self.config_name}")
        self._collect_selections()
        
        # 创建分析类
        print(f"Step 2: Creating analysis class for {self.config_name}")
        self._create_analysis_class()
        
        # 启动分析
        print(f"Step 3: Starting analysis for {self.config_name}")
        self._start_analysis()
    
    def _create_analysis_class(self):
        """创建分析类 - 子类需要重写此方法"""
        raise NotImplementedError("Subclasses must implement _create_analysis_class")
    
    def _start_analysis(self):
        """启动分析"""
        print(f"Checking if analysis class exists for {self.config_name}")
        print(f"hasattr(self, 'cls'): {hasattr(self, 'cls')}")
        if hasattr(self, 'cls'):
            print(f"self.cls is not None: {self.cls is not None}")
            print(f"self.cls type: {type(self.cls)}")
        
        if hasattr(self, 'cls') and self.cls is not None:
            print(f"Creating progress bar for {self.config_name}")
            # 添加进度条
            if not hasattr(self.ui, 'progressBar'):
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
                print(f"Progress bar created for {self.config_name}")
            else:
                print(f"Progress bar already exists for {self.config_name}")
            
            print(f"Creating Worker for {self.config_name}")
            print(f"Worker parameters: frames={self.Info.frame_first}-{self.Info.frame_last}, step={self.Info.step}")
            
            # 创建Worker并启动
            worker = Worker(self.cls, self.Info.frame_first, self.Info.frame_last, self.Info.step)
            worker.progressValueChanged.connect(self.updateProgressBar)
            self.thread = Thread(target=worker.run)
            self.thread.start()
            print(f"{self.config_name} analysis started successfully")
        else:
            print(f"ERROR: {self.config_name} analysis class not created or is None")
            print(f"hasattr(self, 'cls'): {hasattr(self, 'cls')}")
            if hasattr(self, 'cls'):
                print(f"self.cls value: {self.cls}")
    
    def updateProgressBar(self, value):
        """更新进度条"""
        # 检查进度条是否还存在，防止重复删除导致的RuntimeError
        if not hasattr(self.ui, 'progressBar') or self.ui.progressBar is None:
            return  # 如果进度条已被删除，直接返回
        
        self.ui.progressBar.setValue(value)
        if value == 100:
            if not hasattr(self.ui, 'figureTypeWidget'):
                if hasattr(self.ui, 'btnAnalysisRun'):
                    self.ui.btnAnalysisRun.deleteLater()
                self._create_figure_type_selection()
            
            # 计算分析时间并显示成功消息
            end_time = time.time()
            elapsed_time = end_time - self.start_time
            formatted_time = time.strftime('%M:%S', time.gmtime(elapsed_time))
            success_message = (
                'Analysis Completed\n'
                f"Time taken: {formatted_time}\n"
                "The gro file and topol file were saved at:\n"
                f"{self.Info.path_result}"
            )
            create_warn_dialog(success_message, 'Analysis')
            
            # 删除进度条并设置为None，防止重复访问
            if hasattr(self.ui, 'progressBar') and self.ui.progressBar is not None:
                self.ui.progressBar.deleteLater()
                self.ui.progressBar = None  # 设置为None，防止重复访问
    
    def _create_figure_type_selection(self):
        """创建图表类型选择界面"""
        self.ui.figureTypeWidget = QWidget()
        layout = QVBoxLayout(self.ui.figureTypeWidget)
        
        # 添加标题
        title_label = UIItemsMake.make_label('Select Figure Type')
        title_label.setStyleSheet(f"font-size: {AnalysisBtnClick.FONT_SIZE}; font-weight: bold;")
        layout.addWidget(title_label)
        
        # 获取分析类支持的图表类型
        if hasattr(self, 'cls') and self.cls and hasattr(self.cls, 'supported_figure_types'):
            figure_types = self.cls.supported_figure_types
            print(f"Available figure types for {type(self.cls).__name__}: {figure_types}")
        else:
            # 默认图表类型
            figure_types = ['Line Chart', 'Bar Chart']
            print(f"No supported_figure_types found, using default: {figure_types}")
        
        # 创建图表类型按钮
        for figure_type in figure_types:
            btn = UIItemsMake.make_btn(f'Create {figure_type}')
            btn.clicked.connect(lambda checked, ft=figure_type: self._create_figure(ft))
            layout.addWidget(btn)
        
        # 添加到主布局
        self.ui.VLayoutRightMain.addWidget(self.ui.figureTypeWidget)
    
    def _create_figure(self, figure_type):
        """创建指定类型的图表"""
        print(f"Creating {figure_type} for {self.config_name} analysis")
        print(f"Result path: {self.Info.path_result}")
        print(f"Analysis completed in: {time.strftime('%M:%S', time.gmtime(time.time() - self.start_time))}")
        
        # 调用分析类的绘图方法
        if hasattr(self, 'cls') and self.cls:
            print(f"Analysis class: {type(self.cls).__name__}")
            
            # 根据图表类型调用相应的绘图方法
            if figure_type == 'Line Chart':
                if hasattr(self.cls, 'plot_line'):
                    print("Calling plot_line method...")
                    self.cls.plot_line()
                else:
                    print(f"{type(self.cls).__name__} does not support plot_line method")
            elif figure_type == 'Bar Chart':
                if hasattr(self.cls, 'plot_bar'):
                    print("Calling plot_bar method...")
                    self.cls.plot_bar()
                else:
                    print(f"{type(self.cls).__name__} does not support plot_bar method")
            elif figure_type == 'Scatter Plot':
                if hasattr(self.cls, 'plot_scatter'):
                    print("Calling plot_scatter method...")
                    self.cls.plot_scatter()
                else:
                    print(f"{type(self.cls).__name__} does not support plot_scatter method")
            else:
                print(f"Unknown figure type: {figure_type}")
        else:
            print("No analysis class available for plotting")

    

    def AtomsLayoutUnified(self, widgetName: str, widgetLayout: str, labelText: str,
                          groups, btnName: str, stackID, func, connect=False, atom_type='checkbox'):
        """
        创建统一的残基和原子选择界面（类似Density的界面）
        :param widgetName: 创建新的窗口名称
        :param widgetLayout: 创建新的布局的名称
        :param labelText: 标签文本
        :param groups: 储存信息的Group列表
        :param btnName: 按钮名称
        :param stackID: 堆叠窗口ID
        :param func: 按钮点击函数
        :param connect: 是否将按钮设置为全局运行按钮
        """
        widget = UIItemsMake.make_widget()
        layout = QVBoxLayout(widget)
        
        # 设置widget属性
        setattr(self.ui, widgetName, widget)
        setattr(self.ui, widgetLayout, layout)
        
        # 添加标题
        title_label = UIItemsMake.make_label(labelText)
        title_label.setStyleSheet(f"font-size: {AnalysisBtnClick.FONT_SIZE}; font-weight: bold;")
        layout.addWidget(title_label)
        
        # 创建滚动区域
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # 清空groups列表
        groups.clear()
        
        # 为每个残基创建选择界面
        for residue_name in self.Box.residues:
            # 创建残基组
            residue_group = QGroupBox(f"{residue_name}")
            residue_group.setCheckable(True)
            residue_group.setChecked(False)
            residue_layout = QVBoxLayout(residue_group)
            
            # 获取该残基的所有原子
            atoms = self.Box.residues_atoms[residue_name]
            atom_checkboxes = []
            
            for atom_name in atoms:
                if atom_type == 'radio':
                    atom_widget = QRadioButton(atom_name)
                else:  # 默认使用checkbox
                    atom_widget = QCheckBox(atom_name)
                atom_widget.setEnabled(False)  # 默认禁用，只有选中残基时才启用
                atom_layout = QHBoxLayout()
                atom_layout.addWidget(atom_widget)
                atom_layout.addStretch()
                residue_layout.addLayout(atom_layout)
                atom_checkboxes.append(atom_widget)
            
            # 当残基被选中时，启用原子选择
            def make_residue_handler(residue_group, atom_checkboxes):
                def handler():
                    enabled = residue_group.isChecked()
                    for checkbox in atom_checkboxes:
                        checkbox.setEnabled(enabled)
                return handler
            
            residue_group.toggled.connect(make_residue_handler(residue_group, atom_checkboxes))
            
            # 存储残基组信息
            residue_group.atom_checkboxes = atom_checkboxes
            residue_group.residue_name = residue_name
            groups.append(residue_group)
            
            scroll_layout.addWidget(residue_group)
        
        scroll_area.setWidget(scroll_widget)
        scroll_area.setWidgetResizable(True)
        layout.addWidget(scroll_area)
        
        # 添加Next按钮
        btnNext = UIItemsMake.make_btn(btnName)
        layout.addWidget(btnNext)
        
        # 添加到堆叠窗口
        self.ui.stackedWidget_Analysis.addWidget(widget)
        self.ui.stackedWidget_Analysis.setCurrentIndex(stackID)
        
        # 连接按钮事件
        btnNext.clicked.connect(func)
        if connect: 
            setattr(self.ui, 'btnAnalysisRun', btnNext)
            print(f"Button connected to {func.__name__} and set as global run button")
        else:
            print(f"Button connected to {func.__name__} but not set as global run button")

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
                          , filePath=self.Info.path_result
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
                          , filePath=self.Info.path_result
                          , k=self.Info.K)
        else:
            print('暂时不支持%s力场' % self.FF_TYPE)




class AreaHandler(UnifiedAnalysisHandler):
    """Area分析处理器"""
    
    def __init__(self, ui, Box, Info):
        super().__init__(ui, Box, Info, Area, 'Area', atom_type='radio')
    
    def _create_analysis_class(self):
        """创建Area分析类"""
        print("AreaHandler._create_analysis_class() called")
        
        print("Calling self.Info.get_text()")
        self.Info.get_text()
        
        print("Getting AreaHeadAtoms from config")
        area_head_atoms = self.Box.get_config('Area').AreaHeadAtoms
        print(f"AreaHeadAtoms: {area_head_atoms}")
        print(f"AreaHeadAtoms type: {type(area_head_atoms)}")
        print(f"AreaHeadAtoms length: {len(area_head_atoms) if isinstance(area_head_atoms, dict) else 'Not a dict'}")
        
        print("Getting other parameters")
        print(f"universe: {self.Box.u}")
        print(f"filePath: {self.Info.path_result}")
        print(f"k: {self.Info.K}")
        
        print("Creating Area analysis class")
        try:
            self.cls = Area(universe=self.Box.u
                            , residueGroup=area_head_atoms
                            , filePath=self.Info.path_result
                            , k=self.Info.K)
            print("Area analysis class created successfully")
            print(f"self.cls: {self.cls}")
            print(f"self.cls type: {type(self.cls)}")
        except Exception as e:
            print(f"ERROR creating Area analysis class: {e}")
            import traceback
            traceback.print_exc()
            self.cls = None


class GyrationHandler(UnifiedAnalysisHandler):
    """Gyration分析处理器"""
    
    def __init__(self, ui, Box, Info):
        super().__init__(ui, Box, Info, Gyration, 'Gyration', atom_type='radio')
    
    def _create_analysis_class(self):
        """创建Gyration分析类"""
        self.Info.get_text()
        print(f"Gyration analysis starting with {len(self.Box.get_config('Gyration').GyrationHeadAtoms)} residue groups")
        self.cls = Gyration(self.Box.u
                           , self.Box.get_config('Gyration').GyrationHeadAtoms
                           , filePath=self.Info.path_result)
        print("Gyration analysis class created successfully")


class AnisotropyHandler(UnifiedAnalysisHandler):
    """Anisotropy分析处理器"""
    
    def __init__(self, ui, Box, Info):
        super().__init__(ui, Box, Info, Anisotropy, 'Anisotropy', atom_type='radio')
    
    def _create_analysis_class(self):
        """创建Anisotropy分析类"""
        self.Info.get_text()
        print(f"Anisotropy analysis starting with {len(self.Box.get_config('Anisotropy').AnisotropyHeadAtoms)} residue groups")
        self.cls = Anisotropy(self.Box.u
                             , self.Box.get_config('Anisotropy').AnisotropyHeadAtoms
                             , filePath=self.Info.path_result)
        print("Anisotropy analysis class created successfully")






class ClusterLayout(AnalysisLayout):
    def step_1(self):
        self.AtomsLayout('widgetAtomsCluster'
                         , 'VLayoutAtomsCluster'
                         , 'select atoms'
                         , self.Box.get_config('Cluster').ClusterHeadGroups
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
        AnalysisUtils.get_atom(self.Box.get_config('Cluster').ClusterHeadGroups
                               , QCheckBox
                               , self.Box.get_config('Cluster').ClusterHeadAtoms)
        self.Box.get_config('Cluster').ClusterCutoff = AnalysisUtils.get_spin_value(self.ui.widgetCutoffSpinBox)
        self.cls = Cluster(self.Box.u
                           , self.Box.get_config('Cluster').ClusterHeadAtoms
                           , filePath=self.Info.path_result
                           , cutoff=self.Box.get_config('Cluster').ClusterCutoff)

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
    __slots__ = ('MeanCurvatureHeadAtoms', 'MeanCurvatureHeadGroups')

    def __init__(self):
        self.MeanCurvatureHeadAtoms = {}
        self.MeanCurvatureHeadGroups = []


class ASPConfig(BaseConfig):
    __slots__ = ('AnisotropyHeadAtoms', 'AnisotropyHeadGroups')

    def __init__(self):
        self.AnisotropyHeadAtoms = {}
        self.AnisotropyHeadGroups = []


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
    __slots__ = ('GyrationHeadAtoms', 'GyrationHeadGroups')

    def __init__(self):
        self.GyrationHeadAtoms = {}
        self.GyrationHeadGroups = []


class CLConfig(BaseConfig):
    __slots__ = ('ClusterHeadAtoms', 'ClusterHeadGroups', 'ClusterCutoff')

    def __init__(self):
        self.ClusterHeadAtoms = {}
        self.ClusterHeadGroups = []
        self.ClusterCutoff = 12

class DensityRadiusConfig(BaseConfig):
    __slots__ = ('DensityHeadAtoms', 'DensityGroups', 'DensityMW', 'GasGroups', 'GasHeadAtoms', 'DensityRadius', 'DensityNumberSegments')

    def __init__(self):
        self.DensityHeadAtoms = {}
        self.DensityGroups = []
        self.DensityMW = 14
        self.GasGroups = []
        self.GasHeadAtoms = {}
        self.DensityRadius = 50
        self.DensityNumberSegments = 5


class DensityMultiRadiusConfig(BaseConfig):
    __slots__ = ('DensityHeadAtoms', 'DensityGroups', 'DensityMW', 'GasGroups', 'GasHeadAtoms', 'DensityMaxRadius', 'DensityNumberSegments')

    def __init__(self):
        self.DensityHeadAtoms = {}
        self.DensityGroups = []
        self.DensityMW = 14
        self.GasGroups = []
        self.GasHeadAtoms = {}
        self.DensityMaxRadius = 50
        self.DensityNumberSegments = 5

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
        , 'DensityRadius': DensityRadiusConfig
        , 'DensityMultiRadius': DensityMultiRadiusConfig
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


class DensityRadiusHandler:
    """专门处理DensityRadius分析流程的类，不继承AtomsLayout"""
    
    def __init__(self, ui, Box, Info):
        self.ui = ui
        self.Box = Box
        self.Info = Info
        self.start_time = None
        self.current_step = 0
        self.step_widgets = []  # 存储每个步骤的widget
    
    def start(self):
        """开始DensityRadius分析流程"""
        print("DensityRadiusHandler.start() called")
        self.step_1()
    
    def step_1_from_mw_radius(self):
        """从MW和Radius界面进入下一步"""
        print("DensityRadiusHandler.step_1_from_mw_radius() called")
        # 保存MW和Radius参数
        self.Box.get_config('DensityRadius').DensityMW = AnalysisUtils.get_spin_value(self.ui.widgetspinboxSpinBox0)
        self.Box.get_config('DensityRadius').DensityRadius = AnalysisUtils.get_spin_value(self.ui.widgetspinboxSpinBox1)
        
        print(f"MW: {self.Box.get_config('DensityRadius').DensityMW}, Radius: {self.Box.get_config('DensityRadius').DensityRadius}")
        
        # 更新步骤
        self.current_step = 1
        
        # 创建残基和原子选择界面
        self._create_residue_atom_selection_layout()
    
    def step_1(self):
        """保存MW和Radius参数并进入下一步（兼容性方法）"""
        self.step_1_from_mw_radius()
    
    def _create_residue_atom_selection_layout(self):
        """创建残基和原子选择界面"""
        widget = UIItemsMake.make_widget()
        layout = QVBoxLayout(widget)
        
        # 设置widget属性
        self.ui.widgetResidueAtomDensity = widget
        self.ui.VLayoutResidueAtomDensity = layout
        
        # 添加标题
        title_label = UIItemsMake.make_label('Select Residues Group')
        title_label.setStyleSheet(f"font-size: {AnalysisBtnClick.FONT_SIZE}; font-weight: bold;")
        layout.addWidget(title_label)
        
        # 创建滚动区域
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # 为每个残基创建选择界面
        self.Box.get_config('DensityRadius').DensityGroups = []
        for residue_name in self.Box.residues:
            # 创建残基组
            residue_group = QGroupBox(f"{residue_name}")
            residue_group.setCheckable(True)
            residue_group.setChecked(False)
            residue_layout = QVBoxLayout(residue_group)
            
            # 获取该残基的所有原子
            atoms = self.Box.residues_atoms[residue_name]
            atom_checkboxes = []
            
            for atom_name in atoms:
                atom_widget = QCheckBox(atom_name)
                atom_widget.setEnabled(False)  # 默认禁用，只有选中残基时才启用
                atom_layout = QHBoxLayout()
                atom_layout.addWidget(atom_widget)
                atom_layout.addStretch()
                residue_layout.addLayout(atom_layout)
                atom_checkboxes.append(atom_widget)
            
            # 当残基被选中时，启用原子选择
            def make_residue_handler(residue_group, atom_checkboxes):
                def handler():
                    enabled = residue_group.isChecked()
                    for checkbox in atom_checkboxes:
                        checkbox.setEnabled(enabled)
                return handler
            
            residue_group.toggled.connect(make_residue_handler(residue_group, atom_checkboxes))
            
            # 存储残基组信息
            residue_group.atom_checkboxes = atom_checkboxes
            residue_group.residue_name = residue_name
            self.Box.get_config('DensityRadius').DensityGroups.append(residue_group)
            
            scroll_layout.addWidget(residue_group)
        
        scroll_area.setWidget(scroll_widget)
        scroll_area.setWidgetResizable(True)
        layout.addWidget(scroll_area)
        
        # 添加Next按钮
        btnNext = UIItemsMake.make_btn('Next')
        layout.addWidget(btnNext)
        
        # 添加到stackedWidget
        self.ui.stackedWidget_Analysis.addWidget(widget)
        self.ui.stackedWidget_Analysis.setCurrentIndex(1)
        
        # 连接按钮事件
        btnNext.clicked.connect(self.step_2_from_residue_atom)
    
    def step_2_from_residue_atom(self):
        """从残基和原子选择界面进入下一步"""
        print("DensityRadiusHandler.step_2_from_residue_atom() called")
        # 收集选中的残基和原子信息
        self._collect_density_selections()
        
        # 更新步骤
        self.current_step = 2
        
        # 创建气体组分选择界面
        self._create_gas_group_selection_layout()
        print("DensityRadiusHandler.step_2_from_residue_atom() completed")
    
    def step_2(self):
        """进入选择气体组分步骤（兼容性方法）"""
        self.step_2_from_residue_atom()
    
    def _collect_density_selections(self):
        """收集密度分析选中的残基和原子"""
        print("_collect_density_selections() called")
        density_groups = {}
        density_head_atoms = {}
        
        current_groups = self.Box.get_config('DensityRadius').DensityGroups
        print(f"Total residue groups: {len(current_groups)}")
        print(f"Type of DensityGroups: {type(current_groups)}")
        
        # 检查DensityGroups的类型
        if isinstance(current_groups, list):
            # 如果是QGroupBox对象列表
            for residue_group in current_groups:
                if hasattr(residue_group, 'residue_name') and hasattr(residue_group, 'isChecked'):
                    print(f"Checking residue group: {residue_group.residue_name}, checked: {residue_group.isChecked()}")
                    if residue_group.isChecked():
                        residue_name = residue_group.residue_name
                        selected_atoms = []
                        
                        for checkbox in residue_group.atom_checkboxes:
                            if checkbox.isChecked():
                                selected_atoms.append(checkbox.text())
                        
                        print(f"Selected atoms for {residue_name}: {selected_atoms}")
                        if selected_atoms:  # 只有当选择了原子时才记录
                            density_groups[residue_name] = selected_atoms
                            density_head_atoms[residue_name] = selected_atoms
        elif isinstance(current_groups, dict):
            # 如果已经是字典格式，直接使用
            print("DensityGroups is already a dictionary, using existing data")
            density_groups = current_groups.copy()
            density_head_atoms = self.Box.get_config('DensityRadius').DensityHeadAtoms.copy()
        
        # 保存到配置中
        self.Box.get_config('DensityRadius').DensityGroups = density_groups
        self.Box.get_config('DensityRadius').DensityHeadAtoms = density_head_atoms
        
        print(f"Final selected density groups: {density_groups}")
    
    def _create_gas_group_selection_layout(self):
        """创建气体组分选择界面"""
        print("_create_gas_group_selection_layout() called")
        widget = UIItemsMake.make_widget()
        layout = QVBoxLayout(widget)
        
        # 设置widget属性
        self.ui.widgetGasGroupDensity = widget
        self.ui.VLayoutGasGroupDensity = layout
        
        # 添加标题
        title_label = UIItemsMake.make_label('Select Gas Groups')
        title_label.setStyleSheet(f"font-size: {AnalysisBtnClick.FONT_SIZE}; font-weight: bold;")
        layout.addWidget(title_label)
        
        # 创建滚动区域
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # 为每个残基创建选择界面
        self.Box.get_config('DensityRadius').GasGroups = []
        for residue_name in self.Box.residues:
            # 创建残基组
            residue_group = QGroupBox(f"{residue_name}")
            residue_group.setCheckable(True)
            residue_group.setChecked(False)
            residue_layout = QVBoxLayout(residue_group)
            
            # 获取该残基的所有原子
            atoms = self.Box.residues_atoms[residue_name]
            atom_checkboxes = []
            
            for atom_name in atoms:
                atom_widget = QCheckBox(atom_name)
                atom_widget.setEnabled(False)  # 默认禁用，只有选中残基时才启用
                atom_layout = QHBoxLayout()
                atom_layout.addWidget(atom_widget)
                atom_layout.addStretch()
                residue_layout.addLayout(atom_layout)
                atom_checkboxes.append(atom_widget)
            
            # 当残基被选中时，启用原子选择
            def make_residue_handler(residue_group, atom_checkboxes):
                def handler():
                    enabled = residue_group.isChecked()
                    for checkbox in atom_checkboxes:
                        checkbox.setEnabled(enabled)
                return handler
            
            residue_group.toggled.connect(make_residue_handler(residue_group, atom_checkboxes))
            
            # 存储残基组信息
            residue_group.atom_checkboxes = atom_checkboxes
            residue_group.residue_name = residue_name
            self.Box.get_config('DensityRadius').GasGroups.append(residue_group)
            
            scroll_layout.addWidget(residue_group)
        
        scroll_area.setWidget(scroll_widget)
        scroll_area.setWidgetResizable(True)
        layout.addWidget(scroll_area)
        
        # 添加Run按钮（直接运行分析）
        btnRun = UIItemsMake.make_btn('Run!')
        layout.addWidget(btnRun)
        
        # 添加到stackedWidget
        self.ui.stackedWidget_Analysis.addWidget(widget)
        self.ui.stackedWidget_Analysis.setCurrentIndex(2)
        
        # 连接按钮事件
        btnRun.clicked.connect(self.run)
    
    def _collect_gas_selections(self):
        """收集气体组分选中的残基和原子"""
        gas_groups = {}
        gas_head_atoms = {}
        
        current_gas_groups = self.Box.get_config('DensityRadius').GasGroups
        print(f"Type of GasGroups: {type(current_gas_groups)}")
        
        # 检查GasGroups的类型
        if isinstance(current_gas_groups, list):
            # 如果是QGroupBox对象列表
            for residue_group in current_gas_groups:
                if hasattr(residue_group, 'residue_name') and hasattr(residue_group, 'isChecked'):
                    if residue_group.isChecked():
                        residue_name = residue_group.residue_name
                        selected_atoms = []
                        
                        for checkbox in residue_group.atom_checkboxes:
                            if checkbox.isChecked():
                                selected_atoms.append(checkbox.text())
                        
                        if selected_atoms:  # 只有当选择了原子时才记录
                            gas_groups[residue_name] = selected_atoms
                            gas_head_atoms[residue_name] = selected_atoms
        elif isinstance(current_gas_groups, dict):
            # 如果已经是字典格式，直接使用
            print("GasGroups is already a dictionary, using existing data")
            gas_groups = current_gas_groups.copy()
            gas_head_atoms = self.Box.get_config('DensityRadius').GasHeadAtoms.copy()
        
        # 保存到配置中
        self.Box.get_config('DensityRadius').GasGroups = gas_groups
        self.Box.get_config('DensityRadius').GasHeadAtoms = gas_head_atoms
        
        print(f"Selected gas groups: {gas_groups}")
    
    def run(self):
        """执行密度分析"""
        # 记录开始时间
        self.start_time = time.time()
        
        # 收集选中的气体组分信息
        self._collect_gas_selections()
        
        # 检查是否选择了气体组分
        if not self.Box.get_config('DensityRadius').GasGroups:
            create_warn_dialog("Please select at least one gas group!", "Error")
            return
        
        # 获取分析参数
        self.Info.get_text()
        
        # 添加进度条
        self._add_progress_bar()
        
        # 创建Density分析实例
        self.cls = Density(self.Box.u
                           , self.Box.get_config('DensityRadius').DensityHeadAtoms
                           , self.Box.get_config('DensityRadius').GasHeadAtoms
                           , self.Box.get_config('DensityRadius').DensityMW
                           , self.Box.get_config('DensityRadius').DensityRadius
                           , filePath=self.Info.path_result)
        
        # 启动分析
        self._start_analysis()
    
    def _add_progress_bar(self):
        """添加进度条"""
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
    
    def _start_analysis(self):
        """启动分析"""
        worker = Worker(self.cls, self.Info.frame_first, self.Info.frame_last, self.Info.step)
        worker.progressValueChanged.connect(self.updateProgressBar)
        self.thread = Thread(target=worker.run)
        self.thread.start()
    
    def updateProgressBar(self, value):
        """更新进度条"""
        # 检查进度条是否还存在，防止重复删除导致的RuntimeError
        if not hasattr(self.ui, 'progressBar') or self.ui.progressBar is None:
            return  # 如果进度条已被删除，直接返回
        
        self.ui.progressBar.setValue(value)
        if value == 100:
            if not hasattr(self.ui, 'figureTypeWidget'):
                if hasattr(self.ui, 'btnAnalysisRun'):
                    self.ui.btnAnalysisRun.deleteLater()
                self._create_figure_type_selection()
            
            # 计算分析时间并显示成功消息
            end_time = time.time()
            elapsed_time = end_time - self.start_time
            formatted_time = time.strftime('%M:%S', time.gmtime(elapsed_time))
            success_message = (
                'Analysis Completed\n'
                f"Time taken: {formatted_time}\n"
                "The gro file and topol file were saved at:\n"
                f"{self.Info.path_result}"
            )
            create_warn_dialog(success_message, 'Analysis')
            
            # 删除进度条并设置为None，防止重复访问
            if hasattr(self.ui, 'progressBar') and self.ui.progressBar is not None:
                self.ui.progressBar.deleteLater()
                self.ui.progressBar = None  # 设置为None，防止重复访问
    
    def _create_figure_type_selection(self):
        """创建图表类型选择界面 - 参考UnifiedAnalysisHandler的实现"""
        if not hasattr(self, 'cls') or self.cls is None:
            print(f"No analysis class available for DensityRadius")
            return
        
        # 创建图表类型选择widget
        self.ui.figureTypeWidget = QWidget()
        layout = QVBoxLayout(self.ui.figureTypeWidget)
        
        # 添加标题
        title_label = UIItemsMake.make_label('Select Figure Type')
        title_label.setStyleSheet(f"font-size: {AnalysisBtnClick.FONT_SIZE}; font-weight: bold;")
        layout.addWidget(title_label)
        
        # 获取分析类支持的图表类型 - 参考UnifiedAnalysisHandler的实现
        if hasattr(self.cls, 'supported_figure_types'):
            figure_types = self.cls.supported_figure_types
            print(f"Available figure types for {type(self.cls).__name__}: {figure_types}")
        else:
            # 默认图表类型
            figure_types = ['Line Chart', 'Bar Chart']
            print(f"No supported_figure_types found, using default: {figure_types}")
        
        # 创建图表类型按钮
        for figure_type in figure_types:
            btn = UIItemsMake.make_btn(f'Create {figure_type}')
            btn.clicked.connect(lambda checked, ft=figure_type: self._create_figure(ft))
            layout.addWidget(btn)
        
        # 添加到主布局
        self.ui.VLayoutRightMain.addWidget(self.ui.figureTypeWidget)

    def _create_figure(self, figure_type):
        """创建图表 - 参考UnifiedAnalysisHandler的实现"""
        if hasattr(self, 'cls') and self.cls is not None:
            if figure_type == 'Line Chart':
                self.cls.plot_line()
            elif figure_type == 'Bar Chart':
                self.cls.plot_bar()
            elif figure_type == 'Scatter Chart':
                self.cls.plot_scatter()
            elif figure_type == 'Heatmap':
                self.cls.plot_heatmap()
            elif figure_type == '3D Surface':
                self.cls.plot_3d_surface()
            else:
                print(f"Unsupported figure type: {figure_type}")
        else:
            print(f"Analysis class not available for DensityRadius")
    
    def makeFigure(self):
        """创建图表"""
        self.window = AnalysisFigureWidget(self.ui
                                           , residues=self.Box.residue_click
                                           , runMethod=self.Info.method_analysis
                                           , resultPath=self.Info.path_result)
        self.window.show()
    
    def go_back(self):
        """返回上一步"""
        print(f"DensityRadiusHandler.go_back() called, current_step: {self.current_step}")
        if self.current_step > 0:
            # 在返回前保存当前步骤的选择
            if self.current_step == 1:
                # 保存残基和原子选择
                self._save_current_density_selections()
            elif self.current_step == 2:
                # 保存气体组分选择
                self._save_current_gas_selections()
            
            self.current_step -= 1
            if self.current_step == 0:
                # 返回第一步：MW和Radius输入
                self.ui.stackedWidget_Analysis.setCurrentIndex(0)
            elif self.current_step == 1:
                # 返回第二步：残基和原子选择
                self.ui.stackedWidget_Analysis.setCurrentIndex(1)
            elif self.current_step == 2:
                # 返回第三步：气体组分选择
                self.ui.stackedWidget_Analysis.setCurrentIndex(2)
        else:
            print("Already at first step, cannot go back")
    
    
    def _create_mw_radius_layout(self):
        """创建MW和Radius输入界面（用于refresh）"""
        # 清空当前界面
        current_widget = self.ui.stackedWidget_Analysis.widget(0)
        if current_widget:
            current_widget.deleteLater()
        
        # 重新创建界面
        widget = UIItemsMake.make_widget()
        layout = QVBoxLayout(widget)
        self.ui.widgetspinbox = widget
        self.ui.VLayoutspinbox = layout
        
        # 创建MW SpinBox
        label_mw = UIItemsMake.make_label('MW(g/mol)')
        spin_box_mw = UIItemsMake.make_spin_box(14, 0, 1000000)
        self.ui.widgetspinboxLabel0 = label_mw
        self.ui.widgetspinboxSpinBox0 = spin_box_mw
        layout.addWidget(label_mw)
        layout.addWidget(spin_box_mw)
        
        # 创建Radius SpinBox
        label_radius = UIItemsMake.make_label('Radius(A)')
        spin_box_radius = UIItemsMake.make_spin_box(50, 0, 1000000)
        self.ui.widgetspinboxLabel1 = label_radius
        self.ui.widgetspinboxSpinBox1 = spin_box_radius
        layout.addWidget(label_radius)
        layout.addWidget(spin_box_radius)
        
        # 创建Next按钮
        btnNext = UIItemsMake.make_btn('Next')
        layout.addWidget(btnNext)
        
        # 替换第一个widget
        self.ui.stackedWidget_Analysis.insertWidget(0, widget)
        self.ui.stackedWidget_Analysis.setCurrentIndex(0)
        
        # 重新连接按钮事件
        btnNext.clicked.connect(self.step_1_from_mw_radius)
    
    def _save_current_density_selections(self):
        """保存当前残基和原子选择"""
        print("_save_current_density_selections() called")
        density_groups = {}
        density_head_atoms = {}
        
        current_groups = self.Box.get_config('DensityRadius').DensityGroups
        if isinstance(current_groups, list):
            # 如果是QGroupBox对象列表，收集选择
            for residue_group in current_groups:
                if hasattr(residue_group, 'residue_name') and hasattr(residue_group, 'isChecked'):
                    if residue_group.isChecked():
                        residue_name = residue_group.residue_name
                        selected_atoms = []
                        
                        for checkbox in residue_group.atom_checkboxes:
                            if checkbox.isChecked():
                                selected_atoms.append(checkbox.text())
                        
                        if selected_atoms:
                            density_groups[residue_name] = selected_atoms
                            density_head_atoms[residue_name] = selected_atoms
        
        # 保存到配置中
        self.Box.get_config('DensityRadius').DensityGroups = density_groups
        self.Box.get_config('DensityRadius').DensityHeadAtoms = density_head_atoms
        
        print(f"Saved density groups: {density_groups}")
    
    def _save_current_gas_selections(self):
        """保存当前气体组分选择"""
        print("_save_current_gas_selections() called")
        gas_groups = {}
        gas_head_atoms = {}
        
        current_gas_groups = self.Box.get_config('DensityRadius').GasGroups
        if isinstance(current_gas_groups, list):
            # 如果是QGroupBox对象列表，收集选择
            for residue_group in current_gas_groups:
                if hasattr(residue_group, 'residue_name') and hasattr(residue_group, 'isChecked'):
                    if residue_group.isChecked():
                        residue_name = residue_group.residue_name
                        selected_atoms = []
                        
                        for checkbox in residue_group.atom_checkboxes:
                            if checkbox.isChecked():
                                selected_atoms.append(checkbox.text())
                        
                        if selected_atoms:
                            gas_groups[residue_name] = selected_atoms
                            gas_head_atoms[residue_name] = selected_atoms
        
        # 保存到配置中
        self.Box.get_config('DensityRadius').GasGroups = gas_groups
        self.Box.get_config('DensityRadius').GasHeadAtoms = gas_head_atoms
        
        print(f"Saved gas groups: {gas_groups}")
    
    def _refresh_mw_radius_layout(self):
        """刷新MW和Radius输入界面"""
        print("_refresh_mw_radius_layout() called")
        self._create_mw_radius_layout()
    
    def _refresh_residue_atom_layout(self):
        """刷新残基和原子选择界面"""
        print("_refresh_residue_atom_layout() called")
        # 清空当前界面
        current_widget = self.ui.stackedWidget_Analysis.widget(1)
        if current_widget:
            current_widget.deleteLater()
        
        # 重新创建界面
        self._create_residue_atom_selection_layout()
    
    def _refresh_gas_group_layout(self):
        """刷新气体组分选择界面"""
        print("_refresh_gas_group_layout() called")
        # 清空当前界面
        current_widget = self.ui.stackedWidget_Analysis.widget(2)
        if current_widget:
            current_widget.deleteLater()
        
        # 重新创建界面
        self._create_gas_group_selection_layout()


class DensityMultiRadiusHandler:
    """专门处理DensityMultiRadius分析流程的类，不继承AtomsLayout"""
    
    def __init__(self, ui, Box, Info):
        self.ui = ui
        self.Box = Box
        self.Info = Info
        self.current_step = 0
        self.step_widgets = []
    
    def start(self):
        """开始DensityMultiRadius分析流程"""
        print("DensityMultiRadiusHandler.start() called")
        self.current_step = 0
        self._create_mw_maxradius_segments_layout()
    
    def step_1_from_mw_maxradius_segments(self):
        """从MW、MaxRadius、NumberSegments输入进入下一步"""
        print("DensityMultiRadiusHandler.step_1_from_mw_maxradius_segments() called")
        
        # 获取用户输入的参数
        mw = self.ui.widgetspinboxSpinBox0.value()
        max_radius = self.ui.widgetspinboxSpinBox1.value()
        number_segments = self.ui.widgetspinboxSpinBox2.value()
        
        print(f"MW: {mw}, MaxRadius: {max_radius}, NumberSegments: {number_segments}")
        
        # 保存参数到配置
        config = self.Box.get_config('DensityMultiRadius')
        config.DensityMW = mw
        config.DensityMaxRadius = max_radius
        config.DensityNumberSegments = number_segments
        
        # 进入下一步
        self.current_step = 1
        self._create_residue_atom_selection_layout()
    
    def step_1(self):
        """兼容性方法"""
        self.step_1_from_mw_maxradius_segments()
    
    def _create_mw_maxradius_segments_layout(self):
        """创建MW、MaxRadius、NumberSegments输入界面"""
        print("_create_mw_maxradius_segments_layout() called")
        
        # 创建界面
        widget = UIItemsMake.make_widget()
        layout = QVBoxLayout(widget)
        
        # 创建MW SpinBox
        label_mw = UIItemsMake.make_label('MW(g/mol)')
        spin_box_mw = UIItemsMake.make_spin_box(14, 0, 1000000)
        layout.addWidget(label_mw)
        layout.addWidget(spin_box_mw)
        
        # 创建MaxRadius SpinBox
        label_max_radius = UIItemsMake.make_label('MaxRadius(A)')
        spin_box_max_radius = UIItemsMake.make_spin_box(50, 0, 1000000)
        layout.addWidget(label_max_radius)
        layout.addWidget(spin_box_max_radius)
        
        # 创建NumberSegments SpinBox
        label_number_segments = UIItemsMake.make_label('NumberSegments')
        spin_box_number_segments = UIItemsMake.make_spin_box(5, 1, 1000)
        layout.addWidget(label_number_segments)
        layout.addWidget(spin_box_number_segments)
        
        # 创建Next按钮
        btnNext = UIItemsMake.make_btn('Next')
        layout.addWidget(btnNext)
        
        # 保存UI组件引用
        self.ui.widgetspinbox = widget
        self.ui.VLayoutspinbox = layout
        self.ui.widgetspinboxLabel0 = label_mw
        self.ui.widgetspinboxSpinBox0 = spin_box_mw
        self.ui.widgetspinboxLabel1 = label_max_radius
        self.ui.widgetspinboxSpinBox1 = spin_box_max_radius
        self.ui.widgetspinboxLabel2 = label_number_segments
        self.ui.widgetspinboxSpinBox2 = spin_box_number_segments
        
        # 添加到stackedWidget
        self.ui.stackedWidget_Analysis.addWidget(widget)
        self.ui.stackedWidget_Analysis.setCurrentIndex(0)
        
        # 连接按钮事件
        btnNext.clicked.connect(lambda: self.step_1_from_mw_maxradius_segments())
    
    def _create_residue_atom_selection_layout(self):
        """创建残基和原子选择界面"""
        widget = UIItemsMake.make_widget()
        layout = QVBoxLayout(widget)
        
        # 设置widget属性
        self.ui.widgetResidueAtomDensityMulti = widget
        self.ui.VLayoutResidueAtomDensityMulti = layout
        
        # 添加标题
        title_label = UIItemsMake.make_label('Select Residues Group')
        title_label.setStyleSheet(f"font-size: {AnalysisBtnClick.FONT_SIZE}; font-weight: bold;")
        layout.addWidget(title_label)
        
        # 创建滚动区域
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # 为每个残基创建选择界面
        self.Box.get_config('DensityMultiRadius').DensityGroups = []
        for residue_name in self.Box.residues:
            # 创建残基组
            residue_group = QGroupBox(f"{residue_name}")
            residue_group.setCheckable(True)
            residue_group.setChecked(False)
            residue_layout = QVBoxLayout(residue_group)
            
            # 获取该残基的所有原子
            atoms = self.Box.residues_atoms[residue_name]
            atom_checkboxes = []
            
            for atom_name in atoms:
                atom_widget = QCheckBox(atom_name)
                atom_widget.setEnabled(False)  # 默认禁用，只有选中残基时才启用
                atom_layout = QHBoxLayout()
                atom_layout.addWidget(atom_widget)
                atom_layout.addStretch()
                residue_layout.addLayout(atom_layout)
                atom_checkboxes.append(atom_widget)
            
            # 当残基被选中时，启用原子选择
            def make_residue_handler(residue_group, atom_checkboxes):
                def handler():
                    enabled = residue_group.isChecked()
                    for checkbox in atom_checkboxes:
                        checkbox.setEnabled(enabled)
                return handler
            
            residue_group.toggled.connect(make_residue_handler(residue_group, atom_checkboxes))
            
            # 存储残基组信息
            residue_group.atom_checkboxes = atom_checkboxes
            residue_group.residue_name = residue_name
            self.Box.get_config('DensityMultiRadius').DensityGroups.append(residue_group)
            
            scroll_layout.addWidget(residue_group)
        
        scroll_area.setWidget(scroll_widget)
        scroll_area.setWidgetResizable(True)
        layout.addWidget(scroll_area)
        
        # 添加Next按钮
        btnNext = UIItemsMake.make_btn('Next')
        layout.addWidget(btnNext)
        
        # 添加到stackedWidget
        self.ui.stackedWidget_Analysis.addWidget(widget)
        self.ui.stackedWidget_Analysis.setCurrentIndex(1)
        
        # 连接按钮事件
        btnNext.clicked.connect(self.step_2_from_residue_atom)
    
    def step_2_from_residue_atom(self):
        """从残基和原子选择进入下一步"""
        print("DensityMultiRadiusHandler.step_2_from_residue_atom() called")
        
        # 收集选择
        self._collect_density_selections()
        
        # 进入下一步
        self.current_step = 2
        self._create_gas_group_selection_layout()
    
    def step_2(self):
        """兼容性方法"""
        self.step_2_from_residue_atom()
    
    def _collect_density_selections(self):
        """收集密度分析选中的残基和原子"""
        print("_collect_density_selections() called")
        density_groups = {}
        density_head_atoms = {}
        
        current_groups = self.Box.get_config('DensityMultiRadius').DensityGroups
        print(f"Total residue groups: {len(current_groups)}")
        print(f"Type of DensityGroups: {type(current_groups)}")
        
        # 检查DensityGroups的类型
        if isinstance(current_groups, list):
            # 如果是QGroupBox列表，收集选择
            for residue_group in current_groups:
                if residue_group.isChecked():
                    residue_name = residue_group.residue_name
                    selected_atoms = []
                    for atom_checkbox in residue_group.atom_checkboxes:
                        if atom_checkbox.isChecked():
                            selected_atoms.append(atom_checkbox.text())
                    
                    if selected_atoms:
                        density_groups[residue_name] = selected_atoms
                        density_head_atoms[residue_name] = selected_atoms
        elif isinstance(current_groups, dict):
            # 如果已经是字典格式，直接使用
            print("DensityGroups is already a dictionary, using existing data")
            density_groups = current_groups.copy()
            density_head_atoms = self.Box.get_config('DensityMultiRadius').DensityHeadAtoms.copy()
        
        # 保存到配置中
        self.Box.get_config('DensityMultiRadius').DensityGroups = density_groups
        self.Box.get_config('DensityMultiRadius').DensityHeadAtoms = density_head_atoms
        
        print(f"Final selected density groups: {density_groups}")
    
    def _create_gas_group_selection_layout(self):
        """创建气体组分选择界面"""
        widget = UIItemsMake.make_widget()
        layout = QVBoxLayout(widget)
        
        # 设置widget属性
        self.ui.widgetGasGroupDensityMulti = widget
        self.ui.VLayoutGasGroupDensityMulti = layout
        
        # 添加标题
        title_label = UIItemsMake.make_label('Select Gas group')
        title_label.setStyleSheet(f"font-size: {AnalysisBtnClick.FONT_SIZE}; font-weight: bold;")
        layout.addWidget(title_label)
        
        # 创建滚动区域
        scroll_area = QScrollArea()
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        
        # 为每个残基创建选择界面
        self.Box.get_config('DensityMultiRadius').GasGroups = []
        for residue_name in self.Box.residues:
            # 创建残基组
            residue_group = QGroupBox(f"{residue_name}")
            residue_group.setCheckable(True)
            residue_group.setChecked(False)
            residue_layout = QVBoxLayout(residue_group)
            
            # 获取该残基的所有原子
            atoms = self.Box.residues_atoms[residue_name]
            atom_checkboxes = []
            
            for atom_name in atoms:
                atom_widget = QCheckBox(atom_name)
                atom_widget.setEnabled(False)  # 默认禁用，只有选中残基时才启用
                atom_layout = QHBoxLayout()
                atom_layout.addWidget(atom_widget)
                atom_layout.addStretch()
                residue_layout.addLayout(atom_layout)
                atom_checkboxes.append(atom_widget)
            
            # 当残基被选中时，启用原子选择
            def make_residue_handler(residue_group, atom_checkboxes):
                def handler():
                    enabled = residue_group.isChecked()
                    for checkbox in atom_checkboxes:
                        checkbox.setEnabled(enabled)
                return handler
            
            residue_group.toggled.connect(make_residue_handler(residue_group, atom_checkboxes))
            
            # 存储残基组信息
            residue_group.atom_checkboxes = atom_checkboxes
            residue_group.residue_name = residue_name
            self.Box.get_config('DensityMultiRadius').GasGroups.append(residue_group)
            
            scroll_layout.addWidget(residue_group)
        
        scroll_area.setWidget(scroll_widget)
        scroll_area.setWidgetResizable(True)
        layout.addWidget(scroll_area)
        
        # 添加Run按钮（直接运行分析）
        btnRun = UIItemsMake.make_btn('Run!')
        layout.addWidget(btnRun)
        
        # 添加到stackedWidget
        self.ui.stackedWidget_Analysis.addWidget(widget)
        self.ui.stackedWidget_Analysis.setCurrentIndex(2)
        
        # 连接按钮事件
        btnRun.clicked.connect(self.run)
    
    def _collect_gas_selections(self):
        """收集选中的气体组分信息"""
        print("_collect_gas_selections() called")
        gas_groups = {}
        gas_head_atoms = {}
        
        current_groups = self.Box.get_config('DensityMultiRadius').GasGroups
        print(f"Total gas groups: {len(current_groups)}")
        print(f"Type of GasGroups: {type(current_groups)}")
        
        # 检查GasGroups的类型
        if isinstance(current_groups, list):
            # 如果是QGroupBox列表，收集选择
            for residue_group in current_groups:
                if residue_group.isChecked():
                    residue_name = residue_group.residue_name
                    selected_atoms = []
                    for atom_checkbox in residue_group.atom_checkboxes:
                        if atom_checkbox.isChecked():
                            selected_atoms.append(atom_checkbox.text())
                    
                    if selected_atoms:
                        gas_groups[residue_name] = selected_atoms
                        gas_head_atoms[residue_name] = selected_atoms
        elif isinstance(current_groups, dict):
            # 如果已经是字典格式，直接使用
            print("GasGroups is already a dictionary, using existing data")
            gas_groups = current_groups.copy()
            gas_head_atoms = self.Box.get_config('DensityMultiRadius').GasHeadAtoms.copy()
        
        # 保存到配置中
        self.Box.get_config('DensityMultiRadius').GasGroups = gas_groups
        self.Box.get_config('DensityMultiRadius').GasHeadAtoms = gas_head_atoms
        
        print(f"Final selected gas groups: {gas_groups}")
    
    def run(self):
        """执行密度分析"""
        # 记录开始时间
        self.start_time = time.time()
        
        # 收集选中的气体组分信息
        self._collect_gas_selections()
        
        # 检查是否选择了气体组分
        if not self.Box.get_config('DensityMultiRadius').GasGroups:
            create_warn_dialog("Please select at least one gas group!", "Error")
            return
        
        # 获取分析参数
        self.Info.get_text()
        
        # 添加进度条
        self._add_progress_bar()
        
        # 创建DensityMultiRadius分析实例
        self.cls = DensityMultiRadius(self.Box.u
                           , self.Box.get_config('DensityMultiRadius').DensityHeadAtoms
                           , self.Box.get_config('DensityMultiRadius').GasHeadAtoms
                           , self.Box.get_config('DensityMultiRadius').DensityMW
                           , self.Box.get_config('DensityMultiRadius').DensityMaxRadius
                           , self.Box.get_config('DensityMultiRadius').DensityNumberSegments
                           , filePath=self.Info.path_result)
        
        # 启动分析
        self._start_analysis()
    
    def _add_progress_bar(self):
        """添加进度条"""
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
    
    def _start_analysis(self):
        """启动分析"""
        worker = Worker(self.cls, self.Info.frame_first, self.Info.frame_last, self.Info.step)
        worker.progressValueChanged.connect(self.updateProgressBar)
        self.thread = Thread(target=worker.run)
        self.thread.start()
    
    def updateProgressBar(self, value):
        """更新进度条"""
        # 检查进度条是否还存在，防止重复删除导致的RuntimeError
        if not hasattr(self.ui, 'progressBar') or self.ui.progressBar is None:
            return  # 如果进度条已被删除，直接返回
        
        self.ui.progressBar.setValue(value)
        if value == 100:
            if not hasattr(self.ui, 'figureTypeWidget'):
                if hasattr(self.ui, 'btnAnalysisRun'):
                    self.ui.btnAnalysisRun.deleteLater()
                self._create_figure_type_selection()
            
            # 计算分析时间并显示成功消息
            end_time = time.time()
            elapsed_time = end_time - self.start_time
            formatted_time = time.strftime('%M:%S', time.gmtime(elapsed_time))
            success_message = (
                'Analysis Completed\n'
                f"Time taken: {formatted_time}\n"
                "The gro file and topol file were saved at:\n"
                f"{self.Info.path_result}"
            )
            create_warn_dialog(success_message, 'Analysis')
            
            # 删除进度条并设置为None，防止重复访问
            if hasattr(self.ui, 'progressBar') and self.ui.progressBar is not None:
                self.ui.progressBar.deleteLater()
                self.ui.progressBar = None  # 设置为None，防止重复访问
    
    def _create_figure_type_selection(self):
        """创建图表类型选择界面 - 参考UnifiedAnalysisHandler的实现"""
        if not hasattr(self, 'cls') or self.cls is None:
            print(f"No analysis class available for DensityMultiRadius")
            return
        
        # 创建图表类型选择widget
        self.ui.figureTypeWidget = QWidget()
        layout = QVBoxLayout(self.ui.figureTypeWidget)
        
        # 添加标题
        title_label = UIItemsMake.make_label('Select Figure Type')
        title_label.setStyleSheet(f"font-size: {AnalysisBtnClick.FONT_SIZE}; font-weight: bold;")
        layout.addWidget(title_label)
        
        # 获取分析类支持的图表类型 - 参考UnifiedAnalysisHandler的实现
        if hasattr(self.cls, 'supported_figure_types'):
            figure_types = self.cls.supported_figure_types
            print(f"Available figure types for {type(self.cls).__name__}: {figure_types}")
        else:
            # 默认图表类型
            figure_types = ['Line Chart', 'Bar Chart']
            print(f"No supported_figure_types found, using default: {figure_types}")
        
        # 创建图表类型按钮
        for figure_type in figure_types:
            btn = UIItemsMake.make_btn(f'Create {figure_type}')
            btn.clicked.connect(lambda checked, ft=figure_type: self._create_figure(ft))
            layout.addWidget(btn)
        
        # 添加到主布局
        self.ui.VLayoutRightMain.addWidget(self.ui.figureTypeWidget)

    def _create_figure(self, figure_type):
        """创建图表 - 参考UnifiedAnalysisHandler的实现"""
        if hasattr(self, 'cls') and self.cls is not None:
            if figure_type == 'Line Chart':
                self.cls.plot_line()
            elif figure_type == 'Bar Chart':
                self.cls.plot_bar()
            elif figure_type == 'Scatter Chart':
                self.cls.plot_scatter()
            elif figure_type == 'Heatmap':
                self.cls.plot_heatmap()
            elif figure_type == '3D Surface':
                self.cls.plot_3d_surface()
            else:
                print(f"Unsupported figure type: {figure_type}")
        else:
            print(f"Analysis class not available for DensityMultiRadius")
    
    def _make_figure(self):
        """创建图表"""
        if hasattr(self, 'cls') and self.cls:
            self.cls.make_figure()
    
    def go_back(self):
        """返回上一步"""
        print(f"DensityMultiRadiusHandler.go_back() called, current_step: {self.current_step}")
        
        if self.current_step > 0:
            # 保存当前步骤的选择
            if self.current_step == 1:
                self._save_current_density_selections()
            elif self.current_step == 2:
                self._save_current_gas_selections()
            
            self.current_step -= 1
            
            if self.current_step == 0:
                # 返回第一步：MW、MaxRadius、NumberSegments输入
                self.ui.stackedWidget_Analysis.setCurrentIndex(0)
            elif self.current_step == 1:
                # 返回第二步：残基和原子选择
                self.ui.stackedWidget_Analysis.setCurrentIndex(1)
        else:
            print("Already at first step, cannot go back")
    
    def _save_current_density_selections(self):
        """保存当前残基和原子选择"""
        print("_save_current_density_selections() called")
        
        config = self.Box.get_config('DensityMultiRadius')
        
        if isinstance(config.DensityGroups, list):
            # 收集当前选择
            density_selections = {}
            for group_box in config.DensityGroups:
                group_data = group_box.group_data
                residue_name = group_data['residue_name']
                residue_checkbox = group_data['residue_checkbox']
                atom_checkboxes = group_data['atom_checkboxes']
                atoms = group_data['atoms']
                
                if residue_checkbox.isChecked():
                    selected_atoms = []
                    for i, atom_checkbox in enumerate(atom_checkboxes):
                        if atom_checkbox.isChecked():
                            selected_atoms.append(atoms[i])
                    
                    if selected_atoms:
                        density_selections[residue_name] = selected_atoms
            
            config.DensityGroups = density_selections
            print(f"Saved density selections: {density_selections}")
    
    def _save_current_gas_selections(self):
        """保存当前气体组分选择"""
        print("_save_current_gas_selections() called")
        
        config = self.Box.get_config('DensityMultiRadius')
        
        if isinstance(config.GasGroups, list):
            # 收集当前选择
            gas_selections = {}
            for group_box in config.GasGroups:
                group_data = group_box.group_data
                residue_name = group_data['residue_name']
                residue_checkbox = group_data['residue_checkbox']
                atom_checkboxes = group_data['atom_checkboxes']
                atoms = group_data['atoms']
                
                if residue_checkbox.isChecked():
                    selected_atoms = []
                    for i, atom_checkbox in enumerate(atom_checkboxes):
                        if atom_checkbox.isChecked():
                            selected_atoms.append(atoms[i])
                    
                    if selected_atoms:
                        gas_selections[residue_name] = selected_atoms
            
            config.GasGroups = gas_selections
            print(f"Saved gas selections: {gas_selections}")
    



# class PCALayout(AnalysisLayout):
#     def step_1(self):
#         self.AtomsLayout('widgetAtomsPCA'
#                          , 'VLayoutAtomsPCA'
#                          , 'select head atom'
#                          , self.Box.get_config('PCA').PCAHeadGroups
#                          , QRadioButton
#                          , 'Run!'
#                          , 1
#                          , self.run
#                          , connect=True)

#     @AnalysisLayout._addProgressBar
#     def run(self):
#         self.Info.get_text()
#         AnalysisUtils.get_atom(self.Box.get_config('PCA').PCAHeadGroups
#                                , QRadioButton
#                                , self.Box.get_config('PCA').PCAHeadAtoms)
#         self.cls = PCA(self.Box.u
#                               , self.Box.get_config('PCA').PCAHeadAtoms
#                               , filePath=self.Info.path_result)


# class RDLayout(AnalysisLayout):
#     def step_1(self):
#         self.AtomsLayout('widgetAtomsRD'
#                          , 'VLayoutAtomsRD'
#                          , 'select head atom'
#                          , self.Box.get_config('RadialDistribution').RDHeadGroups
#                          , QRadioButton
#                          , 'Run!'
#                          , 1
#                          , self.run
#                          , connect=True)

#     @AnalysisLayout._addProgressBar
#     def run(self):
#         self.Info.get_text()
#         AnalysisUtils.get_atom(self.Box.get_config('RadialDistribution').RDHeadGroups
#                                , QRadioButton
#                                , self.Box.get_config('RadialDistribution').RDHeadAtoms)
#         self.cls = CalRad(self.Box.u
#                           , self.Box.get_config('RadialDistribution').RDHeadAtoms
#                           , filePath=self.Info.path_result)


# class PressureLayout(AnalysisLayout):
#     def step_1(self):
#         self.AtomsLayout('widgetBallResiduesPressure'
#                          , 'VLayoutBallResiduesPressure'
#                          , 'select head atom'
#                          , self.Box.get_config('Pressure').PressureHeadGroups
#                          , QRadioButton
#                          , 'Next!'
#                          , 1
#                          , func=self.step_2
#                          , connect=True)

#     def step_2(self):
#         AnalysisUtils.get_atom(self.Box.get_config('Pressure').PressureHeadGroups
#                                , QRadioButton
#                                , self.Box.get_config('Pressure').PressureHeadAtoms)
#         self.AtomsLayout('widgetGasResidues'
#                          , 'VLayoutGasResidues'
#                          , 'select gas residue'
#                          , self.Box.get_config('Pressure').GasGroups
#                          , QRadioButton
#                          , 'Next!'
#                          , 2
#                          , func=self.step_3
#                          , connect=True)

#     def step_3(self):

#         self.AtomsLayout('widgetGasAtoms'
#                          , 'VLayoutGasAtoms'
#                          , 'select gas atom'
#                          , self.Box.get_config('Pressure').GasHeadAtoms
#                          , QRadioButton
#                          , 'Run!'
#                          , 3
#                          , self.step_4
#                          , connect=True)

#     def step_4(self):
#         AnalysisUtils.get_atom(self.Box.get_config('Pressure').GasGroups
#                                , QRadioButton
#                                , self.Box.get_config('Pressure').GasHeadAtoms)
#         self.SpinLayout('widgetNcircles'
#                         , 'VLayoutNcircles'
#                         , 'select num of circle'
#                         , 10
#                         , 0
#                         , 1000000
#                         , 'RUN!'
#                         , 4
#                         , self.run
#                         , connect=True)

#     @AnalysisLayout._addProgressBar
#     def run(self):
#         self.Info.get_text()
#         self.cls = Pressure(self.Box.u
#                             , self.Box.get_config('Pressure').PressureHeadAtoms
#                             , self.Box.get_config('Pressure').GasHeadAtoms
#                             , n_circles=AnalysisUtils.get_spin_value(self.ui.widgetNcirclesSpinBox)
#                             , filePath=self.Info.path_result)


# class MeanCurvatureLayout(AnalysisLayout):
#     def step_1(self):
#         self.AtomsLayout('widgetAtomsMeanCurvature'
#                          , 'VLayoutAtomsMeanCurvature'
#                          , 'select head atom'
#                          , self.Box.get_config('MeanCurvature').MCHeadGroups
#                          , QRadioButton
#                          , 'Run!'
#                          , 1
#                          , self.run
#                          , connect=True)

#     @AnalysisLayout._addProgressBar
#     def run(self):
#         self.Info.get_text()
#         AnalysisUtils.get_atom(self.Box.get_config('MeanCurvature').MCHeadGroups
#                                , QRadioButton
#                                , self.Box.get_config('MeanCurvature').MeanCurvatureHeadAtoms
#                                )
#         self.cls = Curvature(self.Box.u
#                              , self.Box.get_config('MeanCurvature').MeanCurvatureHeadAtoms
#                              , filePath=self.Info.path_result
#                              , k=self.Info.K
#                              , method='mean')


# class NClusterLayout(AnalysisLayout):
#     def step_1(self):
#         self.AtomsLayout('widgetAtomsNCluster'
#                          , 'VLayoutAtomsNCluster'
#                          , 'select atoms'
#                          , self.Box.get_config('NCluster').NCLGroups
#                          , QCheckBox
#                          , 'Next'
#                          , 1
#                          , self.step_2)

#     def step_2(self):
#         self.SpinLayout('widgetNCutoff'
#                         , 'VLayoutNCutoff'
#                         , ['select cutoff value(A)', 'select cutoff number']
#                         , [12, 10]
#                         , [0, 0]
#                         , [1000000, 1000000]
#                         , 'RUN!'
#                         , 2
#                         , self.run
#                         , connect=True
#                         , num=2)

#     @AnalysisLayout._addProgressBar
#     def run(self):
#         self.Info.get_text()
#         AnalysisUtils.get_atom(self.Box.get_config('NCluster').NCLGroups
#                                , QCheckBox
#                                , self.Box.get_config('NCluster').NCLHeadAtoms)
#         self.Box.get_config('NCluster').NCLCutoff = AnalysisUtils.get_spin_value(self.ui.widgetNCutoffSpinBox0)
#         self.Box.get_config('NCluster').NCutoff = AnalysisUtils.get_spin_value(self.ui.widgetNCutoffSpinBox1)
#         self.cls = NCluster(self.Box.u
#                            , self.Box.get_config('NCluster').NCLHeadAtoms
#                            , filePath=self.Info.path_result
#                            , cutoff=self.Box.get_config('NCluster').NCLCutoff
#                            , N_cutoff=self.Box.get_config('NCluster').NCutoff)


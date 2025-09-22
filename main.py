# ///////////////////////////////////////////////////////////////
# LNB-Analysis
# BY: RenXinYU
# PROJECT MADE WITH: Qt Designer and PySide6
# V: 1.0.0
#
#
# ///////////////////////////////////////////////////////////////
import time
import sys
import os
import platform
import configparser
from modules.vmd_control import VMDCommands, VMDTcp
# IMPORT / GUI AND MODULES AND WIDGETS
# ///////////////////////////////////////////////////////////////

from modules import *
from widgets import *
from figure import *
from analysis import *
from modules.file_detection import FileTypeDetector, DynamicTabManager
from modules.theme_manager import ThemeManager, add_theme_button_to_top_menu

os.environ["QT_FONT_DPI"] = "96"  # FIX Problem for High DPI and Scale above 100%

# SET AS GLOBAL WIDGETS
# ///////////////////////////////////////////////////////////////
widgets = None

def load_config():
    """加载配置文件"""
    config = configparser.ConfigParser()
    config_file = "config.ini"
    
    # 如果配置文件不存在，创建默认配置
    if not os.path.exists(config_file):
        create_default_config(config_file)
    
    config.read(config_file, encoding='utf-8')
    return config

def create_default_config(config_file):
    """创建默认配置文件"""
    default_config = """[VMD]
# VMD路径配置
# Windows用户请修改为您的VMD安装路径
# macOS用户请修改为您的VMD安装路径
# Linux用户请修改为您的VMD安装路径

# Windows示例路径:
# vmd_path = C:/Program Files/VMD/vmd.exe
# vmd_path = C:/VMD/vmd.exe

# macOS示例路径:
# vmd_path = /Applications/VMD.app/Contents/MacOS/VMD
# vmd_path = /usr/local/bin/vmd

# Linux示例路径:
# vmd_path = /usr/local/bin/vmd
# vmd_path = /opt/vmd/vmd

# 默认路径 (请根据您的系统修改)
vmd_path = C:/Program Files/VMD/vmd.exe

[Analysis]
# 分析模块默认配置
default_parallel = true
default_n_jobs = -1
default_k_value = 15

[UI]
# 界面配置
theme = dark
language = zh_CN
"""
    with open(config_file, 'w', encoding='utf-8') as f:
        f.write(default_config)

def get_vmd_path():
    """获取VMD路径"""
    config = load_config()
    vmd_path = config.get('VMD', 'vmd_path', fallback='C:/Program Files/VMD/vmd.exe')
    
    # 检查路径是否存在
    if not os.path.exists(vmd_path):
        print(f"警告: VMD路径不存在: {vmd_path}")
        print("请在 config.ini 文件中修改正确的VMD路径")
        return None
    
    return vmd_path


class MainWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)

        # SET AS GLOBAL WIDGETS
        # ///////////////////////////////////////////////////////////////
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
    # +++ 添加以下代码以解决背景颜色问题 +++
        # 强制所有QPushButton都填充自己的背景，以覆盖系统原生样式。
        # 这是解决macOS等系统上背景色不生效问题的关键。
        # for button in self.findChildren(QPushButton):
        #     button.setAutoFillBackground(True)

        # 初始化 VMDControlPanel
        global widgets
        widgets = self.ui
        # 初始化 VMD 相关
        self.rctl_path = "./remote_ctl.tcl"
        self.vmd_path = get_vmd_path()  # 从配置文件读取VMD路径
        self.vmd = None
        self.connected = False
        self.data = None
        self.valid_comments = None

        # 设置表格选择模式
        self.ui.vmd_tablewidget.setSelectionMode(QTableWidget.ExtendedSelection)

        # 启用拖放
        self.setAcceptDrops(True)

        # 绑定事件
        self.ui.vmd_btn_start.clicked.connect(self.pushStartVMD)
        self.ui.vmd_btn_stop.clicked.connect(self.pushStopVMD)
        self.ui.vmd_tablewidget.selectionModel().selectionChanged.connect(self.onSelectionChanged)

        # 初始化 UI 状态
        self.ui.vmd_btn_stop.setEnabled(False)
        self.ui.vmd_label.setText("Click 'Start VMD' to launch VMD, then drag and drop a CSV file")
        # USE CUSTOM TITLE BAR | USE AS "False" FOR MAC OR LINUX
        # ///////////////////////////////////////////////////////////////
        Settings.ENABLE_CUSTOM_TITLE_BAR = True

        # TOGGLE MENU
        # ///////////////////////////////////////////////////////////////
        widgets.toggleButton.clicked.connect(lambda: UIFunctions.toggleMenu(self, True))

        # SET UI DEFINITIONS
        UIFunctions.uiDefinitions(self)


        # LEFT MENUS
        widgets.btn_home.clicked.connect(self.buttonLeftClick)
        widgets.btn_generate.clicked.connect(self.buttonLeftClick)
        widgets.btn_figure.clicked.connect(self.buttonLeftClick)
        widgets.btn_analysis.clicked.connect(self.buttonLeftClick)
        widgets.btn_data_process.clicked.connect(self.buttonLeftClick)
        self.initialSettings()
        
        # 初始化文件检测功能
        self.setup_file_detection()
        
        # 初始化主题切换功能
        self.setup_theme_switch()
        
        # 确保figure页面显示信息文本框
        self.ensure_figure_widget_display()

        # EXTRA RIGHT BOX

        def openCloseRightBox():
            UIFunctions.toggleRightBox(self, True)

        def openCloseGeneBox():
            UIFunctions.toggleGeneRightBox(self, True)

        def openCloseFigureColorBox():
            if self.ui.figure_extra == 1:
                UIFunctions.toggleFigureColorBox(self, True)
            elif self.ui.figure_extra == 0:
                UIFunctions.toggleFigureColorBox(self, False)

        def openCloseFigureShapeBox():
            UIFunctions.toggleFigureShapeBox(self, True)

        def openCloseFigureExtra(index):
            if index == 0 or index == 1:
                UIFunctions.toggleFigureShapeBox(self, False)
            elif index == 2:
                UIFunctions.toggleFigureColorBox(self, False)

        self.ui.btn_language.clicked.connect(lambda: AppFunctions.btnLanguageClick(self.ui))
        # //////////////////////////Generation/////////////////////////////////////
        # btnGeneration Click
        # self.ui.btn_gene_run
        self.ui.btn_gene_path.clicked.connect(lambda: BtnGetPath.run(self.ui.edit_gene_path, 'gene_gro'))
        self.ui.btn_gene_run.clicked.connect(lambda: BtnGeneClick(self.ui))

        self.ui.btn_gene_lipid.clicked.connect(openCloseGeneBox)
        self.ui.btn_gene_lipid.clicked.connect(lambda: lipidsSelect(self.ui))
        # //////////////////////////Analysis/////////////////////////////////////
        # btnTools CLICK
        # ///////////////////////////////////////////////////////////////
        self.ui.btnSructure.clicked.connect(lambda: BtnGetPath.run(self.ui.editStructure, 'analysis_gro'))
        self.ui.btnTrajectory.clicked.connect(lambda: BtnGetPath.run(self.ui.editTrajectory, 'analysis_xtc'))
        self.ui.btnResult.clicked.connect(lambda: BtnGetPath.run(self.ui.editResult, 'analysis_result'))

        self.ui.btnNext.clicked.connect(openCloseRightBox)
        self.ui.btnNext.clicked.connect(lambda: NextClick(self.ui))

        # //////////////////////////Figure/////////////////////////////////////
        # btnFigure Click
        # ///////////////////////////////////////////////////////////////
        self.ui.FigureColorLayout = None
        self.ui.FigureShapeWidget = None
        self.ui.btn_figure_run.clicked.connect(lambda: FigurePage.figureBtnMakeFigure(self.ui))
        self.ui.figure_btn_path.clicked.connect(lambda: BtnGetPath.run(self.ui.figure_edit_path, 'figure_xlsx'))
        # self.ui.btn_gene_run.setStyleSheet("background-color: red; color: white;")

        # SHOW APP
        # ///////////////////////////////////////////////////////////////
        self.show()

        # SET CUSTOM THEME
        # ///////////////////////////////////////////////////////////////
        # useCustomTheme = True
        # themeFile = "themes/py_dracula_light.qss"

        # # SET THEME AND HACKS
        # if useCustomTheme:
        #     # LOAD AND APPLY STYLE
        #     UIFunctions.theme(self, themeFile, True)

            # SET HACKS
            # AppFunctions.setThemeHack(self)

        # SET HOME PAGE AND SELECT MENU
        widgets.stackedWidget.setCurrentWidget(widgets.page_home)
        widgets.btn_home.setStyleSheet(UIFunctions.selectMenu(widgets.btn_home.styleSheet()))

    # VMD
    def dragEnterEvent(self, event: QDragEnterEvent):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event: QDropEvent):
        urls = event.mimeData().urls()
        if not urls:
            return
        file_path = urls[0].toLocalFile()
        if file_path.lower().endswith('.csv'):
            self.loadCSV(file_path)
        else:
            self.ui.vmd_label.setText("Please drop a CSV file")

    def loadCSV(self, csv_path):
        if not os.path.exists(csv_path):
            self.ui.vmd_label.setText("CSV file not found!")
            return
        try:
            self.valid_comments, self.data = read_excel(csv_path)
            if self.data is not None:
                self.data.rename(columns={'Resid': 'resid', 'Resname': 'resname'}, inplace=True)
                # self.data = self.data.drop(columns=['resname', 'coordinations'], errors='ignore')
                frame_cols = [col for col in self.data.columns if col != 'resid']
                self.displayData(frame_cols)
                self.ui.vmd_label.setText(f"CSV loaded successfully. Valid comment: {self.valid_comments}")
            else:
                self.ui.vmd_label.setText("Failed to load CSV data")
        except Exception as e:
            self.ui.vmd_label.setText(f"Error loading CSV: {e}")

    def displayData(self, frame_cols):
        self.ui.vmd_tablewidget.clear()
        self.ui.vmd_tablewidget.setRowCount(len(self.data))
        self.ui.vmd_tablewidget.setColumnCount(len(frame_cols) + 1)
        self.ui.vmd_tablewidget.setHorizontalHeaderLabels(['resid'] + frame_cols)

        for i, row in self.data.iterrows():
            self.ui.vmd_tablewidget.setItem(i, 0, QTableWidgetItem(str(row['resid'])))
            for j, frame in enumerate(frame_cols):
                self.ui.vmd_tablewidget.setItem(i, j + 1, QTableWidgetItem(str(row[frame])))

    def pushStartVMD(self):
        try:
            self.vmd = VMDTcp(self.rctl_path, self.vmd_path)
            self.ui.vmd_btn_start.setEnabled(False)
            self.ui.vmd_btn_stop.setEnabled(True)
            response = self.vmd.start()
            if response == -1:
                self.ui.vmd_label.setText("Could not connect to VMD!")
                self.ui.vmd_btn_start.setEnabled(True)
                self.ui.vmd_btn_stop.setEnabled(False)
            else:
                self.ui.vmd_label.setText("VMD started and connected")
                self.connected = True
        except FileNotFoundError as e:
            self.ui.vmd_label.setText(str(e))

    def pushStopVMD(self):
        if self.vmd:
            self.vmd.stop()
        self.connected = False
        self.ui.vmd_btn_start.setEnabled(True)
        self.ui.vmd_btn_stop.setEnabled(False)
        self.ui.vmd_label.setText("VMD stopped")

    def onSelectionChanged(self):
        if not self.connected:
            self.ui.vmd_label.setText("VMD is not connected")
            return

        selected_items = self.ui.vmd_tablewidget.selectedItems()
        if not selected_items:
            return

        columns = set(item.column() for item in selected_items)
        if len(columns) != 1:
            self.ui.vmd_label.setText("Please select cells in the same frame (column)")
            return

        column = columns.pop()
        if column == 0:
            self.ui.vmd_label.setText("Please select cells in a frame column (not resid)")
            return

        frame = self.ui.vmd_tablewidget.horizontalHeaderItem(column).text()

        resids = []
        for item in selected_items:
            row = item.row()
            resid = self.ui.vmd_tablewidget.item(row, 0).text()
            if resid not in resids:
                resids.append(resid)

        self.vmd.send_command(VMDCommands.gotoFrame(frame))
        self.vmd.send_command(VMDCommands.highlightResid(resids))
        self.ui.vmd_label.setText(f"Showing resids {', '.join(resids)} at frame {frame}")

    # INITIAL SETTINGS
    # Post here your directions for your main UI
    def initialSettings(self):

        # Info
        self.ui.FigureInfo = None
        self.ui.GenerationInfo = None
        self.ui.AnalysisInfo = None
        # UI_GLOBAL
        self.ui.figure_extra = 1
        # ////////////////////////Generation////////////////////////////////////////
        # # box
        # x
        self.ui.spin_box_x.setMinimum(0)
        self.ui.spin_box_x.setMaximum(1e9)
        # y
        self.ui.spin_box_y.setMinimum(0)
        self.ui.spin_box_y.setMaximum(1e9)
        # z
        self.ui.spin_box_z.setMinimum(0)
        self.ui.spin_box_z.setMaximum(1e9)

        # # lnb
        # gas
        self.ui.spin_gas_density.setMinimum(0)
        self.ui.spin_gas_density.setMaximum(1e9)
        # area
        self.ui.spin_area_5.setMinimum(0)
        self.ui.spin_area_5.setMaximum(1e9)

        # # solvent
        # salt
        self.ui.spin_salt.setMinimum(0)
        self.ui.spin_salt.setMaximum(1e9)

        # # path
        self.ui.edit_gene_path.setReadOnly(True)

        # //////////////////////////Analysis/////////////////////////////////////
        # # restriction of lineEdit
        #  first frame
        self.ui.editFirstFrame.setMinimum(0)
        self.ui.editFirstFrame.setMaximum(1e9)

        # last frame
        self.ui.editLastFrame.setMinimum(-1)
        self.ui.editLastFrame.setMaximum(1e9)
        self.ui.editLastFrame.setValue(-1)

        # step
        self.ui.editStep.setMinimum(1)
        self.ui.editStep.setMaximum(100000)

        # k
        self.ui.editK.setMinimum(3)
        self.ui.editK.setMaximum(1e9)
        self.ui.editK.setValue(21)

        # path
        self.ui.editStructure.setReadOnly(True)  # Path of gro
        self.ui.editTrajectory.setReadOnly(True)  # Path of xtx
        self.ui.editResult.setReadOnly(True)  # Path of Save



    # BUTTONS CLICK
    # Post here your functions for clicked buttons
    # ///////////////////////////////////////////////////////////////
    def buttonLeftClick(self):
        # GET BUTTON CLICKED
        btn = self.sender()
        btnName = btn.objectName()

        # SHOW HOME PAGE
        if btnName == "btn_home":
            widgets.stackedWidget.setCurrentWidget(widgets.page_home)
            # 用来改变侧边按钮点击后的显示
            UIFunctions.resetStyle(self, btnName)
            btn.setStyleSheet(UIFunctions.selectMenu(btn.styleSheet()))

        # SHOW WIDGETS PAGE
        if btnName == "btn_generate":
            widgets.stackedWidget.setCurrentWidget(widgets.page_generation)
            UIFunctions.resetStyle(self, btnName)
            btn.setStyleSheet(UIFunctions.selectMenu(btn.styleSheet()))

        # SHOW NEW PAGE
        if btnName == "btn_figure":
            widgets.stackedWidget.setCurrentWidget(widgets.page_figure)  # SET PAGE
            UIFunctions.resetStyle(self, btnName)  # RESET ANOTHERS BUTTONS SELECTED
            btn.setStyleSheet(UIFunctions.selectMenu(btn.styleSheet()))  # SELECT MENU

        if btnName == "btn_analysis":
            widgets.stackedWidget.setCurrentWidget(widgets.page_analysis)
            UIFunctions.resetStyle(self, btnName)
            btn.setStyleSheet(UIFunctions.selectMenu(btn.styleSheet()))

        if btnName == 'btn_data_process':
            widgets.stackedWidget.setCurrentWidget(widgets.page_vmd)
            UIFunctions.resetStyle(self, btnName)
            btn.setStyleSheet(UIFunctions.selectMenu(btn.styleSheet()))

    # RESIZE EVENTS
    # ///////////////////////////////////////////////////////////////
    def resizeEvent(self, event):
        # Update Size Grips
        UIFunctions.resize_grips(self)

    # MOUSE CLICK EVENTS
    # ///////////////////////////////////////////////////////////////
    def mousePressEvent(self, event):
        # SET DRAG POS WINDOW
        self.dragPos = event.globalPos()
        #
        # PRINT MOUSE EVENTS
        # if event.buttons() == Qt.LeftButton:
        #     print('Mouse click: LEFT CLICK')
        # if event.buttons() == Qt.RightButton:
        #     print('Mouse click: RIGHT CLICK')

    def setup_file_detection(self):
        """设置文件检测功能"""
        # 创建动态tab管理器，传入UI实例
        self.tab_manager = DynamicTabManager(self.ui)
        
        # 将tab_manager设置到ui上，以便参数读取器可以访问
        self.ui.tab_manager = self.tab_manager
        
        # 设置文件路径输入框为可编辑
        self.ui.figure_edit_path.setReadOnly(True)
        
        # 添加文件选择功能到路径输入框
        self.ui.figure_edit_path.textChanged.connect(self.on_file_path_changed)
    
    def setup_theme_switch(self):
        """设置主题切换功能"""
        try:
            # 将主题切换按钮添加到顶部菜单
            self.theme_manager = add_theme_button_to_top_menu(self.ui, self)
            print("主题切换功能已成功添加到顶部菜单")
        except Exception as e:
            print(f"添加主题切换功能失败: {e}")
            # 如果添加失败，至少创建主题管理器
            self.theme_manager = ThemeManager(self)
    
    def ensure_figure_widget_display(self):
        """确保figure页面默认显示信息文本框"""
        try:
            # 强制清理figure页面并重新设置
            self._force_cleanup_figure_page()
            
            # tabWidget_lipids已被删除，无需隐藏
            
            # 添加信息文本框到文件导入区域和RUN按钮之间
            self._add_info_textbox_between_import_and_run()
            
            # 确保widget_2可见
            if hasattr(self.ui, 'widget_2'):
                self.ui.widget_2.setVisible(True)
                self.ui.widget_2.show()
                self.ui.widget_2.raise_()  # 确保在顶层
                print("已确保widget_2可见")
                
        except Exception as e:
            print(f"确保figure widget显示失败: {e}")
    
    def _add_info_textbox_between_import_and_run(self):
        """在文件导入区域和RUN按钮之间添加信息文本框"""
        try:
            from PySide6.QtWidgets import QTextEdit
            
            # 创建信息文本框（合并标题和内容）
            info_textbox = QTextEdit()
            info_textbox.setObjectName("figure_info_textbox")
            info_textbox.setMinimumHeight(300)
            info_textbox.setMaximumHeight(350)
            info_textbox.setReadOnly(True)
            info_textbox.setStyleSheet("""
                QTextEdit {
                    font: 12pt "华文细黑";
                    background-color: rgba(33, 37, 43, 0.9);
                    border: 2px solid rgb(189, 147, 249);
                    border-radius: 8px;
                    padding: 15px;
                    color: white;
                    margin: 5px;
                }
                QTextEdit:focus {
                    border: 2px solid rgb(255, 121, 198);
                }
            """)
            
            # 设置默认内容（包含标题）
            default_content = """📊 Supported Plot Types

📈 LIPIDS Analysis:
• Line Chart: Time series analysis of lipid properties
• Bar Chart: Statistical comparison of lipid groups  
• Scatter Plot: Correlation analysis between lipid parameters

🫧 BUBBLE Analysis:
• Line Chart: Bubble size evolution over time
• Bar Chart: Bubble distribution statistics

💡 Usage: Import CSV file with TYPE:Lipids or TYPE:Bubble in line 4"""
            
            info_textbox.setText(default_content)
            
            # 将信息组件添加到widget_2的布局中，在frame_8和btn_figure_run之间
            if hasattr(self.ui, 'widget_2'):
                layout = self.ui.widget_2.layout()
                if layout:
                    # 找到btn_figure_run的位置
                    run_button_index = -1
                    for i in range(layout.count()):
                        item = layout.itemAt(i)
                        if item and item.widget() == self.ui.btn_figure_run:
                            run_button_index = i
                            break
                    
                    if run_button_index >= 0:
                        # 在RUN按钮之前插入信息文本框
                        layout.insertWidget(run_button_index, info_textbox)
                        print("已在文件导入区域和RUN按钮之间添加信息文本框")
                    else:
                        # 如果找不到RUN按钮，添加到末尾
                        layout.addWidget(info_textbox)
                        print("已添加信息文本框到布局末尾")
                    
        except Exception as e:
            print(f"添加信息文本框失败: {e}")
    
    
    
    def _force_cleanup_figure_page(self):
        """强制清理figure页面，移除任何可能的table widget"""
        try:
            from PySide6.QtWidgets import QTableWidget
            
            # 遍历figure页面的所有子widget
            all_children = self.ui.page_figure.findChildren(QTableWidget)
            for child in all_children:
                # 删除任何QTableWidget（除了我们不想删除的）
                if child.objectName() != 'info_textbox':  # 保护信息文本框不被删除
                    child.setParent(None)  # 从父widget中移除
                    child.deleteLater()    # 标记为删除
                    print(f"已移除figure页面中的table widget: {child.objectName()}")
            
            # 确保页面布局正确
            layout = self.ui.page_figure.layout()
            if layout:
                # 清理布局中可能存在的无效项目
                for i in reversed(range(layout.count())):
                    item = layout.itemAt(i)
                    if item and item.widget():
                        widget = item.widget()
                        # 如果是QTableWidget且不是我们需要的，从布局中移除
                        if isinstance(widget, QTableWidget):
                            layout.removeWidget(widget)
                            print(f"已从布局中移除table widget: {widget.objectName()}")
                
                # 确保核心widget在布局中
                if hasattr(self.ui, 'widget_2'):
                    widget_in_layout = False
                    for i in range(layout.count()):
                        item = layout.itemAt(i)
                        if item and item.widget() == self.ui.widget_2:
                            widget_in_layout = True
                            break
                    
                    if not widget_in_layout:
                        layout.insertWidget(0, self.ui.widget_2)  # 插入到第一个位置
                        print("已重新添加widget_2到figure页面布局")
                        
        except Exception as e:
            print(f"强制清理figure页面失败: {e}")

    def on_file_path_changed(self):
        """当文件路径改变时的回调"""
        file_path = self.ui.figure_edit_path.text().strip()
        if file_path and os.path.exists(file_path):
            # 自动检测文件类型
            self.auto_detect_file_type(file_path)

    def auto_detect_file_type(self, file_path):
        """自动检测文件类型并创建相应的TabWidget"""
        try:
            # 检测文件类型
            file_type = FileTypeDetector.detect_file_type(file_path)
            print(f"检测到文件类型: {file_type}")
            
            # 移除信息文本框
            self._remove_info_textbox()
            
            # 确保FigureInfo被正确设置（重要！）
            self._ensure_figure_info_for_file(file_path)
            
            # 创建并替换TabWidget
            tab_widget = self.tab_manager.replace_tab_widget(file_type, file_path)
            
            # 显示成功消息
            from PySide6.QtWidgets import QMessageBox
            tab_name = self.tab_manager._get_tab_name(file_type)
            QMessageBox.information(self, "检测成功", 
                                  f"文件类型检测完成！\n文件: {os.path.basename(file_path)}\n类型: {tab_name}")
            
        except Exception as e:
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.critical(self, "错误", f"检测文件类型时发生错误: {str(e)}")
    
    def _ensure_figure_info_for_file(self, file_path):
        """确保FigureInfo被正确设置"""
        try:
            # 导入create_parameter_reader函数
            from modules.Fuctions_Figure import create_parameter_reader
            
            # 检查是否需要重新创建FigureInfo
            if (self.ui.FigureInfo is None or 
                not hasattr(self.ui.FigureInfo, 'path_figure') or 
                self.ui.FigureInfo.path_figure != file_path):
                
                print(f"正在为文件创建参数读取器: {file_path}")
                
                # 创建参数读取器
                self.ui.FigureInfo = create_parameter_reader(self.ui)
                
                print("参数读取器创建成功")
            else:
                print("参数读取器已存在且文件路径未改变")
                
        except Exception as e:
            print(f"创建参数读取器失败: {e}")
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "警告", f"文件信息读取失败: {str(e)}")
    
    def _remove_info_textbox(self):
        """移除信息文本框"""
        try:
            if hasattr(self.ui, 'widget_2'):
                layout = self.ui.widget_2.layout()
                if layout:
                    # 查找并移除信息文本框
                    for i in range(layout.count()):
                        item = layout.itemAt(i)
                        if item and item.widget():
                            widget = item.widget()
                            if widget.objectName() == "figure_info_textbox":
                                layout.removeWidget(widget)
                                widget.setParent(None)
                                widget.deleteLater()
                                print("已移除信息文本框")
                                break
                        
        except Exception as e:
            print(f"移除信息文本框失败: {e}")


if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon("icon.ico"))
    window = MainWindow()
    sys.exit(app.exec())

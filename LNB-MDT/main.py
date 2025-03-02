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

# IMPORT / GUI AND MODULES AND WIDGETS
# ///////////////////////////////////////////////////////////////

from modules import *
from widgets import *
from figure import *
from analysis import *

os.environ["QT_FONT_DPI"] = "96"  # FIX Problem for High DPI and Scale above 100%

# SET AS GLOBAL WIDGETS
# ///////////////////////////////////////////////////////////////
widgets = None


class MainWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)

        # SET AS GLOBAL WIDGETS
        # ///////////////////////////////////////////////////////////////
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        global widgets
        widgets = self.ui

        # USE CUSTOM TITLE BAR | USE AS "False" FOR MAC OR LINUX
        # ///////////////////////////////////////////////////////////////
        Settings.ENABLE_CUSTOM_TITLE_BAR = True

        # TOGGLE MENU
        # ///////////////////////////////////////////////////////////////
        widgets.toggleButton.clicked.connect(lambda: UIFunctions.toggleMenu(self, True))

        # SET UI DEFINITIONS
        # ///////////////////////////////////////////////////////////////
        UIFunctions.uiDefinitions(self)

        for child in widgets.extraRightBox.children():
            child.deleteLater()

        # LEFT MENUS
        widgets.btn_home.clicked.connect(self.buttonLeftClick)
        widgets.btn_generate.clicked.connect(self.buttonLeftClick)
        widgets.btn_figure.clicked.connect(self.buttonLeftClick)
        widgets.btn_analysis.clicked.connect(self.buttonLeftClick)
        widgets.btn_data_process.clicked.connect(self.buttonLeftClick)
        self.initialSettings()

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

        self.ui.figure_line_btn_color_2.clicked.connect(lambda: FigurePage.figureBtnColor(self.ui))
        self.ui.figure_line_btn_color_2.clicked.connect(openCloseFigureColorBox)
        # self.ui.figure_bar_btn_trend_2.clicked.connect(lambda: FigurePage.figure_trend_line_color_btn)
        self.ui.figure_bar_btn_color_2.clicked.connect(lambda: FigurePage.figureBtnColor(self.ui))
        self.ui.figure_bar_btn_color_2.clicked.connect(openCloseFigureColorBox)
        self.ui.figure_scatter_btn_shape_2.clicked.connect(openCloseFigureShapeBox)
        self.ui.figure_scatter_btn_shape_2.clicked.connect(lambda: FigurePage.figureBtnShape(self.ui))
        self.ui.tabWidget.currentChanged.connect(openCloseFigureExtra)
        # //////////////////////////DataProcess/////////////////////////////////////
        self.ui.dataprocess_btn_path.clicked.connect(lambda: BtnGetPath.run(self.ui.dataprocess_edit_path
                                                                            , 'multiple_files'
                                                                            , self.ui.dataprocess_edit_files))
        # self.ui.
        self.ui.dataprocess_btn_save.clicked.connect(lambda: BtnGetPath.run(self.ui.dataprocess_edit_save
                                                                            , 'data_save'
                                                                            ))
        self.ui.btn_data_process_run.clicked.connect(lambda: merge_files_lipid(self.ui.dataprocess_edit_files.toPlainText()
                                                                               , self.ui.dataprocess_edit_save.text()
                                                                               , ui_info=[self.ui.data_process_franes_edit.text(),
                                                                                          self.ui.data_process_spinbox_start.value(),
                                                                                          self.ui.data_process_spinbox_end.value(),
                                                                                          self.ui.data_process_spinbox_step.value()]))
        # SHOW APP
        # ///////////////////////////////////////////////////////////////
        self.show()

        # SET CUSTOM THEME
        # ///////////////////////////////////////////////////////////////
        useCustomTheme = False
        themeFile = "themes\py_dracula_light.qss"

        # SET THEME AND HACKS
        if useCustomTheme:
            # LOAD AND APPLY STYLE
            UIFunctions.theme(self, themeFile, True)

            # SET HACKS
            AppFunctions.setThemeHack(self)

        # SET HOME PAGE AND SELECT MENU
        # ///////////////////////////////////////////////////////////////
        widgets.stackedWidget.setCurrentWidget(widgets.page_home)
        widgets.btn_home.setStyleSheet(UIFunctions.selectMenu(widgets.btn_home.styleSheet()))

        # Top of the main UI
        self.ui.contentTopBg.setStyleSheet('background-color:#343B48;')
        # stackWidget
        self.ui.pagesContainer.setStyleSheet("background-color: #2C313C;")
        # Bottom of the main UI
        self.ui.bottomBar.setStyleSheet('background-color:#343B48;')

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
        # //////////////////////////DataProcess/////////////////////////////////////
        self.ui.data_process_spinbox_start.setValue(0)
        self.ui.data_process_spinbox_start.setMaximum(100000000)
        self.ui.data_process_spinbox_end.setValue(0)
        self.ui.data_process_spinbox_end.setMaximum(100000000)
        self.ui.data_process_spinbox_step.setMinimum(1)


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
            widgets.stackedWidget.setCurrentWidget(widgets.page_dataprocess)
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


if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    app.setWindowIcon(QIcon("icon.ico"))
    window = MainWindow()
    sys.exit(app.exec())

# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'main.ui'
##
## Created by: Qt User Interface Compiler version 6.8.1
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor,
    QFont, QFontDatabase, QGradient, QIcon,
    QImage, QKeySequence, QLinearGradient, QPainter,
    QPalette, QPixmap, QRadialGradient, QTransform)
from PySide6.QtWidgets import (QApplication, QButtonGroup, QComboBox, QDoubleSpinBox,
    QFrame, QGridLayout, QHBoxLayout, QHeaderView,
    QLabel, QLineEdit, QMainWindow, QPushButton,
    QRadioButton, QSizePolicy, QSpinBox, QStackedWidget,
    QTabWidget, QTableWidget, QTableWidgetItem, QVBoxLayout,
    QWidget)
from .resources_rc import *

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(800, 638)
        MainWindow.setMinimumSize(QSize(0, 0))
        MainWindow.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")
        self.styleSheet = QWidget(MainWindow)
        self.styleSheet.setObjectName(u"styleSheet")
        self.styleSheet.setEnabled(True)
        font = QFont()
        font.setFamilies([u"Segoe UI"])
        font.setPointSize(10)
        font.setBold(False)
        font.setItalic(False)
        self.styleSheet.setFont(font)
        self.styleSheet.setStyleSheet(u"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"\n"
"SET APP STYLESHEET - FULL STYLES HERE\n"
"DARK THEME - DRACULA COLOR BASED\n"
"\n"
"///////////////////////////////////////////////////////////////////////////////////////////////// */\n"
"\n"
"QWidget{\n"
"   color: rgb(221, 221, 221);\n"
"   font: 10pt \"Segoe UI\";\n"
"}\n"
"\n"
"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"Tooltip */\n"
"QToolTip {\n"
"   color: #ffffff;\n"
"   background-color: rgba(33, 37, 43, 180);\n"
"   border: 1px solid rgb(44, 49, 58);\n"
"   background-image: none;\n"
"   background-position: left center;\n"
"    background-repeat: no-repeat;\n"
"   border: none;\n"
"   border-left: 2px solid rgb(255, 121, 198);\n"
"   text-align: left;\n"
"   padding-left: 8px;\n"
"   margin: 0px;\n"
"}\n"
"\n"
"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"Bg App */\n"
"#bgApp {   \n"
"   background"
                        "-color: rgb(40, 44, 52);\n"
"   border: 1px solid rgb(44, 49, 58);\n"
"}\n"
"\n"
"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"Left Menu */\n"
"#leftMenuBg {  \n"
"   background-color: rgb(33, 37, 43);\n"
"}\n"
"#topLogo {\n"
"   background-color: rgb(33, 37, 43);\n"
"   background-image: url(:/images/images/images/PyDracula.png);\n"
"   background-position: centered;\n"
"   background-repeat: no-repeat;\n"
"}\n"
"#titleLeftApp { font: 63 12pt \"Segoe UI Semibold\"; }\n"
"#titleLeftDescription { font: 8pt \"Segoe UI\"; color: rgb(189, 147, 249); }\n"
"\n"
"/* MENUS */\n"
"#topMenu .QPushButton {    \n"
"   background-position: left center;\n"
"    background-repeat: no-repeat;\n"
"   border: none;\n"
"   border-left: 22px solid transparent;\n"
"   background-color: transparent;\n"
"   text-align: left;\n"
"   padding-left: 44px;\n"
"}\n"
"#topMenu .QPushButton:hover {\n"
"   background-color: rgb(40, 44, 52);\n"
"}\n"
"#topMenu .QPushButton:pressed {    \n"
"   background-color: rgb(18"
                        "9, 147, 249);\n"
"   color: rgb(255, 255, 255);\n"
"}\n"
"#bottomMenu .QPushButton { \n"
"   background-position: left center;\n"
"    background-repeat: no-repeat;\n"
"   border: none;\n"
"   border-left: 20px solid transparent;\n"
"   background-color:transparent;\n"
"   text-align: left;\n"
"   padding-left: 44px;\n"
"}\n"
"#bottomMenu .QPushButton:hover {\n"
"   background-color: rgb(40, 44, 52);\n"
"}\n"
"#bottomMenu .QPushButton:pressed { \n"
"   background-color: rgb(189, 147, 249);\n"
"   color: rgb(255, 255, 255);\n"
"}\n"
"#leftMenuFrame{\n"
"   border-top: 3px solid rgb(44, 49, 58);\n"
"}\n"
"\n"
"/* Toggle Button */\n"
"#toggleButton {\n"
"   background-position: left center;\n"
"    background-repeat: no-repeat;\n"
"   border: none;\n"
"   border-left: 20px solid transparent;\n"
"   background-color: rgb(37, 41, 48);\n"
"   text-align: left;\n"
"   padding-left: 44px;\n"
"   color: rgb(113, 126, 149);\n"
"}\n"
"#toggleButton:hover {\n"
"   background-color: rgb(40, 44, 52);\n"
"}\n"
"#toggleButton:pressed {\n"
"   background-color: rgb("
                        "189, 147, 249);\n"
"}\n"
"\n"
"/* Title Menu */\n"
"#titleRightInfo { padding-left: 10px; }\n"
"\n"
"\n"
"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"Extra Tab */\n"
"#extraLeftBox {    \n"
"   background-color: rgb(44, 49, 58);\n"
"}\n"
"#extraTopBg{   \n"
"   background-color: rgb(189, 147, 249)\n"
"}\n"
"\n"
"/* Icon */\n"
"#extraIcon {\n"
"   background-position: center;\n"
"   background-repeat: no-repeat;\n"
"   background-image: url(:/icons/images/icons/icon_settings.png);\n"
"}\n"
"\n"
"/* Label */\n"
"#extraLabel { color: rgb(255, 255, 255); }\n"
"\n"
"/* Btn Close */\n"
"#extraCloseColumnBtn { background-color: rgba(255, 255, 255, 0); border: none;  border-radius: 5px; }\n"
"#extraCloseColumnBtn:hover { background-color: rgb(196, 161, 249); border-style: solid; border-radius: 4px; }\n"
"#extraCloseColumnBtn:pressed { background-color: rgb(180, 141, 238); border-style: solid; border-radius: 4px; }\n"
"\n"
"/* Extra Content */\n"
"#extraContent{\n"
"   border"
                        "-top: 3px solid rgb(40, 44, 52);\n"
"}\n"
"\n"
"/* Extra Top Menus */\n"
"#extraTopMenu .QPushButton {\n"
"background-position: left center;\n"
"    background-repeat: no-repeat;\n"
"   border: none;\n"
"   border-left: 22px solid transparent;\n"
"   background-color:transparent;\n"
"   text-align: left;\n"
"   padding-left: 44px;\n"
"}\n"
"#extraTopMenu .QPushButton:hover {\n"
"   background-color: rgb(40, 44, 52);\n"
"}\n"
"#extraTopMenu .QPushButton:pressed {   \n"
"   background-color: rgb(189, 147, 249);\n"
"   color: rgb(255, 255, 255);\n"
"}\n"
"\n"
"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"Content App */\n"
"#contentTopBg{ \n"
"   background-color: rgb(33, 37, 43);\n"
"}\n"
"#contentBottom{\n"
"   border-top: 3px solid rgb(44, 49, 58);\n"
"}\n"
"\n"
"/* Top Buttons */\n"
"#rightButtons .QPushButton { background-color: rgba(255, 255, 255, 0); border: none;  border-radius: 5px; }\n"
"#rightButtons .QPushButton:hover { background-color: rgb(44, 49, 57); border-sty"
                        "le: solid; border-radius: 4px; }\n"
"#rightButtons .QPushButton:pressed { background-color: rgb(23, 26, 30); border-style: solid; border-radius: 4px; }\n"
"\n"
"/* Theme Settings */\n"
"#extraRightBox { background-color: rgb(44, 49, 58); }\n"
"#themeSettingsTopDetail { background-color: rgb(189, 147, 249); }\n"
"\n"
"/* Bottom Bar */\n"
"#bottomBar { background-color: rgb(44, 49, 58); }\n"
"#bottomBar QLabel { font-size: 11px; color: rgb(113, 126, 149); padding-left: 10px; padding-right: 10px; padding-bottom: 2px; }\n"
"\n"
"/* CONTENT SETTINGS */\n"
"/* MENUS */\n"
"#contentSettings .QPushButton {    \n"
"   background-position: left center;\n"
"    background-repeat: no-repeat;\n"
"   border: none;\n"
"   border-left: 22px solid transparent;\n"
"   background-color:transparent;\n"
"   text-align: left;\n"
"   padding-left: 44px;\n"
"}\n"
"#contentSettings .QPushButton:hover {\n"
"   background-color: rgb(40, 44, 52);\n"
"}\n"
"#contentSettings .QPushButton:pressed {    \n"
"   background-color: rgb(189, 147, 249);\n"
"   color: rgb"
                        "(255, 255, 255);\n"
"}\n"
"\n"
"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"QTableWidget */\n"
"QTableWidget { \n"
"   background-color: transparent;\n"
"   padding: 10px;\n"
"   border-radius: 5px;\n"
"   gridline-color: rgb(44, 49, 58);\n"
"   border-bottom: 1px solid rgb(44, 49, 60);\n"
"}\n"
"QTableWidget::item{\n"
"   border-color: rgb(44, 49, 60);\n"
"   padding-left: 5px;\n"
"   padding-right: 5px;\n"
"   gridline-color: rgb(44, 49, 60);\n"
"}\n"
"QTableWidget::item:selected{\n"
"   background-color: rgb(189, 147, 249);\n"
"}\n"
"QHeaderView::section{\n"
"   background-color: rgb(33, 37, 43);\n"
"   max-width: 30px;\n"
"   border: 1px solid rgb(44, 49, 58);\n"
"   border-style: none;\n"
"    border-bottom: 1px solid rgb(44, 49, 60);\n"
"    border-right: 1px solid rgb(44, 49, 60);\n"
"}\n"
"QTableWidget::horizontalHeader {   \n"
"   background-color: rgb(33, 37, 43);\n"
"}\n"
"QHeaderView::section:horizontal\n"
"{\n"
"    border: 1px solid rgb(33, 37, 43);\n"
"   background-co"
                        "lor: rgb(33, 37, 43);\n"
"   padding: 3px;\n"
"   border-top-left-radius: 7px;\n"
"    border-top-right-radius: 7px;\n"
"}\n"
"QHeaderView::section:vertical\n"
"{\n"
"    border: 1px solid rgb(44, 49, 60);\n"
"}\n"
"\n"
"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"LineEdit */\n"
"QLineEdit {\n"
"   background-color: rgb(33, 37, 43);\n"
"   border-radius: 5px;\n"
"   border: 2px solid rgb(33, 37, 43);\n"
"   padding-left: 10px;\n"
"   selection-color: rgb(255, 255, 255);\n"
"   selection-background-color: rgb(255, 121, 198);\n"
"}\n"
"QLineEdit:hover {\n"
"   border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QLineEdit:focus {\n"
"   border: 2px solid rgb(91, 101, 124);\n"
"}\n"
"\n"
"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"PlainTextEdit */\n"
"QPlainTextEdit {\n"
"   background-color: rgb(27, 29, 35);\n"
"   border-radius: 5px;\n"
"   padding: 10px;\n"
"   selection-color: rgb(255, 255, 255);\n"
"   selection-background-c"
                        "olor: rgb(255, 121, 198);\n"
"}\n"
"QPlainTextEdit  QScrollBar:vertical {\n"
"    width: 8px;\n"
" }\n"
"QPlainTextEdit  QScrollBar:horizontal {\n"
"    height: 8px;\n"
" }\n"
"QPlainTextEdit:hover {\n"
"   border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QPlainTextEdit:focus {\n"
"   border: 2px solid rgb(91, 101, 124);\n"
"}\n"
"\n"
"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"ScrollBars */\n"
"QScrollBar:horizontal {\n"
"    border: none;\n"
"    background: rgb(52, 59, 72);\n"
"    height: 8px;\n"
"    margin: 0px 21px 0 21px;\n"
"   border-radius: 0px;\n"
"}\n"
"QScrollBar::handle:horizontal {\n"
"    background: rgb(189, 147, 249);\n"
"    min-width: 25px;\n"
"   border-radius: 4px\n"
"}\n"
"QScrollBar::add-line:horizontal {\n"
"    border: none;\n"
"    background: rgb(55, 63, 77);\n"
"    width: 20px;\n"
"   border-top-right-radius: 4px;\n"
"    border-bottom-right-radius: 4px;\n"
"    subcontrol-position: right;\n"
"    subcontrol-origin: margin;\n"
"}\n"
""
                        "QScrollBar::sub-line:horizontal {\n"
"    border: none;\n"
"    background: rgb(55, 63, 77);\n"
"    width: 20px;\n"
"   border-top-left-radius: 4px;\n"
"    border-bottom-left-radius: 4px;\n"
"    subcontrol-position: left;\n"
"    subcontrol-origin: margin;\n"
"}\n"
"QScrollBar::up-arrow:horizontal, QScrollBar::down-arrow:horizontal\n"
"{\n"
"     background: none;\n"
"}\n"
"QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal\n"
"{\n"
"     background: none;\n"
"}\n"
" QScrollBar:vertical {\n"
"   border: none;\n"
"    background: rgb(52, 59, 72);\n"
"    width: 8px;\n"
"    margin: 21px 0 21px 0;\n"
"   border-radius: 0px;\n"
" }\n"
" QScrollBar::handle:vertical { \n"
"   background: rgb(189, 147, 249);\n"
"    min-height: 25px;\n"
"   border-radius: 4px\n"
" }\n"
" QScrollBar::add-line:vertical {\n"
"     border: none;\n"
"    background: rgb(55, 63, 77);\n"
"     height: 20px;\n"
"   border-bottom-left-radius: 4px;\n"
"    border-bottom-right-radius: 4px;\n"
"     subcontrol-position: bottom;\n"
"     su"
                        "bcontrol-origin: margin;\n"
" }\n"
" QScrollBar::sub-line:vertical {\n"
"   border: none;\n"
"    background: rgb(55, 63, 77);\n"
"     height: 20px;\n"
"   border-top-left-radius: 4px;\n"
"    border-top-right-radius: 4px;\n"
"     subcontrol-position: top;\n"
"     subcontrol-origin: margin;\n"
" }\n"
" QScrollBar::up-arrow:vertical, QScrollBar::down-arrow:vertical {\n"
"     background: none;\n"
" }\n"
"\n"
" QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical {\n"
"     background: none;\n"
" }\n"
"\n"
"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"CheckBox */\n"
"QCheckBox::indicator {\n"
"    border: 3px solid rgb(52, 59, 72);\n"
"   width: 15px;\n"
"   height: 15px;\n"
"   border-radius: 10px;\n"
"    background: rgb(44, 49, 60);\n"
"}\n"
"QCheckBox::indicator:hover {\n"
"    border: 3px solid rgb(58, 66, 81);\n"
"}\n"
"QCheckBox::indicator:checked {\n"
"    background: 3px solid rgb(52, 59, 72);\n"
"   border: 3px solid rgb(52, 59, 72); \n"
"   back"
                        "ground-image: url(:/icons/images/icons/cil-check-alt.png);\n"
"}\n"
"\n"
"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"RadioButton */\n"
"QRadioButton::indicator {\n"
"    border: 3px solid rgb(52, 59, 72);\n"
"   width: 15px;\n"
"   height: 15px;\n"
"   border-radius: 10px;\n"
"    background: rgb(44, 49, 60);\n"
"}\n"
"QRadioButton::indicator:hover {\n"
"    border: 3px solid rgb(58, 66, 81);\n"
"}\n"
"QRadioButton::indicator:checked {\n"
"    background: 3px solid rgb(94, 106, 130);\n"
"   border: 3px solid rgb(52, 59, 72); \n"
"}\n"
"\n"
"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"ComboBox */\n"
"QComboBox{\n"
"   background-color: rgb(27, 29, 35);\n"
"   border-radius: 5px;\n"
"   border: 2px solid rgb(33, 37, 43);\n"
"   padding: 5px;\n"
"   padding-left: 10px;\n"
"}\n"
"QComboBox:hover{\n"
"   border: 2px solid rgb(64, 71, 88);\n"
"}\n"
"QComboBox::drop-down {\n"
"   subcontrol-origin: padding;\n"
"   subco"
                        "ntrol-position: top right;\n"
"   width: 25px; \n"
"   border-left-width: 3px;\n"
"   border-left-color: rgba(39, 44, 54, 150);\n"
"   border-left-style: solid;\n"
"   border-top-right-radius: 3px;\n"
"   border-bottom-right-radius: 3px;   \n"
"   background-image: url(:/icons/images/icons/cil-arrow-bottom.png);\n"
"   background-position: center;\n"
"   background-repeat: no-reperat;\n"
" }\n"
"QComboBox QAbstractItemView {\n"
"   color: rgb(255, 121, 198); \n"
"   background-color: rgb(33, 37, 43);\n"
"   padding: 10px;\n"
"   selection-background-color: rgb(39, 44, 54);\n"
"}\n"
"\n"
"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"Sliders */\n"
"QSlider::groove:horizontal {\n"
"    border-radius: 5px;\n"
"    height: 10px;\n"
"   margin: 0px;\n"
"   background-color: rgb(52, 59, 72);\n"
"}\n"
"QSlider::groove:horizontal:hover {\n"
"   background-color: rgb(55, 62, 76);\n"
"}\n"
"QSlider::handle:horizontal {\n"
"    background-color: rgb(189, 147, 249);\n"
"    border: none;\n"
"    h"
                        "eight: 10px;\n"
"    width: 10px;\n"
"    margin: 0px;\n"
"   border-radius: 5px;\n"
"}\n"
"QSlider::handle:horizontal:hover {\n"
"    background-color: rgb(195, 155, 255);\n"
"}\n"
"QSlider::handle:horizontal:pressed {\n"
"    background-color: rgb(255, 121, 198);\n"
"}\n"
"\n"
"QSlider::groove:vertical {\n"
"    border-radius: 5px;\n"
"    width: 10px;\n"
"    margin: 0px;\n"
"   background-color: rgb(52, 59, 72);\n"
"}\n"
"QSlider::groove:vertical:hover {\n"
"   background-color: rgb(55, 62, 76);\n"
"}\n"
"QSlider::handle:vertical {\n"
"    background-color: rgb(189, 147, 249);\n"
"   border: none;\n"
"    height: 10px;\n"
"    width: 10px;\n"
"    margin: 0px;\n"
"   border-radius: 5px;\n"
"}\n"
"QSlider::handle:vertical:hover {\n"
"    background-color: rgb(195, 155, 255);\n"
"}\n"
"QSlider::handle:vertical:pressed {\n"
"    background-color: rgb(255, 121, 198);\n"
"}\n"
"\n"
"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"CommandLinkButton */\n"
"QCommandLi"
                        "nkButton { \n"
"   color: rgb(255, 121, 198);\n"
"   border-radius: 5px;\n"
"   padding: 5px;\n"
"   color: rgb(255, 170, 255);\n"
"}\n"
"QCommandLinkButton:hover { \n"
"   color: rgb(255, 170, 255);\n"
"   background-color: rgb(44, 49, 60);\n"
"}\n"
"QCommandLinkButton:pressed {   \n"
"   color: rgb(189, 147, 249);\n"
"   background-color: rgb(52, 58, 71);\n"
"}\n"
"\n"
"/* /////////////////////////////////////////////////////////////////////////////////////////////////\n"
"Button */\n"
"#pagesContainer QPushButton {\n"
"   border: 2px solid rgb(52, 59, 72);\n"
"   border-radius: 5px;    \n"
"   background-color: rgb(52, 59, 72);\n"
"}\n"
"#pagesContainer QPushButton:hover {\n"
"   background-color: rgb(57, 65, 80);\n"
"   border: 2px solid rgb(61, 70, 86);\n"
"}\n"
"#pagesContainer QPushButton:pressed {  \n"
"   background-color: rgb(35, 40, 49);\n"
"   border: 2px solid rgb(43, 50, 61);\n"
"}\n"
"\n"
"")
        self.verticalLayout_23 = QVBoxLayout(self.styleSheet)
        self.verticalLayout_23.setSpacing(0)
        self.verticalLayout_23.setObjectName(u"verticalLayout_23")
        self.verticalLayout_23.setContentsMargins(10, 10, 10, 10)
        self.bgApp = QFrame(self.styleSheet)
        self.bgApp.setObjectName(u"bgApp")
        self.bgApp.setStyleSheet(u"")
        self.bgApp.setFrameShape(QFrame.NoFrame)
        self.bgApp.setFrameShadow(QFrame.Raised)
        self.appLayout = QHBoxLayout(self.bgApp)
        self.appLayout.setSpacing(0)
        self.appLayout.setObjectName(u"appLayout")
        self.appLayout.setContentsMargins(0, 0, 0, 0)
        self.leftMenuBg = QFrame(self.bgApp)
        self.leftMenuBg.setObjectName(u"leftMenuBg")
        self.leftMenuBg.setMinimumSize(QSize(60, 0))
        self.leftMenuBg.setMaximumSize(QSize(60, 16777215))
        self.leftMenuBg.setFrameShape(QFrame.NoFrame)
        self.leftMenuBg.setFrameShadow(QFrame.Raised)
        self.verticalLayout_3 = QVBoxLayout(self.leftMenuBg)
        self.verticalLayout_3.setSpacing(0)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.leftMenuFrame = QFrame(self.leftMenuBg)
        self.leftMenuFrame.setObjectName(u"leftMenuFrame")
        self.leftMenuFrame.setFrameShape(QFrame.NoFrame)
        self.leftMenuFrame.setFrameShadow(QFrame.Raised)
        self.verticalMenuLayout = QVBoxLayout(self.leftMenuFrame)
        self.verticalMenuLayout.setSpacing(0)
        self.verticalMenuLayout.setObjectName(u"verticalMenuLayout")
        self.verticalMenuLayout.setContentsMargins(0, 0, 0, 0)
        self.toggleBox = QFrame(self.leftMenuFrame)
        self.toggleBox.setObjectName(u"toggleBox")
        self.toggleBox.setMaximumSize(QSize(16777215, 45))
        self.toggleBox.setFrameShape(QFrame.NoFrame)
        self.toggleBox.setFrameShadow(QFrame.Raised)
        self.verticalLayout_4 = QVBoxLayout(self.toggleBox)
        self.verticalLayout_4.setSpacing(0)
        self.verticalLayout_4.setObjectName(u"verticalLayout_4")
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.toggleButton = QPushButton(self.toggleBox)
        self.toggleButton.setObjectName(u"toggleButton")
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.toggleButton.sizePolicy().hasHeightForWidth())
        self.toggleButton.setSizePolicy(sizePolicy)
        self.toggleButton.setMinimumSize(QSize(0, 45))
        font1 = QFont()
        font1.setFamilies([u"\u534e\u6587\u7ec6\u9ed1"])
        font1.setPointSize(14)
        font1.setBold(False)
        font1.setItalic(False)
        self.toggleButton.setFont(font1)
        self.toggleButton.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.toggleButton.setLayoutDirection(Qt.LeftToRight)
        self.toggleButton.setStyleSheet(u"background-image: url(:/icons/images/icons/icon_menu.png);\n"
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"")

        self.verticalLayout_4.addWidget(self.toggleButton)


        self.verticalMenuLayout.addWidget(self.toggleBox)

        self.topMenu = QFrame(self.leftMenuFrame)
        self.topMenu.setObjectName(u"topMenu")
        self.topMenu.setFrameShape(QFrame.NoFrame)
        self.topMenu.setFrameShadow(QFrame.Raised)
        self.verticalLayout_8 = QVBoxLayout(self.topMenu)
        self.verticalLayout_8.setSpacing(0)
        self.verticalLayout_8.setObjectName(u"verticalLayout_8")
        self.verticalLayout_8.setContentsMargins(0, 0, 0, 0)
        self.btn_home = QPushButton(self.topMenu)
        self.buttonGroup = QButtonGroup(MainWindow)
        self.buttonGroup.setObjectName(u"buttonGroup")
        self.buttonGroup.addButton(self.btn_home)
        self.btn_home.setObjectName(u"btn_home")
        sizePolicy.setHeightForWidth(self.btn_home.sizePolicy().hasHeightForWidth())
        self.btn_home.setSizePolicy(sizePolicy)
        self.btn_home.setMinimumSize(QSize(0, 45))
        self.btn_home.setFont(font1)
        self.btn_home.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.btn_home.setLayoutDirection(Qt.LeftToRight)
        self.btn_home.setStyleSheet(u"background-image: url(:/icons/images/icons/cil-home.png);\n"
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")
        self.btn_home.setCheckable(True)

        self.verticalLayout_8.addWidget(self.btn_home)

        self.btn_generate = QPushButton(self.topMenu)
        self.buttonGroup.addButton(self.btn_generate)
        self.btn_generate.setObjectName(u"btn_generate")
        sizePolicy.setHeightForWidth(self.btn_generate.sizePolicy().hasHeightForWidth())
        self.btn_generate.setSizePolicy(sizePolicy)
        self.btn_generate.setMinimumSize(QSize(0, 45))
        self.btn_generate.setFont(font1)
        self.btn_generate.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.btn_generate.setLayoutDirection(Qt.LeftToRight)
        self.btn_generate.setStyleSheet(u"background-image: url(:/icons/images/icons/cil-description.png);\n"
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")
        self.btn_generate.setCheckable(True)

        self.verticalLayout_8.addWidget(self.btn_generate)

        self.btn_figure = QPushButton(self.topMenu)
        self.buttonGroup.addButton(self.btn_figure)
        self.btn_figure.setObjectName(u"btn_figure")
        sizePolicy.setHeightForWidth(self.btn_figure.sizePolicy().hasHeightForWidth())
        self.btn_figure.setSizePolicy(sizePolicy)
        self.btn_figure.setMinimumSize(QSize(0, 45))
        self.btn_figure.setFont(font1)
        self.btn_figure.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.btn_figure.setLayoutDirection(Qt.LeftToRight)
        self.btn_figure.setStyleSheet(u"background-image: url(:/icons/images/icons/cil-image1.png);\n"
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")
        self.btn_figure.setCheckable(True)

        self.verticalLayout_8.addWidget(self.btn_figure)

        self.btn_analysis = QPushButton(self.topMenu)
        self.buttonGroup.addButton(self.btn_analysis)
        self.btn_analysis.setObjectName(u"btn_analysis")
        sizePolicy.setHeightForWidth(self.btn_analysis.sizePolicy().hasHeightForWidth())
        self.btn_analysis.setSizePolicy(sizePolicy)
        self.btn_analysis.setMinimumSize(QSize(0, 45))
        self.btn_analysis.setFont(font1)
        self.btn_analysis.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.btn_analysis.setLayoutDirection(Qt.LeftToRight)
        self.btn_analysis.setStyleSheet(u"background-image: url(:/icons/images/icons/cil-chart-line.png);\n"
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")
        self.btn_analysis.setCheckable(True)

        self.verticalLayout_8.addWidget(self.btn_analysis)

        self.btn_data_process = QPushButton(self.topMenu)
        self.buttonGroup.addButton(self.btn_data_process)
        self.btn_data_process.setObjectName(u"btn_data_process")
        sizePolicy.setHeightForWidth(self.btn_data_process.sizePolicy().hasHeightForWidth())
        self.btn_data_process.setSizePolicy(sizePolicy)
        self.btn_data_process.setMinimumSize(QSize(0, 45))
        self.btn_data_process.setFont(font1)
        self.btn_data_process.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.btn_data_process.setLayoutDirection(Qt.LeftToRight)
        self.btn_data_process.setStyleSheet(u"background-image:url(:/icons/images/icons/cil-window-restore.png);\n"
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")
        self.btn_data_process.setCheckable(True)

        self.verticalLayout_8.addWidget(self.btn_data_process)


        self.verticalMenuLayout.addWidget(self.topMenu, 0, Qt.AlignTop)

        self.bottomMenu = QFrame(self.leftMenuFrame)
        self.bottomMenu.setObjectName(u"bottomMenu")
        self.bottomMenu.setFrameShape(QFrame.NoFrame)
        self.bottomMenu.setFrameShadow(QFrame.Raised)
        self.verticalLayout_9 = QVBoxLayout(self.bottomMenu)
        self.verticalLayout_9.setSpacing(0)
        self.verticalLayout_9.setObjectName(u"verticalLayout_9")
        self.verticalLayout_9.setContentsMargins(0, 0, 0, 0)

        self.verticalMenuLayout.addWidget(self.bottomMenu, 0, Qt.AlignBottom)


        self.verticalLayout_3.addWidget(self.leftMenuFrame)


        self.appLayout.addWidget(self.leftMenuBg)

        self.contentBox = QFrame(self.bgApp)
        self.contentBox.setObjectName(u"contentBox")
        self.contentBox.setFrameShape(QFrame.NoFrame)
        self.contentBox.setFrameShadow(QFrame.Raised)
        self.verticalLayout_2 = QVBoxLayout(self.contentBox)
        self.verticalLayout_2.setSpacing(0)
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.contentTopBg = QFrame(self.contentBox)
        self.contentTopBg.setObjectName(u"contentTopBg")
        self.contentTopBg.setMinimumSize(QSize(0, 50))
        self.contentTopBg.setMaximumSize(QSize(16777215, 50))
        self.contentTopBg.setFrameShape(QFrame.NoFrame)
        self.contentTopBg.setFrameShadow(QFrame.Raised)
        self.horizontalLayout = QHBoxLayout(self.contentTopBg)
        self.horizontalLayout.setSpacing(0)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.horizontalLayout.setContentsMargins(0, 0, 10, 0)
        self.leftBox = QFrame(self.contentTopBg)
        self.leftBox.setObjectName(u"leftBox")
        sizePolicy1 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.leftBox.sizePolicy().hasHeightForWidth())
        self.leftBox.setSizePolicy(sizePolicy1)
        self.leftBox.setFrameShape(QFrame.NoFrame)
        self.leftBox.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_3 = QHBoxLayout(self.leftBox)
        self.horizontalLayout_3.setSpacing(0)
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.horizontalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.titleRightInfo = QLabel(self.leftBox)
        self.titleRightInfo.setObjectName(u"titleRightInfo")
        sizePolicy2 = QSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.titleRightInfo.sizePolicy().hasHeightForWidth())
        self.titleRightInfo.setSizePolicy(sizePolicy2)
        self.titleRightInfo.setMaximumSize(QSize(16777215, 45))
        font2 = QFont()
        font2.setFamilies([u"\u534e\u6587\u7ec6\u9ed1"])
        font2.setPointSize(16)
        font2.setBold(False)
        font2.setItalic(False)
        self.titleRightInfo.setFont(font2)
        self.titleRightInfo.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")
        self.titleRightInfo.setAlignment(Qt.AlignLeading|Qt.AlignLeft|Qt.AlignVCenter)

        self.horizontalLayout_3.addWidget(self.titleRightInfo)


        self.horizontalLayout.addWidget(self.leftBox)

        self.rightButtons = QFrame(self.contentTopBg)
        self.rightButtons.setObjectName(u"rightButtons")
        self.rightButtons.setMinimumSize(QSize(0, 28))
        self.rightButtons.setFrameShape(QFrame.NoFrame)
        self.rightButtons.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_2 = QHBoxLayout(self.rightButtons)
        self.horizontalLayout_2.setSpacing(5)
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.btn_language = QPushButton(self.rightButtons)
        self.btn_language.setObjectName(u"btn_language")
        self.btn_language.setMinimumSize(QSize(0, 0))
        self.btn_language.setMaximumSize(QSize(16000, 16000))
        self.btn_language.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.btn_language.setStyleSheet(u"\n"
"color:white;\n"
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.horizontalLayout_2.addWidget(self.btn_language)

        self.minimizeAppBtn = QPushButton(self.rightButtons)
        self.minimizeAppBtn.setObjectName(u"minimizeAppBtn")
        self.minimizeAppBtn.setMinimumSize(QSize(28, 28))
        self.minimizeAppBtn.setMaximumSize(QSize(28, 28))
        self.minimizeAppBtn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        icon = QIcon()
        icon.addFile(u":/icons/images/icons/icon_minimize.png", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.minimizeAppBtn.setIcon(icon)
        self.minimizeAppBtn.setIconSize(QSize(20, 20))

        self.horizontalLayout_2.addWidget(self.minimizeAppBtn)

        self.maximizeRestoreAppBtn = QPushButton(self.rightButtons)
        self.maximizeRestoreAppBtn.setObjectName(u"maximizeRestoreAppBtn")
        self.maximizeRestoreAppBtn.setMinimumSize(QSize(28, 28))
        self.maximizeRestoreAppBtn.setMaximumSize(QSize(28, 28))
        font3 = QFont()
        font3.setFamilies([u"Segoe UI"])
        font3.setPointSize(10)
        font3.setBold(False)
        font3.setItalic(False)
        font3.setStyleStrategy(QFont.PreferDefault)
        self.maximizeRestoreAppBtn.setFont(font3)
        self.maximizeRestoreAppBtn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        icon1 = QIcon()
        icon1.addFile(u":/icons/images/icons/icon_maximize.png", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.maximizeRestoreAppBtn.setIcon(icon1)
        self.maximizeRestoreAppBtn.setIconSize(QSize(20, 20))

        self.horizontalLayout_2.addWidget(self.maximizeRestoreAppBtn)

        self.closeAppBtn = QPushButton(self.rightButtons)
        self.closeAppBtn.setObjectName(u"closeAppBtn")
        self.closeAppBtn.setMinimumSize(QSize(28, 28))
        self.closeAppBtn.setMaximumSize(QSize(28, 28))
        self.closeAppBtn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        icon2 = QIcon()
        icon2.addFile(u":/icons/images/icons/icon_close.png", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.closeAppBtn.setIcon(icon2)
        self.closeAppBtn.setIconSize(QSize(20, 20))

        self.horizontalLayout_2.addWidget(self.closeAppBtn)


        self.horizontalLayout.addWidget(self.rightButtons, 0, Qt.AlignRight)


        self.verticalLayout_2.addWidget(self.contentTopBg)

        self.contentBottom = QFrame(self.contentBox)
        self.contentBottom.setObjectName(u"contentBottom")
        self.contentBottom.setFrameShape(QFrame.NoFrame)
        self.contentBottom.setFrameShadow(QFrame.Raised)
        self.verticalLayout_6 = QVBoxLayout(self.contentBottom)
        self.verticalLayout_6.setSpacing(0)
        self.verticalLayout_6.setObjectName(u"verticalLayout_6")
        self.verticalLayout_6.setContentsMargins(0, 0, 0, 0)
        self.content = QFrame(self.contentBottom)
        self.content.setObjectName(u"content")
        self.content.setFrameShape(QFrame.NoFrame)
        self.content.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_4 = QHBoxLayout(self.content)
        self.horizontalLayout_4.setSpacing(0)
        self.horizontalLayout_4.setObjectName(u"horizontalLayout_4")
        self.horizontalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.pagesContainer = QFrame(self.content)
        self.pagesContainer.setObjectName(u"pagesContainer")
        self.pagesContainer.setStyleSheet(u"")
        self.pagesContainer.setFrameShape(QFrame.NoFrame)
        self.pagesContainer.setFrameShadow(QFrame.Raised)
        self.verticalLayout_14 = QVBoxLayout(self.pagesContainer)
        self.verticalLayout_14.setObjectName(u"verticalLayout_14")
        self.stackedWidget = QStackedWidget(self.pagesContainer)
        self.stackedWidget.setObjectName(u"stackedWidget")
        palette = QPalette()
        brush = QBrush(QColor(221, 221, 221, 255))
        brush.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.WindowText, brush)
        brush1 = QBrush(QColor(0, 0, 0, 0))
        brush1.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Button, brush1)
        palette.setBrush(QPalette.Active, QPalette.Text, brush)
        palette.setBrush(QPalette.Active, QPalette.ButtonText, brush)
        palette.setBrush(QPalette.Active, QPalette.Base, brush1)
        palette.setBrush(QPalette.Active, QPalette.Window, brush1)
        palette.setBrush(QPalette.Inactive, QPalette.WindowText, brush)
        palette.setBrush(QPalette.Inactive, QPalette.Button, brush1)
        palette.setBrush(QPalette.Inactive, QPalette.Text, brush)
        palette.setBrush(QPalette.Inactive, QPalette.ButtonText, brush)
        palette.setBrush(QPalette.Inactive, QPalette.Base, brush1)
        palette.setBrush(QPalette.Inactive, QPalette.Window, brush1)
        palette.setBrush(QPalette.Disabled, QPalette.WindowText, brush)
        palette.setBrush(QPalette.Disabled, QPalette.Button, brush1)
        palette.setBrush(QPalette.Disabled, QPalette.Text, brush)
        palette.setBrush(QPalette.Disabled, QPalette.ButtonText, brush)
        palette.setBrush(QPalette.Disabled, QPalette.Base, brush1)
        palette.setBrush(QPalette.Disabled, QPalette.Window, brush1)
        self.stackedWidget.setPalette(palette)
        self.stackedWidget.setStyleSheet(u"background: transparent;")
        self.page_home = QWidget()
        self.page_home.setObjectName(u"page_home")
        self.page_home.setEnabled(True)
        self.page_home.setStyleSheet(u"")
        self.verticalLayout_26 = QVBoxLayout(self.page_home)
        self.verticalLayout_26.setObjectName(u"verticalLayout_26")
        self.label_5 = QLabel(self.page_home)
        self.label_5.setObjectName(u"label_5")
        self.label_5.setMinimumSize(QSize(0, 300))
        self.label_5.setStyleSheet(u"color: rgb(189,147,249);")

        self.verticalLayout_26.addWidget(self.label_5)

        self.stackedWidget.addWidget(self.page_home)
        self.page_analysis = QWidget()
        self.page_analysis.setObjectName(u"page_analysis")
        palette1 = QPalette()
        palette1.setBrush(QPalette.Active, QPalette.WindowText, brush)
        palette1.setBrush(QPalette.Active, QPalette.Button, brush1)
        palette1.setBrush(QPalette.Active, QPalette.Text, brush)
        palette1.setBrush(QPalette.Active, QPalette.ButtonText, brush)
        palette1.setBrush(QPalette.Active, QPalette.Base, brush1)
        palette1.setBrush(QPalette.Active, QPalette.Window, brush1)
        palette1.setBrush(QPalette.Inactive, QPalette.WindowText, brush)
        palette1.setBrush(QPalette.Inactive, QPalette.Button, brush1)
        palette1.setBrush(QPalette.Inactive, QPalette.Text, brush)
        palette1.setBrush(QPalette.Inactive, QPalette.ButtonText, brush)
        palette1.setBrush(QPalette.Inactive, QPalette.Base, brush1)
        palette1.setBrush(QPalette.Inactive, QPalette.Window, brush1)
        palette1.setBrush(QPalette.Disabled, QPalette.WindowText, brush)
        palette1.setBrush(QPalette.Disabled, QPalette.Button, brush1)
        palette1.setBrush(QPalette.Disabled, QPalette.Text, brush)
        palette1.setBrush(QPalette.Disabled, QPalette.ButtonText, brush)
        palette1.setBrush(QPalette.Disabled, QPalette.Base, brush1)
        palette1.setBrush(QPalette.Disabled, QPalette.Window, brush1)
        self.page_analysis.setPalette(palette1)
        self.horizontalLayout_11 = QHBoxLayout(self.page_analysis)
        self.horizontalLayout_11.setObjectName(u"horizontalLayout_11")
        self.widget_3 = QWidget(self.page_analysis)
        self.widget_3.setObjectName(u"widget_3")
        self.verticalLayout_10 = QVBoxLayout(self.widget_3)
        self.verticalLayout_10.setObjectName(u"verticalLayout_10")
        self.gridLayout_6 = QGridLayout()
        self.gridLayout_6.setObjectName(u"gridLayout_6")
        self.gridLayout_6.setContentsMargins(-1, -1, -1, 0)
        self.btnTrajectory = QPushButton(self.widget_3)
        self.btnTrajectory.setObjectName(u"btnTrajectory")
        self.btnTrajectory.setMinimumSize(QSize(75, 50))
        self.btnTrajectory.setMaximumSize(QSize(100, 16777215))
        font4 = QFont()
        font4.setFamilies([u"\u534e\u6587\u7ec6\u9ed1"])
        font4.setPointSize(12)
        font4.setBold(False)
        font4.setItalic(False)
        self.btnTrajectory.setFont(font4)
        self.btnTrajectory.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.btnTrajectory.setStyleSheet(u"background-color: rgb(52, 59, 72);\n"
"text-align:left;\n"
"font: 12pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"\n"
"")
        icon3 = QIcon()
        icon3.addFile(u":/icons/images/icons/cil-folder-open.png", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.btnTrajectory.setIcon(icon3)

        self.gridLayout_6.addWidget(self.btnTrajectory, 1, 1, 1, 1)

        self.editStructure = QLineEdit(self.widget_3)
        self.editStructure.setObjectName(u"editStructure")
        self.editStructure.setMinimumSize(QSize(100, 50))
        self.editStructure.setMaximumSize(QSize(16000, 200))
        self.editStructure.setStyleSheet(u"background-color: rgb(33, 37, 43);\n"
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_6.addWidget(self.editStructure, 0, 0, 1, 1)

        self.editTrajectory = QLineEdit(self.widget_3)
        self.editTrajectory.setObjectName(u"editTrajectory")
        self.editTrajectory.setMinimumSize(QSize(100, 50))
        self.editTrajectory.setMaximumSize(QSize(16777215, 200))
        self.editTrajectory.setStyleSheet(u"background-color: rgb(33, 37, 43);\n"
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_6.addWidget(self.editTrajectory, 1, 0, 1, 1)

        self.btnSructure = QPushButton(self.widget_3)
        self.btnSructure.setObjectName(u"btnSructure")
        self.btnSructure.setMinimumSize(QSize(75, 50))
        self.btnSructure.setMaximumSize(QSize(100, 16777215))
        self.btnSructure.setFont(font4)
        self.btnSructure.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.btnSructure.setStyleSheet(u"background-color: rgb(52, 59, 72);\n"
"text-align:left;\n"
"font: 12pt \"\u534e\u6587\u7ec6\u9ed1\";")
        self.btnSructure.setIcon(icon3)
        self.btnSructure.setCheckable(False)

        self.gridLayout_6.addWidget(self.btnSructure, 0, 1, 1, 1)

        self.btnResult = QPushButton(self.widget_3)
        self.btnResult.setObjectName(u"btnResult")
        self.btnResult.setMinimumSize(QSize(75, 50))
        self.btnResult.setMaximumSize(QSize(100, 16777215))
        self.btnResult.setFont(font4)
        self.btnResult.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.btnResult.setStyleSheet(u"background-color: rgb(52, 59, 72);\n"
"text-align:left;\n"
"font: 12pt \"\u534e\u6587\u7ec6\u9ed1\";")
        self.btnResult.setIcon(icon3)

        self.gridLayout_6.addWidget(self.btnResult, 2, 1, 1, 1)

        self.editResult = QLineEdit(self.widget_3)
        self.editResult.setObjectName(u"editResult")
        self.editResult.setMinimumSize(QSize(100, 50))
        self.editResult.setMaximumSize(QSize(16777215, 200))
        self.editResult.setStyleSheet(u"background-color: rgb(33, 37, 43);\n"
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_6.addWidget(self.editResult, 2, 0, 1, 1)


        self.verticalLayout_10.addLayout(self.gridLayout_6)

        self.row1_height = QFrame(self.widget_3)
        self.row1_height.setObjectName(u"row1_height")
        self.row1_height.setFrameShape(QFrame.NoFrame)
        self.row1_height.setFrameShadow(QFrame.Raised)
        self.verticalLayout_24 = QVBoxLayout(self.row1_height)
        self.verticalLayout_24.setObjectName(u"verticalLayout_24")
        self.row2_height_2 = QFrame(self.row1_height)
        self.row2_height_2.setObjectName(u"row2_height_2")
        self.row2_height_2.setFrameShape(QFrame.StyledPanel)
        self.row2_height_2.setFrameShadow(QFrame.Raised)
        self.verticalLayout_25 = QVBoxLayout(self.row2_height_2)
        self.verticalLayout_25.setObjectName(u"verticalLayout_25")
        self.gridLayout_7 = QGridLayout()
        self.gridLayout_7.setObjectName(u"gridLayout_7")
        self.comboBoxMethod = QComboBox(self.row2_height_2)
        self.comboBoxMethod.addItem("")
        self.comboBoxMethod.addItem("")
        self.comboBoxMethod.addItem("")
        self.comboBoxMethod.addItem("")
        self.comboBoxMethod.addItem("")
        self.comboBoxMethod.addItem("")
        self.comboBoxMethod.addItem("")
        self.comboBoxMethod.addItem("")
        self.comboBoxMethod.addItem("")
        self.comboBoxMethod.addItem("")
        self.comboBoxMethod.addItem("")
        self.comboBoxMethod.setObjectName(u"comboBoxMethod")
        self.comboBoxMethod.setMinimumSize(QSize(0, 40))
        font5 = QFont()
        font5.setFamilies([u"\u534e\u6587\u7ec6\u9ed1"])
        font5.setPointSize(12)
        font5.setBold(True)
        font5.setItalic(False)
        self.comboBoxMethod.setFont(font5)
        self.comboBoxMethod.setAutoFillBackground(False)
        self.comboBoxMethod.setStyleSheet(u"background-color:rgb(40,44,52);\n"
"font: 12pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"font-weight:bold;\n"
"border: 1px solid white;")
        self.comboBoxMethod.setIconSize(QSize(16, 16))
        self.comboBoxMethod.setFrame(True)

        self.gridLayout_7.addWidget(self.comboBoxMethod, 2, 5, 2, 1)

        self.editLastFrame = QSpinBox(self.row2_height_2)
        self.editLastFrame.setObjectName(u"editLastFrame")
        self.editLastFrame.setMinimumSize(QSize(0, 40))
        self.editLastFrame.setStyleSheet(u"QSpinBox {\n"
"    font: 12pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QSpinBox::up-button {\n"
"    height: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 30px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QSpinBox::down-button {\n"
"    height: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 30px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"")

        self.gridLayout_7.addWidget(self.editLastFrame, 2, 1, 2, 1)

        self.editStep = QSpinBox(self.row2_height_2)
        self.editStep.setObjectName(u"editStep")
        self.editStep.setMinimumSize(QSize(0, 40))
        self.editStep.setStyleSheet(u"QSpinBox {\n"
"    font: 12pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QSpinBox::up-button {\n"
"    height: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 30px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QSpinBox::down-button {\n"
"    height: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 30px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"")

        self.gridLayout_7.addWidget(self.editStep, 2, 2, 2, 1)

        self.line = QFrame(self.row2_height_2)
        self.line.setObjectName(u"line")
        self.line.setFrameShape(QFrame.Shape.VLine)
        self.line.setFrameShadow(QFrame.Shadow.Sunken)

        self.gridLayout_7.addWidget(self.line, 0, 3, 4, 1)

        self.analysis_label_frame = QLabel(self.row2_height_2)
        self.analysis_label_frame.setObjectName(u"analysis_label_frame")
        self.analysis_label_frame.setMinimumSize(QSize(0, 0))
        self.analysis_label_frame.setMaximumSize(QSize(16777215, 50))
        self.analysis_label_frame.setStyleSheet(u"font: 12pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_7.addWidget(self.analysis_label_frame, 0, 0, 1, 3)

        self.editFirstFrame = QSpinBox(self.row2_height_2)
        self.editFirstFrame.setObjectName(u"editFirstFrame")
        self.editFirstFrame.setMinimumSize(QSize(0, 40))
        self.editFirstFrame.setStyleSheet(u"QSpinBox {\n"
"    font: 12pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color: rgb(33, 37, 43);\n"
"}\n"
"\n"
"QSpinBox::up-button {\n"
"    height: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 30px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QSpinBox::down-button {\n"
"    height: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 30px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"")

        self.gridLayout_7.addWidget(self.editFirstFrame, 2, 0, 2, 1)

        self.analysis_label_method = QLabel(self.row2_height_2)
        self.analysis_label_method.setObjectName(u"analysis_label_method")
        self.analysis_label_method.setMinimumSize(QSize(0, 40))
        self.analysis_label_method.setMaximumSize(QSize(16777215, 100))
        self.analysis_label_method.setStyleSheet(u"font: 12pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_7.addWidget(self.analysis_label_method, 0, 5, 2, 1)

        self.analysis_label_k = QLabel(self.row2_height_2)
        self.analysis_label_k.setObjectName(u"analysis_label_k")
        self.analysis_label_k.setMinimumSize(QSize(0, 40))
        self.analysis_label_k.setMaximumSize(QSize(16777215, 100))
        self.analysis_label_k.setStyleSheet(u"font: 12pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_7.addWidget(self.analysis_label_k, 0, 4, 2, 1)

        self.analysis_label_first = QLabel(self.row2_height_2)
        self.analysis_label_first.setObjectName(u"analysis_label_first")
        self.analysis_label_first.setMaximumSize(QSize(16777215, 50))
        self.analysis_label_first.setStyleSheet(u"font-family:Verdana;\n"
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"color:white;")

        self.gridLayout_7.addWidget(self.analysis_label_first, 1, 0, 1, 1)

        self.analysis_label_last = QLabel(self.row2_height_2)
        self.analysis_label_last.setObjectName(u"analysis_label_last")
        self.analysis_label_last.setMaximumSize(QSize(16777215, 50))
        self.analysis_label_last.setStyleSheet(u"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_7.addWidget(self.analysis_label_last, 1, 1, 1, 1)

        self.analysis_label_step = QLabel(self.row2_height_2)
        self.analysis_label_step.setObjectName(u"analysis_label_step")
        self.analysis_label_step.setMaximumSize(QSize(16777215, 50))
        self.analysis_label_step.setStyleSheet(u"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_7.addWidget(self.analysis_label_step, 1, 2, 1, 1)

        self.editK = QSpinBox(self.row2_height_2)
        self.editK.setObjectName(u"editK")
        self.editK.setMinimumSize(QSize(0, 40))
        self.editK.setStyleSheet(u"QSpinBox {\n"
"    font: 12pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QSpinBox::up-button {\n"
"    height: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 30px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QSpinBox::down-button {\n"
"    height: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 30px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"")

        self.gridLayout_7.addWidget(self.editK, 2, 4, 2, 1)


        self.verticalLayout_25.addLayout(self.gridLayout_7)

        self.btnNext = QPushButton(self.row2_height_2)
        self.btnNext.setObjectName(u"btnNext")
        self.btnNext.setMinimumSize(QSize(150, 45))
        self.btnNext.setMaximumSize(QSize(16000, 50))
        self.btnNext.setFont(font)
        self.btnNext.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.btnNext.setLayoutDirection(Qt.LeftToRight)
        self.btnNext.setStyleSheet(u"background-color: rgb(189,147,249);\n"
"")
        icon4 = QIcon()
        icon4.addFile(u":/icons/images/icons/cil-chevron-double-right.png", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.btnNext.setIcon(icon4)

        self.verticalLayout_25.addWidget(self.btnNext)


        self.verticalLayout_24.addWidget(self.row2_height_2)


        self.verticalLayout_10.addWidget(self.row1_height)


        self.horizontalLayout_11.addWidget(self.widget_3)

        self.extraRightBox = QFrame(self.page_analysis)
        self.extraRightBox.setObjectName(u"extraRightBox")
        self.extraRightBox.setMaximumSize(QSize(0, 16777215))
        self.extraRightBox.setFrameShape(QFrame.StyledPanel)
        self.extraRightBox.setFrameShadow(QFrame.Raised)

        self.horizontalLayout_11.addWidget(self.extraRightBox)

        self.stackedWidget.addWidget(self.page_analysis)
        self.page_generation = QWidget()
        self.page_generation.setObjectName(u"page_generation")
        self.page_generation.setStyleSheet(u"b")
        self.horizontalLayout_6 = QHBoxLayout(self.page_generation)
        self.horizontalLayout_6.setObjectName(u"horizontalLayout_6")
        self.widget = QWidget(self.page_generation)
        self.widget.setObjectName(u"widget")
        self.verticalLayout = QVBoxLayout(self.widget)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.frame = QFrame(self.widget)
        self.frame.setObjectName(u"frame")
        self.frame.setFrameShape(QFrame.StyledPanel)
        self.frame.setFrameShadow(QFrame.Raised)
        self.verticalLayout_5 = QVBoxLayout(self.frame)
        self.verticalLayout_5.setObjectName(u"verticalLayout_5")
        self.label_gene_label = QLabel(self.frame)
        self.label_gene_label.setObjectName(u"label_gene_label")
        self.label_gene_label.setMaximumSize(QSize(16777215, 35))
        self.label_gene_label.setFont(font2)
        self.label_gene_label.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.verticalLayout_5.addWidget(self.label_gene_label)

        self.frame_5 = QFrame(self.frame)
        self.frame_5.setObjectName(u"frame_5")
        self.frame_5.setFrameShape(QFrame.StyledPanel)
        self.frame_5.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_7 = QHBoxLayout(self.frame_5)
        self.horizontalLayout_7.setObjectName(u"horizontalLayout_7")
        self.edit_gene_path = QLineEdit(self.frame_5)
        self.edit_gene_path.setObjectName(u"edit_gene_path")
        self.edit_gene_path.setMinimumSize(QSize(0, 45))
        self.edit_gene_path.setStyleSheet(u"background-color: rgb(33, 37, 43);\n"
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.horizontalLayout_7.addWidget(self.edit_gene_path)

        self.btn_gene_path = QPushButton(self.frame_5)
        self.btn_gene_path.setObjectName(u"btn_gene_path")
        self.btn_gene_path.setMinimumSize(QSize(150, 45))
        self.btn_gene_path.setFont(font4)
        self.btn_gene_path.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.btn_gene_path.setStyleSheet(u"background-color: rgb(52, 59, 72);\n"
"text-align:left;\n"
"font: 12pt \"\u534e\u6587\u7ec6\u9ed1\";")
        self.btn_gene_path.setIcon(icon3)

        self.horizontalLayout_7.addWidget(self.btn_gene_path)


        self.verticalLayout_5.addWidget(self.frame_5)


        self.verticalLayout.addWidget(self.frame)

        self.frame_2 = QFrame(self.widget)
        self.frame_2.setObjectName(u"frame_2")
        self.frame_2.setFrameShape(QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QFrame.Raised)
        self.gridLayout_2 = QGridLayout(self.frame_2)
        self.gridLayout_2.setObjectName(u"gridLayout_2")
        self.gridLayout = QGridLayout()
        self.gridLayout.setObjectName(u"gridLayout")
        self.label_gene_y = QLabel(self.frame_2)
        self.label_gene_y.setObjectName(u"label_gene_y")
        self.label_gene_y.setMinimumSize(QSize(0, 0))
        self.label_gene_y.setMaximumSize(QSize(16777215, 30))
        self.label_gene_y.setStyleSheet(u"font: 15pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout.addWidget(self.label_gene_y, 0, 1, 1, 1)

        self.label_gene_x = QLabel(self.frame_2)
        self.label_gene_x.setObjectName(u"label_gene_x")
        self.label_gene_x.setMinimumSize(QSize(0, 0))
        self.label_gene_x.setMaximumSize(QSize(16777215, 30))
        self.label_gene_x.setStyleSheet(u"font: 15pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout.addWidget(self.label_gene_x, 0, 0, 1, 1)

        self.label_gene_z = QLabel(self.frame_2)
        self.label_gene_z.setObjectName(u"label_gene_z")
        self.label_gene_z.setMinimumSize(QSize(0, 0))
        self.label_gene_z.setStyleSheet(u"font: 15pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout.addWidget(self.label_gene_z, 0, 2, 1, 1)

        self.spin_box_x = QDoubleSpinBox(self.frame_2)
        self.spin_box_x.setObjectName(u"spin_box_x")
        self.spin_box_x.setMinimumSize(QSize(0, 30))
        self.spin_box_x.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"")
        self.spin_box_x.setValue(20.000000000000000)

        self.gridLayout.addWidget(self.spin_box_x, 1, 0, 1, 1)

        self.spin_box_y = QDoubleSpinBox(self.frame_2)
        self.spin_box_y.setObjectName(u"spin_box_y")
        self.spin_box_y.setMinimumSize(QSize(0, 30))
        self.spin_box_y.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"")
        self.spin_box_y.setValue(20.000000000000000)

        self.gridLayout.addWidget(self.spin_box_y, 1, 1, 1, 1)

        self.spin_box_z = QDoubleSpinBox(self.frame_2)
        self.spin_box_z.setObjectName(u"spin_box_z")
        self.spin_box_z.setMinimumSize(QSize(0, 30))
        self.spin_box_z.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"")
        self.spin_box_z.setValue(20.000000000000000)

        self.gridLayout.addWidget(self.spin_box_z, 1, 2, 1, 1)


        self.gridLayout_2.addLayout(self.gridLayout, 0, 0, 1, 1)


        self.verticalLayout.addWidget(self.frame_2)

        self.frame_3 = QFrame(self.widget)
        self.frame_3.setObjectName(u"frame_3")
        self.frame_3.setFrameShape(QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QFrame.Raised)
        self.gridLayout_8 = QGridLayout(self.frame_3)
        self.gridLayout_8.setObjectName(u"gridLayout_8")
        self.label_gene_gas = QLabel(self.frame_3)
        self.label_gene_gas.setObjectName(u"label_gene_gas")
        self.label_gene_gas.setMaximumSize(QSize(16777215, 30))
        self.label_gene_gas.setStyleSheet(u"font: 15pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_8.addWidget(self.label_gene_gas, 0, 0, 1, 1)

        self.como_gas = QComboBox(self.frame_3)
        self.como_gas.addItem("")
        self.como_gas.setObjectName(u"como_gas")
        self.como_gas.setMinimumSize(QSize(30, 0))
        self.como_gas.setStyleSheet(u"background-color:rgb(40,44,52);\n"
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"font-weight:bold;\n"
"border: 1px solid white;")

        self.gridLayout_8.addWidget(self.como_gas, 1, 0, 1, 1)

        self.label_gene_gden = QLabel(self.frame_3)
        self.label_gene_gden.setObjectName(u"label_gene_gden")
        self.label_gene_gden.setStyleSheet(u"font: 15pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_8.addWidget(self.label_gene_gden, 0, 1, 1, 1)

        self.label_gene_a = QLabel(self.frame_3)
        self.label_gene_a.setObjectName(u"label_gene_a")
        self.label_gene_a.setStyleSheet(u"font: 15pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_8.addWidget(self.label_gene_a, 0, 2, 1, 1)

        self.spin_gas_density = QSpinBox(self.frame_3)
        self.spin_gas_density.setObjectName(u"spin_gas_density")
        self.spin_gas_density.setMinimumSize(QSize(0, 30))
        self.spin_gas_density.setStyleSheet(u"QSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"")
        self.spin_gas_density.setMaximum(100000)
        self.spin_gas_density.setValue(145)

        self.gridLayout_8.addWidget(self.spin_gas_density, 1, 1, 1, 1)

        self.spin_area_5 = QDoubleSpinBox(self.frame_3)
        self.spin_area_5.setObjectName(u"spin_area_5")
        self.spin_area_5.setMinimumSize(QSize(0, 30))
        self.spin_area_5.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"")
        self.spin_area_5.setValue(1.000000000000000)

        self.gridLayout_8.addWidget(self.spin_area_5, 1, 2, 1, 1)

        self.label_gene_r = QLabel(self.frame_3)
        self.label_gene_r.setObjectName(u"label_gene_r")
        self.label_gene_r.setStyleSheet(u"font: 15pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_8.addWidget(self.label_gene_r, 0, 3, 1, 1)

        self.spin_r = QDoubleSpinBox(self.frame_3)
        self.spin_r.setObjectName(u"spin_r")
        self.spin_r.setMinimumSize(QSize(0, 30))
        self.spin_r.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"")
        self.spin_r.setValue(5.000000000000000)

        self.gridLayout_8.addWidget(self.spin_r, 1, 3, 1, 1)

        self.btn_gene_lipid = QPushButton(self.frame_3)
        self.btn_gene_lipid.setObjectName(u"btn_gene_lipid")
        self.btn_gene_lipid.setMinimumSize(QSize(0, 80))
        self.btn_gene_lipid.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.btn_gene_lipid.setStyleSheet(u"font: 18pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"background-color: rgb(189,147,249);\n"
"color:white;")

        self.gridLayout_8.addWidget(self.btn_gene_lipid, 0, 4, 2, 1)


        self.verticalLayout.addWidget(self.frame_3)

        self.frame_4 = QFrame(self.widget)
        self.frame_4.setObjectName(u"frame_4")
        self.frame_4.setFrameShape(QFrame.StyledPanel)
        self.frame_4.setFrameShadow(QFrame.Raised)
        self.gridLayout_10 = QGridLayout(self.frame_4)
        self.gridLayout_10.setObjectName(u"gridLayout_10")
        self.spin_salt = QDoubleSpinBox(self.frame_4)
        self.spin_salt.setObjectName(u"spin_salt")
        self.spin_salt.setMinimumSize(QSize(0, 30))
        self.spin_salt.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"")
        self.spin_salt.setValue(0.150000000000000)

        self.gridLayout_10.addWidget(self.spin_salt, 2, 1, 1, 1)

        self.label_gene_sol = QLabel(self.frame_4)
        self.label_gene_sol.setObjectName(u"label_gene_sol")
        self.label_gene_sol.setMaximumSize(QSize(16777215, 30))
        self.label_gene_sol.setStyleSheet(u"font: 15pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_10.addWidget(self.label_gene_sol, 0, 0, 1, 1)

        self.label_gene_salt = QLabel(self.frame_4)
        self.label_gene_salt.setObjectName(u"label_gene_salt")
        self.label_gene_salt.setStyleSheet(u"font: 15pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_10.addWidget(self.label_gene_salt, 0, 1, 1, 1)

        self.como_solvent = QComboBox(self.frame_4)
        self.como_solvent.addItem("")
        self.como_solvent.setObjectName(u"como_solvent")
        self.como_solvent.setMinimumSize(QSize(0, 30))
        self.como_solvent.setStyleSheet(u"background-color:rgb(40,44,52);\n"
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"font-weight:bold;\n"
"border: 1px solid white;")

        self.gridLayout_10.addWidget(self.como_solvent, 2, 0, 1, 1)


        self.verticalLayout.addWidget(self.frame_4)

        self.btn_gene_run = QPushButton(self.widget)
        self.btn_gene_run.setObjectName(u"btn_gene_run")
        self.btn_gene_run.setMinimumSize(QSize(0, 30))
        self.btn_gene_run.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.btn_gene_run.setStyleSheet(u"background-color: rgb(189,147,249);\n"
"color:white;\n"
"font-size:18pt;")

        self.verticalLayout.addWidget(self.btn_gene_run)


        self.horizontalLayout_6.addWidget(self.widget)

        self.Gene_extraBox = QFrame(self.page_generation)
        self.Gene_extraBox.setObjectName(u"Gene_extraBox")
        self.Gene_extraBox.setMaximumSize(QSize(0, 16777215))
        self.Gene_extraBox.setFrameShape(QFrame.StyledPanel)
        self.Gene_extraBox.setFrameShadow(QFrame.Raised)

        self.horizontalLayout_6.addWidget(self.Gene_extraBox)

        self.stackedWidget.addWidget(self.page_generation)
        self.page_figure = QWidget()
        self.page_figure.setObjectName(u"page_figure")
        self.horizontalLayout_8 = QHBoxLayout(self.page_figure)
        self.horizontalLayout_8.setObjectName(u"horizontalLayout_8")
        self.widget_2 = QWidget(self.page_figure)
        self.widget_2.setObjectName(u"widget_2")
        self.widget_2.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")
        self.verticalLayout_11 = QVBoxLayout(self.widget_2)
        self.verticalLayout_11.setObjectName(u"verticalLayout_11")
        self.frame_8 = QFrame(self.widget_2)
        self.frame_8.setObjectName(u"frame_8")
        self.frame_8.setFrameShape(QFrame.StyledPanel)
        self.frame_8.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_9 = QHBoxLayout(self.frame_8)
        self.horizontalLayout_9.setObjectName(u"horizontalLayout_9")
        self.figure_edit_path = QLineEdit(self.frame_8)
        self.figure_edit_path.setObjectName(u"figure_edit_path")
        self.figure_edit_path.setMinimumSize(QSize(0, 45))
        self.figure_edit_path.setStyleSheet(u"background-color: rgb(33, 37, 43);\n"
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")
        self.figure_edit_path.setReadOnly(True)

        self.horizontalLayout_9.addWidget(self.figure_edit_path)

        self.figure_btn_path = QPushButton(self.frame_8)
        self.figure_btn_path.setObjectName(u"figure_btn_path")
        self.figure_btn_path.setMinimumSize(QSize(80, 45))
        self.figure_btn_path.setStyleSheet(u"background-color: rgb(52, 59, 72);\n"
"text-align:left;\n"
"font: 12pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.horizontalLayout_9.addWidget(self.figure_btn_path)


        self.verticalLayout_11.addWidget(self.frame_8)

        self.tabWidget = QTabWidget(self.widget_2)
        self.tabWidget.setObjectName(u"tabWidget")
        self.tabWidget.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.tabWidget.setStyleSheet(u"QTabBar::tab {\n"
"    background: lightgray;\n"
"    border: 2px solid #C4C4C3;\n"
"    border-bottom-color: #C4C4C3; /* same as the pane color */\n"
"    border-top-left-radius: 5px;\n"
"    border-top-right-radius: 5px;\n"
"    min-width: 16ex;\n"
"    padding: 2px;\n"
"   font: 18pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"   color:balck;\n"
"}\n"
"\n"
"QTabBar::tab:selected {\n"
"    background: lightblue;\n"
"}\n"
"\n"
"QTabBar::tab:hover {\n"
"    background: pink;\n"
"}\n"
"")
        self.tab_4 = QWidget()
        self.tab_4.setObjectName(u"tab_4")
        self.gridLayout_9 = QGridLayout(self.tab_4)
        self.gridLayout_9.setObjectName(u"gridLayout_9")
        self.figure_line_label_legend = QLabel(self.tab_4)
        self.figure_line_label_legend.setObjectName(u"figure_line_label_legend")
        self.figure_line_label_legend.setStyleSheet(u"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_9.addWidget(self.figure_line_label_legend, 2, 0, 1, 1)

        self.figure_line_btn_color_2 = QPushButton(self.tab_4)
        self.figure_line_btn_color_2.setObjectName(u"figure_line_btn_color_2")
        self.figure_line_btn_color_2.setEnabled(True)
        self.figure_line_btn_color_2.setMinimumSize(QSize(0, 38))
        self.figure_line_btn_color_2.setMaximumSize(QSize(16000000, 16777215))
        self.figure_line_btn_color_2.setBaseSize(QSize(0, 0))
        self.figure_line_btn_color_2.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.figure_line_btn_color_2.setStyleSheet(u"background-color: rgb(189,147,249);\n"
"color:white;\n"
"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_9.addWidget(self.figure_line_btn_color_2, 8, 1, 1, 1)

        self.figure_line_spin_marker_size_2 = QDoubleSpinBox(self.tab_4)
        self.figure_line_spin_marker_size_2.setObjectName(u"figure_line_spin_marker_size_2")
        self.figure_line_spin_marker_size_2.setMinimumSize(QSize(0, 30))
        self.figure_line_spin_marker_size_2.setMaximumSize(QSize(200, 16777215))
        self.figure_line_spin_marker_size_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")

        self.gridLayout_9.addWidget(self.figure_line_spin_marker_size_2, 7, 1, 1, 1)

        self.figure_line_label_x = QLabel(self.tab_4)
        self.figure_line_label_x.setObjectName(u"figure_line_label_x")
        self.figure_line_label_x.setStyleSheet(u"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_9.addWidget(self.figure_line_label_x, 3, 0, 1, 1)

        self.figure_line_label_color = QLabel(self.tab_4)
        self.figure_line_label_color.setObjectName(u"figure_line_label_color")
        self.figure_line_label_color.setStyleSheet(u"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_9.addWidget(self.figure_line_label_color, 8, 0, 1, 1)

        self.figure_line_label_y_range = QLabel(self.tab_4)
        self.figure_line_label_y_range.setObjectName(u"figure_line_label_y_range")
        self.figure_line_label_y_range.setStyleSheet(u"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_9.addWidget(self.figure_line_label_y_range, 6, 0, 1, 1)

        self.figure_line_label_x_range = QLabel(self.tab_4)
        self.figure_line_label_x_range.setObjectName(u"figure_line_label_x_range")
        self.figure_line_label_x_range.setStyleSheet(u"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_9.addWidget(self.figure_line_label_x_range, 5, 0, 1, 1)

        self.figure_line_label_y = QLabel(self.tab_4)
        self.figure_line_label_y.setObjectName(u"figure_line_label_y")
        self.figure_line_label_y.setStyleSheet(u"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_9.addWidget(self.figure_line_label_y, 4, 0, 1, 1)

        self.figure_line_label_axis_tick = QLabel(self.tab_4)
        self.figure_line_label_axis_tick.setObjectName(u"figure_line_label_axis_tick")
        self.figure_line_label_axis_tick.setStyleSheet(u"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_9.addWidget(self.figure_line_label_axis_tick, 0, 0, 1, 1)

        self.figure_line_label_marker = QLabel(self.tab_4)
        self.figure_line_label_marker.setObjectName(u"figure_line_label_marker")
        self.figure_line_label_marker.setStyleSheet(u"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_9.addWidget(self.figure_line_label_marker, 7, 0, 1, 1)

        self.figure_line_label_axis_title = QLabel(self.tab_4)
        self.figure_line_label_axis_title.setObjectName(u"figure_line_label_axis_title")
        self.figure_line_label_axis_title.setMinimumSize(QSize(50, 0))
        self.figure_line_label_axis_title.setStyleSheet(u"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_9.addWidget(self.figure_line_label_axis_title, 1, 0, 1, 1)

        self.figure_line_spin_axis_scale_2 = QDoubleSpinBox(self.tab_4)
        self.figure_line_spin_axis_scale_2.setObjectName(u"figure_line_spin_axis_scale_2")
        self.figure_line_spin_axis_scale_2.setMinimumSize(QSize(0, 30))
        self.figure_line_spin_axis_scale_2.setMaximumSize(QSize(200, 16777215))
        self.figure_line_spin_axis_scale_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")
        self.figure_line_spin_axis_scale_2.setMaximum(100000.000000000000000)
        self.figure_line_spin_axis_scale_2.setValue(10.000000000000000)

        self.gridLayout_9.addWidget(self.figure_line_spin_axis_scale_2, 0, 1, 1, 1)

        self.figure_line_spin_grid_size_2 = QDoubleSpinBox(self.tab_4)
        self.figure_line_spin_grid_size_2.setObjectName(u"figure_line_spin_grid_size_2")
        self.figure_line_spin_grid_size_2.setMinimumSize(QSize(0, 30))
        self.figure_line_spin_grid_size_2.setMaximumSize(QSize(200, 16777215))
        self.figure_line_spin_grid_size_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")
        self.figure_line_spin_grid_size_2.setMaximum(100000.000000000000000)
        self.figure_line_spin_grid_size_2.setValue(10.000000000000000)

        self.gridLayout_9.addWidget(self.figure_line_spin_grid_size_2, 2, 1, 1, 1)

        self.figure_line_spin_axis_text_size_2 = QDoubleSpinBox(self.tab_4)
        self.figure_line_spin_axis_text_size_2.setObjectName(u"figure_line_spin_axis_text_size_2")
        self.figure_line_spin_axis_text_size_2.setMinimumSize(QSize(200, 30))
        self.figure_line_spin_axis_text_size_2.setMaximumSize(QSize(200, 16777215))
        self.figure_line_spin_axis_text_size_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")
        self.figure_line_spin_axis_text_size_2.setMaximum(100000.000000000000000)
        self.figure_line_spin_axis_text_size_2.setValue(11.000000000000000)

        self.gridLayout_9.addWidget(self.figure_line_spin_axis_text_size_2, 1, 1, 1, 1)

        self.figure_line_spin_x_min_2 = QDoubleSpinBox(self.tab_4)
        self.figure_line_spin_x_min_2.setObjectName(u"figure_line_spin_x_min_2")
        self.figure_line_spin_x_min_2.setMinimumSize(QSize(0, 30))
        self.figure_line_spin_x_min_2.setMaximumSize(QSize(200, 16777215))
        self.figure_line_spin_x_min_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")
        self.figure_line_spin_x_min_2.setMinimum(-1000000000.000000000000000)
        self.figure_line_spin_x_min_2.setMaximum(1000000000.000000000000000)

        self.gridLayout_9.addWidget(self.figure_line_spin_x_min_2, 5, 1, 1, 1)

        self.figure_line_edit_y_2 = QLineEdit(self.tab_4)
        self.figure_line_edit_y_2.setObjectName(u"figure_line_edit_y_2")
        self.figure_line_edit_y_2.setMinimumSize(QSize(0, 30))
        self.figure_line_edit_y_2.setMaximumSize(QSize(16000000, 16777215))
        self.figure_line_edit_y_2.setStyleSheet(u"background-color: rgb(33, 37, 43);\n"
"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_9.addWidget(self.figure_line_edit_y_2, 4, 1, 1, 1)

        self.figure_line_spin_y_min_2 = QDoubleSpinBox(self.tab_4)
        self.figure_line_spin_y_min_2.setObjectName(u"figure_line_spin_y_min_2")
        self.figure_line_spin_y_min_2.setMinimumSize(QSize(0, 30))
        self.figure_line_spin_y_min_2.setMaximumSize(QSize(200, 16777215))
        self.figure_line_spin_y_min_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")
        self.figure_line_spin_y_min_2.setMinimum(-1000000000.000000000000000)
        self.figure_line_spin_y_min_2.setMaximum(1000000000.000000000000000)

        self.gridLayout_9.addWidget(self.figure_line_spin_y_min_2, 6, 1, 1, 1)

        self.figure_line_spin_x_max_2 = QDoubleSpinBox(self.tab_4)
        self.figure_line_spin_x_max_2.setObjectName(u"figure_line_spin_x_max_2")
        self.figure_line_spin_x_max_2.setMinimumSize(QSize(200, 30))
        self.figure_line_spin_x_max_2.setMaximumSize(QSize(16000000, 16777215))
        self.figure_line_spin_x_max_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")
        self.figure_line_spin_x_max_2.setMinimum(-1000000000.000000000000000)
        self.figure_line_spin_x_max_2.setMaximum(1000000000.000000000000000)

        self.gridLayout_9.addWidget(self.figure_line_spin_x_max_2, 5, 2, 1, 1)

        self.figure_line_como_marker = QComboBox(self.tab_4)
        self.figure_line_como_marker.addItem("")
        self.figure_line_como_marker.addItem("")
        self.figure_line_como_marker.addItem("")
        self.figure_line_como_marker.addItem("")
        self.figure_line_como_marker.addItem("")
        self.figure_line_como_marker.addItem("")
        self.figure_line_como_marker.addItem("")
        self.figure_line_como_marker.setObjectName(u"figure_line_como_marker")
        self.figure_line_como_marker.setMinimumSize(QSize(0, 30))
        self.figure_line_como_marker.setStyleSheet(u"background-color:rgb(40,44,52);\n"
"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"font-weight:bold;\n"
"border: 1px solid white;")

        self.gridLayout_9.addWidget(self.figure_line_como_marker, 7, 2, 1, 1)

        self.figure_line_edit_x_2 = QLineEdit(self.tab_4)
        self.figure_line_edit_x_2.setObjectName(u"figure_line_edit_x_2")
        self.figure_line_edit_x_2.setMinimumSize(QSize(0, 30))
        self.figure_line_edit_x_2.setMaximumSize(QSize(16777215, 16777215))
        self.figure_line_edit_x_2.setStyleSheet(u"background-color: rgb(33, 37, 43);\n"
"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_9.addWidget(self.figure_line_edit_x_2, 3, 1, 1, 1)

        self.figure_line_spin_y_max_2 = QDoubleSpinBox(self.tab_4)
        self.figure_line_spin_y_max_2.setObjectName(u"figure_line_spin_y_max_2")
        self.figure_line_spin_y_max_2.setMinimumSize(QSize(200, 30))
        self.figure_line_spin_y_max_2.setMaximumSize(QSize(16000000, 16777215))
        self.figure_line_spin_y_max_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")
        self.figure_line_spin_y_max_2.setMinimum(-1000000000.000000000000000)
        self.figure_line_spin_y_max_2.setMaximum(1000000000.000000000000000)

        self.gridLayout_9.addWidget(self.figure_line_spin_y_max_2, 6, 2, 1, 1)

        self.tabWidget.addTab(self.tab_4, "")
        self.tab_5 = QWidget()
        self.tab_5.setObjectName(u"tab_5")
        self.gridLayout_11 = QGridLayout(self.tab_5)
        self.gridLayout_11.setObjectName(u"gridLayout_11")
        self.figure_bar_btn_color_2 = QPushButton(self.tab_5)
        self.figure_bar_btn_color_2.setObjectName(u"figure_bar_btn_color_2")
        self.figure_bar_btn_color_2.setMinimumSize(QSize(0, 30))
        self.figure_bar_btn_color_2.setStyleSheet(u"background-color: rgb(189,147,249);\n"
"color:white;\n"
"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_11.addWidget(self.figure_bar_btn_color_2, 7, 2, 1, 1)

        self.figure_bar_spin_bar_2 = QDoubleSpinBox(self.tab_5)
        self.figure_bar_spin_bar_2.setObjectName(u"figure_bar_spin_bar_2")
        self.figure_bar_spin_bar_2.setMinimumSize(QSize(0, 30))
        self.figure_bar_spin_bar_2.setMaximumSize(QSize(200, 16777215))
        self.figure_bar_spin_bar_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")

        self.gridLayout_11.addWidget(self.figure_bar_spin_bar_2, 6, 1, 1, 1)

        self.figure_bar_spin_trend_2 = QDoubleSpinBox(self.tab_5)
        self.figure_bar_spin_trend_2.setObjectName(u"figure_bar_spin_trend_2")
        self.figure_bar_spin_trend_2.setMinimumSize(QSize(0, 30))
        self.figure_bar_spin_trend_2.setMaximumSize(QSize(200, 16777215))
        self.figure_bar_spin_trend_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")

        self.gridLayout_11.addWidget(self.figure_bar_spin_trend_2, 5, 1, 1, 1)

        self.figure_bar_btn_trend_2 = QPushButton(self.tab_5)
        self.figure_bar_btn_trend_2.setObjectName(u"figure_bar_btn_trend_2")
        self.figure_bar_btn_trend_2.setMinimumSize(QSize(0, 30))
        self.figure_bar_btn_trend_2.setStyleSheet(u"background-color: rgb(189,147,249);\n"
"color:white;\n"
"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_11.addWidget(self.figure_bar_btn_trend_2, 5, 2, 1, 1)

        self.figure_bar_label_trend = QLabel(self.tab_5)
        self.figure_bar_label_trend.setObjectName(u"figure_bar_label_trend")
        self.figure_bar_label_trend.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_11.addWidget(self.figure_bar_label_trend, 5, 0, 1, 1)

        self.figure_bar_label_y = QLabel(self.tab_5)
        self.figure_bar_label_y.setObjectName(u"figure_bar_label_y")
        self.figure_bar_label_y.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_11.addWidget(self.figure_bar_label_y, 3, 0, 1, 1)

        self.figure_bar_edit_x_2 = QLineEdit(self.tab_5)
        self.figure_bar_edit_x_2.setObjectName(u"figure_bar_edit_x_2")
        self.figure_bar_edit_x_2.setMinimumSize(QSize(0, 30))
        self.figure_bar_edit_x_2.setStyleSheet(u"background-color: rgb(33, 37, 43);\n"
"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_11.addWidget(self.figure_bar_edit_x_2, 2, 1, 1, 1)

        self.figure_bar_spin_axis_scale_2 = QDoubleSpinBox(self.tab_5)
        self.figure_bar_spin_axis_scale_2.setObjectName(u"figure_bar_spin_axis_scale_2")
        self.figure_bar_spin_axis_scale_2.setMinimumSize(QSize(0, 30))
        self.figure_bar_spin_axis_scale_2.setMaximumSize(QSize(200, 16777215))
        self.figure_bar_spin_axis_scale_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")
        self.figure_bar_spin_axis_scale_2.setMaximum(10000.000000000000000)
        self.figure_bar_spin_axis_scale_2.setValue(10.000000000000000)

        self.gridLayout_11.addWidget(self.figure_bar_spin_axis_scale_2, 0, 1, 1, 1)

        self.figure_bar_edit_y_2 = QLineEdit(self.tab_5)
        self.figure_bar_edit_y_2.setObjectName(u"figure_bar_edit_y_2")
        self.figure_bar_edit_y_2.setMinimumSize(QSize(0, 30))
        self.figure_bar_edit_y_2.setStyleSheet(u"background-color: rgb(33, 37, 43);\n"
"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_11.addWidget(self.figure_bar_edit_y_2, 3, 1, 1, 1)

        self.figure_bar_spin_axis_text_size_2 = QDoubleSpinBox(self.tab_5)
        self.figure_bar_spin_axis_text_size_2.setObjectName(u"figure_bar_spin_axis_text_size_2")
        self.figure_bar_spin_axis_text_size_2.setMinimumSize(QSize(0, 30))
        self.figure_bar_spin_axis_text_size_2.setMaximumSize(QSize(200, 16777215))
        self.figure_bar_spin_axis_text_size_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")
        self.figure_bar_spin_axis_text_size_2.setMaximum(10000.000000000000000)
        self.figure_bar_spin_axis_text_size_2.setValue(10.000000000000000)

        self.gridLayout_11.addWidget(self.figure_bar_spin_axis_text_size_2, 1, 1, 1, 1)

        self.figure_bar_label_axis_title = QLabel(self.tab_5)
        self.figure_bar_label_axis_title.setObjectName(u"figure_bar_label_axis_title")
        self.figure_bar_label_axis_title.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_11.addWidget(self.figure_bar_label_axis_title, 1, 0, 1, 1)

        self.figure_bar_label_error = QLabel(self.tab_5)
        self.figure_bar_label_error.setObjectName(u"figure_bar_label_error")
        self.figure_bar_label_error.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_11.addWidget(self.figure_bar_label_error, 8, 0, 1, 1)

        self.figure_bar_label_x = QLabel(self.tab_5)
        self.figure_bar_label_x.setObjectName(u"figure_bar_label_x")
        self.figure_bar_label_x.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_11.addWidget(self.figure_bar_label_x, 2, 0, 1, 1)

        self.figure_bar_label_axis_tick = QLabel(self.tab_5)
        self.figure_bar_label_axis_tick.setObjectName(u"figure_bar_label_axis_tick")
        self.figure_bar_label_axis_tick.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_11.addWidget(self.figure_bar_label_axis_tick, 0, 0, 1, 1)

        self.figure_bar_label_color = QLabel(self.tab_5)
        self.figure_bar_label_color.setObjectName(u"figure_bar_label_color")
        self.figure_bar_label_color.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_11.addWidget(self.figure_bar_label_color, 7, 0, 1, 1)

        self.figure_bar_label_bar = QLabel(self.tab_5)
        self.figure_bar_label_bar.setObjectName(u"figure_bar_label_bar")
        self.figure_bar_label_bar.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_11.addWidget(self.figure_bar_label_bar, 6, 0, 1, 1)

        self.figure_bar_label_y_range = QLabel(self.tab_5)
        self.figure_bar_label_y_range.setObjectName(u"figure_bar_label_y_range")
        self.figure_bar_label_y_range.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_11.addWidget(self.figure_bar_label_y_range, 4, 0, 1, 1)

        self.figure_bar_radio_error_2 = QRadioButton(self.tab_5)
        self.figure_bar_radio_error_2.setObjectName(u"figure_bar_radio_error_2")
        self.figure_bar_radio_error_2.setMinimumSize(QSize(0, 30))
        self.figure_bar_radio_error_2.setStyleSheet(u"font-size: 16pt \"\u534e\u6587\u7ec6\u9ed1\"")

        self.gridLayout_11.addWidget(self.figure_bar_radio_error_2, 8, 1, 1, 1)

        self.radioButton_6 = QRadioButton(self.tab_5)
        self.radioButton_6.setObjectName(u"radioButton_6")
        self.radioButton_6.setMinimumSize(QSize(0, 30))
        self.radioButton_6.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\"")
        self.radioButton_6.setChecked(True)

        self.gridLayout_11.addWidget(self.radioButton_6, 8, 2, 1, 1)

        self.figure_bar_spin_y_min_2 = QDoubleSpinBox(self.tab_5)
        self.figure_bar_spin_y_min_2.setObjectName(u"figure_bar_spin_y_min_2")
        self.figure_bar_spin_y_min_2.setMinimumSize(QSize(200, 30))
        self.figure_bar_spin_y_min_2.setMaximumSize(QSize(200, 16777215))
        self.figure_bar_spin_y_min_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")
        self.figure_bar_spin_y_min_2.setMinimum(-1000000000.000000000000000)
        self.figure_bar_spin_y_min_2.setMaximum(1000000000.000000000000000)

        self.gridLayout_11.addWidget(self.figure_bar_spin_y_min_2, 4, 1, 1, 1)

        self.figure_bar_spin_y_max_2 = QDoubleSpinBox(self.tab_5)
        self.figure_bar_spin_y_max_2.setObjectName(u"figure_bar_spin_y_max_2")
        self.figure_bar_spin_y_max_2.setMinimumSize(QSize(200, 30))
        self.figure_bar_spin_y_max_2.setMaximumSize(QSize(200, 16777215))
        self.figure_bar_spin_y_max_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")
        self.figure_bar_spin_y_max_2.setMinimum(-1000000000.000000000000000)
        self.figure_bar_spin_y_max_2.setMaximum(1000000000.000000000000000)

        self.gridLayout_11.addWidget(self.figure_bar_spin_y_max_2, 4, 2, 1, 1)

        self.tabWidget.addTab(self.tab_5, "")
        self.tab_6 = QWidget()
        self.tab_6.setObjectName(u"tab_6")
        self.gridLayout_12 = QGridLayout(self.tab_6)
        self.gridLayout_12.setObjectName(u"gridLayout_12")
        self.figure_scatter_label_shape_size = QLabel(self.tab_6)
        self.figure_scatter_label_shape_size.setObjectName(u"figure_scatter_label_shape_size")
        self.figure_scatter_label_shape_size.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_12.addWidget(self.figure_scatter_label_shape_size, 4, 0, 1, 1)

        self.figure_scatter_label_legend = QLabel(self.tab_6)
        self.figure_scatter_label_legend.setObjectName(u"figure_scatter_label_legend")
        self.figure_scatter_label_legend.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_12.addWidget(self.figure_scatter_label_legend, 0, 0, 1, 1)

        self.figure_scatter_label_range = QLabel(self.tab_6)
        self.figure_scatter_label_range.setObjectName(u"figure_scatter_label_range")
        self.figure_scatter_label_range.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_12.addWidget(self.figure_scatter_label_range, 1, 0, 1, 1)

        self.figure_scatter_label_color = QLabel(self.tab_6)
        self.figure_scatter_label_color.setObjectName(u"figure_scatter_label_color")
        self.figure_scatter_label_color.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_12.addWidget(self.figure_scatter_label_color, 2, 0, 1, 1)

        self.figure_scatter_label_shape = QLabel(self.tab_6)
        self.figure_scatter_label_shape.setObjectName(u"figure_scatter_label_shape")
        self.figure_scatter_label_shape.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_12.addWidget(self.figure_scatter_label_shape, 3, 0, 1, 1)

        self.figure_scatter_spin_grid_size_2 = QDoubleSpinBox(self.tab_6)
        self.figure_scatter_spin_grid_size_2.setObjectName(u"figure_scatter_spin_grid_size_2")
        self.figure_scatter_spin_grid_size_2.setMinimumSize(QSize(0, 30))
        self.figure_scatter_spin_grid_size_2.setMaximumSize(QSize(200, 16777215))
        self.figure_scatter_spin_grid_size_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")
        self.figure_scatter_spin_grid_size_2.setMaximum(100000.000000000000000)
        self.figure_scatter_spin_grid_size_2.setValue(10.000000000000000)

        self.gridLayout_12.addWidget(self.figure_scatter_spin_grid_size_2, 0, 1, 1, 1)

        self.figure_scatter_color_min_2 = QDoubleSpinBox(self.tab_6)
        self.figure_scatter_color_min_2.setObjectName(u"figure_scatter_color_min_2")
        self.figure_scatter_color_min_2.setMinimumSize(QSize(0, 30))
        self.figure_scatter_color_min_2.setMaximumSize(QSize(200, 16777215))
        self.figure_scatter_color_min_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")

        self.gridLayout_12.addWidget(self.figure_scatter_color_min_2, 1, 1, 1, 1)

        self.figure_scatter_color_max_2 = QDoubleSpinBox(self.tab_6)
        self.figure_scatter_color_max_2.setObjectName(u"figure_scatter_color_max_2")
        self.figure_scatter_color_max_2.setMinimumSize(QSize(0, 30))
        self.figure_scatter_color_max_2.setMaximumSize(QSize(200, 16777215))
        self.figure_scatter_color_max_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")

        self.gridLayout_12.addWidget(self.figure_scatter_color_max_2, 1, 2, 1, 1)

        self.figure_scatter_spin_shape_size_2 = QDoubleSpinBox(self.tab_6)
        self.figure_scatter_spin_shape_size_2.setObjectName(u"figure_scatter_spin_shape_size_2")
        self.figure_scatter_spin_shape_size_2.setMinimumSize(QSize(0, 30))
        self.figure_scatter_spin_shape_size_2.setStyleSheet(u"QDoubleSpinBox {\n"
"    font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"    color: white;\n"
"    border: 1px solid white;\n"
"background-color:rgb(33, 37, 43);\n"
"}\n"
"\n"
"QDoubleSpinBox::up-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: top right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}\n"
"\n"
"QDoubleSpinBox::down-button {\n"
"    height: 10px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u9ad8\u5ea6 */\n"
"    width: 20px; /* \u8bbe\u7f6e\u4e00\u4e2a\u5408\u9002\u7684\u5bbd\u5ea6 */\n"
"    subcontrol-position: bottom right; /* \u4fdd\u6301\u6309\u94ae\u5728\u53f3\u4fa7 */\n"
"}")
        self.figure_scatter_spin_shape_size_2.setMaximum(100000.000000000000000)
        self.figure_scatter_spin_shape_size_2.setValue(10.000000000000000)

        self.gridLayout_12.addWidget(self.figure_scatter_spin_shape_size_2, 4, 1, 1, 1)

        self.figure_scatter_como_color_2 = QComboBox(self.tab_6)
        self.figure_scatter_como_color_2.addItem("")
        self.figure_scatter_como_color_2.addItem("")
        self.figure_scatter_como_color_2.addItem("")
        self.figure_scatter_como_color_2.addItem("")
        self.figure_scatter_como_color_2.addItem("")
        self.figure_scatter_como_color_2.addItem("")
        self.figure_scatter_como_color_2.addItem("")
        self.figure_scatter_como_color_2.addItem("")
        self.figure_scatter_como_color_2.setObjectName(u"figure_scatter_como_color_2")
        self.figure_scatter_como_color_2.setMinimumSize(QSize(0, 30))
        self.figure_scatter_como_color_2.setMaximumSize(QSize(200, 16777215))
        self.figure_scatter_como_color_2.setStyleSheet(u"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"background-color:rgb(33, 37, 43);")

        self.gridLayout_12.addWidget(self.figure_scatter_como_color_2, 2, 1, 1, 1)

        self.figure_scatter_btn_shape_2 = QPushButton(self.tab_6)
        self.figure_scatter_btn_shape_2.setObjectName(u"figure_scatter_btn_shape_2")
        self.figure_scatter_btn_shape_2.setMinimumSize(QSize(0, 30))
        self.figure_scatter_btn_shape_2.setStyleSheet(u"background-color: rgb(189,147,249);\n"
"color:white;\n"
"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_12.addWidget(self.figure_scatter_btn_shape_2, 3, 1, 1, 1)

        self.tabWidget.addTab(self.tab_6, "")

        self.verticalLayout_11.addWidget(self.tabWidget)

        self.btn_figure_run = QPushButton(self.widget_2)
        self.btn_figure_run.setObjectName(u"btn_figure_run")
        self.btn_figure_run.setMinimumSize(QSize(150, 30))
        self.btn_figure_run.setFont(font2)
        self.btn_figure_run.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.btn_figure_run.setLayoutDirection(Qt.LeftToRight)
        self.btn_figure_run.setStyleSheet(u"background-color: rgb(189,147,249);\n"
"color:white;\n"
"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"")

        self.verticalLayout_11.addWidget(self.btn_figure_run)


        self.horizontalLayout_8.addWidget(self.widget_2)

        self.figure_color_extra_box = QFrame(self.page_figure)
        self.figure_color_extra_box.setObjectName(u"figure_color_extra_box")
        self.figure_color_extra_box.setMaximumSize(QSize(0, 16777215))
        self.figure_color_extra_box.setStyleSheet(u"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")
        self.figure_color_extra_box.setFrameShape(QFrame.StyledPanel)
        self.figure_color_extra_box.setFrameShadow(QFrame.Raised)

        self.horizontalLayout_8.addWidget(self.figure_color_extra_box)

        self.figure_shape_extra_box = QFrame(self.page_figure)
        self.figure_shape_extra_box.setObjectName(u"figure_shape_extra_box")
        self.figure_shape_extra_box.setEnabled(True)
        self.figure_shape_extra_box.setMaximumSize(QSize(0, 16777215))
        self.figure_shape_extra_box.setStyleSheet(u"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";")
        self.figure_shape_extra_box.setFrameShape(QFrame.StyledPanel)
        self.figure_shape_extra_box.setFrameShadow(QFrame.Raised)

        self.horizontalLayout_8.addWidget(self.figure_shape_extra_box)

        self.stackedWidget.addWidget(self.page_figure)
        self.page_vmd = QWidget()
        self.page_vmd.setObjectName(u"page_vmd")
        self.verticalLayout_12 = QVBoxLayout(self.page_vmd)
        self.verticalLayout_12.setObjectName(u"verticalLayout_12")
        self.frame_6 = QFrame(self.page_vmd)
        self.frame_6.setObjectName(u"frame_6")
        self.frame_6.setFrameShape(QFrame.StyledPanel)
        self.frame_6.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_10 = QHBoxLayout(self.frame_6)
        self.horizontalLayout_10.setObjectName(u"horizontalLayout_10")
        self.vmd_btn_start = QPushButton(self.frame_6)
        self.vmd_btn_start.setObjectName(u"vmd_btn_start")
        self.vmd_btn_start.setStyleSheet(u"background-color: rgb(189,147,249);\n"
"color:white;\n"
"font-size:18pt;")

        self.horizontalLayout_10.addWidget(self.vmd_btn_start)

        self.vmd_btn_stop = QPushButton(self.frame_6)
        self.vmd_btn_stop.setObjectName(u"vmd_btn_stop")
        self.vmd_btn_stop.setStyleSheet(u"background-color: rgb(189,147,249);\n"
"color:white;\n"
"font-size:18pt;")

        self.horizontalLayout_10.addWidget(self.vmd_btn_stop)


        self.verticalLayout_12.addWidget(self.frame_6)

        self.vmd_tablewidget = QTableWidget(self.page_vmd)
        self.vmd_tablewidget.setObjectName(u"vmd_tablewidget")

        self.verticalLayout_12.addWidget(self.vmd_tablewidget)

        self.vmd_label = QLabel(self.page_vmd)
        self.vmd_label.setObjectName(u"vmd_label")

        self.verticalLayout_12.addWidget(self.vmd_label)

        self.stackedWidget.addWidget(self.page_vmd)

        self.verticalLayout_14.addWidget(self.stackedWidget)


        self.horizontalLayout_4.addWidget(self.pagesContainer)

        self.extraRightBox_3 = QFrame(self.content)
        self.extraRightBox_3.setObjectName(u"extraRightBox_3")
        self.extraRightBox_3.setMinimumSize(QSize(0, 0))
        self.extraRightBox_3.setMaximumSize(QSize(0, 16777215))
        self.extraRightBox_3.setFrameShape(QFrame.NoFrame)
        self.extraRightBox_3.setFrameShadow(QFrame.Raised)
        self.verticalLayout_7 = QVBoxLayout(self.extraRightBox_3)
        self.verticalLayout_7.setSpacing(0)
        self.verticalLayout_7.setObjectName(u"verticalLayout_7")
        self.verticalLayout_7.setContentsMargins(0, 0, 0, 0)

        self.horizontalLayout_4.addWidget(self.extraRightBox_3)


        self.verticalLayout_6.addWidget(self.content)

        self.bottomBar = QFrame(self.contentBottom)
        self.bottomBar.setObjectName(u"bottomBar")
        self.bottomBar.setMinimumSize(QSize(0, 22))
        self.bottomBar.setMaximumSize(QSize(16777215, 22))
        self.bottomBar.setFrameShape(QFrame.NoFrame)
        self.bottomBar.setFrameShadow(QFrame.Raised)
        self.horizontalLayout_5 = QHBoxLayout(self.bottomBar)
        self.horizontalLayout_5.setSpacing(0)
        self.horizontalLayout_5.setObjectName(u"horizontalLayout_5")
        self.horizontalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.label_2 = QLabel(self.bottomBar)
        self.label_2.setObjectName(u"label_2")
        self.label_2.setStyleSheet(u"font: 12pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.horizontalLayout_5.addWidget(self.label_2)

        self.version = QLabel(self.bottomBar)
        self.version.setObjectName(u"version")
        self.version.setStyleSheet(u"font: 12pt \"\u534e\u6587\u7ec6\u9ed1\";")
        self.version.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)

        self.horizontalLayout_5.addWidget(self.version)


        self.verticalLayout_6.addWidget(self.bottomBar)


        self.verticalLayout_2.addWidget(self.contentBottom)


        self.appLayout.addWidget(self.contentBox)


        self.verticalLayout_23.addWidget(self.bgApp)

        MainWindow.setCentralWidget(self.styleSheet)

        self.retranslateUi(MainWindow)

        self.stackedWidget.setCurrentIndex(4)
        self.tabWidget.setCurrentIndex(0)


        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"MainWindow", None))
        self.toggleButton.setText(QCoreApplication.translate("MainWindow", u"Hide", None))
        self.btn_home.setText(QCoreApplication.translate("MainWindow", u"Home", None))
        self.btn_generate.setText(QCoreApplication.translate("MainWindow", u"Generation", None))
        self.btn_figure.setText(QCoreApplication.translate("MainWindow", u"Figure", None))
        self.btn_analysis.setText(QCoreApplication.translate("MainWindow", u"Analysis", None))
        self.btn_data_process.setText(QCoreApplication.translate("MainWindow", u"VMD Integration", None))
        self.titleRightInfo.setText(QCoreApplication.translate("MainWindow", u"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'\u534e\u6587\u7ec6\u9ed1'; font-size:16pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:11pt;\">LNB-MDT</span></p></body></html>", None))
#if QT_CONFIG(tooltip)
        self.btn_language.setToolTip(QCoreApplication.translate("MainWindow", u"Language", None))
#endif // QT_CONFIG(tooltip)
        self.btn_language.setText(QCoreApplication.translate("MainWindow", u"EN", None))
#if QT_CONFIG(tooltip)
        self.minimizeAppBtn.setToolTip(QCoreApplication.translate("MainWindow", u"Minimize", None))
#endif // QT_CONFIG(tooltip)
        self.minimizeAppBtn.setText("")
#if QT_CONFIG(tooltip)
        self.maximizeRestoreAppBtn.setToolTip(QCoreApplication.translate("MainWindow", u"Maximize", None))
#endif // QT_CONFIG(tooltip)
        self.maximizeRestoreAppBtn.setText("")
#if QT_CONFIG(tooltip)
        self.closeAppBtn.setToolTip(QCoreApplication.translate("MainWindow", u"Close", None))
#endif // QT_CONFIG(tooltip)
        self.closeAppBtn.setText("")
        self.label_5.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p align=\"center\"><span style=\" font-size:72pt; font-weight:600;\">LNB-MDT</span></p></body></html>", None))
        self.btnTrajectory.setText(QCoreApplication.translate("MainWindow", u"Trajecory ", None))
        self.editStructure.setText("")
        self.editStructure.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Supported file formats: GRO PDB", None))
        self.editTrajectory.setText("")
        self.editTrajectory.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Supported file formats: XTC, TPR", None))
        self.btnSructure.setText(QCoreApplication.translate("MainWindow", u"Structure ", None))
        self.btnResult.setText(QCoreApplication.translate("MainWindow", u"PATH", None))
        self.editResult.setText("")
        self.editResult.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Please select the result save path", None))
        self.comboBoxMethod.setItemText(0, QCoreApplication.translate("MainWindow", u"Anisotropy", None))
        self.comboBoxMethod.setItemText(1, QCoreApplication.translate("MainWindow", u"Area", None))
        self.comboBoxMethod.setItemText(2, QCoreApplication.translate("MainWindow", u"Cluster", None))
        self.comboBoxMethod.setItemText(3, QCoreApplication.translate("MainWindow", u"Gyration", None))
        self.comboBoxMethod.setItemText(4, QCoreApplication.translate("MainWindow", u"Height", None))
        self.comboBoxMethod.setItemText(5, QCoreApplication.translate("MainWindow", u"MeanCurvature", None))
        self.comboBoxMethod.setItemText(6, QCoreApplication.translate("MainWindow", u"NCluster", None))
        self.comboBoxMethod.setItemText(7, QCoreApplication.translate("MainWindow", u"PCA", None))
        self.comboBoxMethod.setItemText(8, QCoreApplication.translate("MainWindow", u"Pressure", None))
        self.comboBoxMethod.setItemText(9, QCoreApplication.translate("MainWindow", u"RadialDistribution", None))
        self.comboBoxMethod.setItemText(10, QCoreApplication.translate("MainWindow", u"SZ", None))

        self.analysis_label_frame.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" font-size:18pt; font-weight:600;\">Frames</span></p></body></html>", None))
        self.analysis_label_method.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" font-size:16pt; font-weight:600;\">AnalysisMethods</span></p></body></html>", None))
        self.analysis_label_k.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" font-size:16pt; font-weight:600;\">K-Neighboors</span></p></body></html>", None))
        self.analysis_label_first.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" font-size:14pt; font-weight:600;\">First</span></p></body></html>", None))
        self.analysis_label_last.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" font-size:14pt; font-weight:600;\">Last</span></p></body></html>", None))
        self.analysis_label_step.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p><span style=\" font-size:14pt; font-weight:600;\">Step</span></p></body></html>", None))
        self.btnNext.setText("")
        self.label_gene_label.setText(QCoreApplication.translate("MainWindow", u"Generation", None))
        self.edit_gene_path.setText("")
        self.edit_gene_path.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Please select path to save gro and topol file", None))
        self.btn_gene_path.setText(QCoreApplication.translate("MainWindow", u"Save path", None))
        self.label_gene_y.setText(QCoreApplication.translate("MainWindow", u"Box Y Length(nm)", None))
        self.label_gene_x.setText(QCoreApplication.translate("MainWindow", u"Box X Length(nm)", None))
        self.label_gene_z.setText(QCoreApplication.translate("MainWindow", u"Box Z Length(nm)", None))
        self.label_gene_gas.setText(QCoreApplication.translate("MainWindow", u"Gas type", None))
        self.como_gas.setItemText(0, QCoreApplication.translate("MainWindow", u"N2", None))

        self.label_gene_gden.setText(QCoreApplication.translate("MainWindow", u"Gas density", None))
        self.label_gene_a.setText(QCoreApplication.translate("MainWindow", u"Area per lipid", None))
        self.label_gene_r.setText(QCoreApplication.translate("MainWindow", u"Radius of LNB", None))
        self.btn_gene_lipid.setText(QCoreApplication.translate("MainWindow", u"Select Lipid", None))
        self.label_gene_sol.setText(QCoreApplication.translate("MainWindow", u"Solvent Type", None))
        self.label_gene_salt.setText(QCoreApplication.translate("MainWindow", u"Salt concentration", None))
        self.como_solvent.setItemText(0, QCoreApplication.translate("MainWindow", u"W", None))

        self.btn_gene_run.setText(QCoreApplication.translate("MainWindow", u"RUN!", None))
        self.figure_edit_path.setInputMask("")
        self.figure_edit_path.setText("")
        self.figure_edit_path.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Please import Excel file", None))
        self.figure_btn_path.setText(QCoreApplication.translate("MainWindow", u"Select Excel File", None))
        self.figure_line_label_legend.setText(QCoreApplication.translate("MainWindow", u"Legend Size(0=None)", None))
        self.figure_line_btn_color_2.setText(QCoreApplication.translate("MainWindow", u"Select Color", None))
        self.figure_line_label_x.setText(QCoreApplication.translate("MainWindow", u"X-Title", None))
        self.figure_line_label_color.setText(QCoreApplication.translate("MainWindow", u"Color", None))
        self.figure_line_label_y_range.setText(QCoreApplication.translate("MainWindow", u"Y-Range", None))
        self.figure_line_label_x_range.setText(QCoreApplication.translate("MainWindow", u"X-Range", None))
        self.figure_line_label_y.setText(QCoreApplication.translate("MainWindow", u"Y-Title", None))
        self.figure_line_label_axis_tick.setText(QCoreApplication.translate("MainWindow", u"Axis Tick Size", None))
        self.figure_line_label_marker.setText(QCoreApplication.translate("MainWindow", u"Marker Size", None))
        self.figure_line_label_axis_title.setText(QCoreApplication.translate("MainWindow", u"Axis Title Size", None))
        self.figure_line_edit_y_2.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Default if none", None))
        self.figure_line_como_marker.setItemText(0, QCoreApplication.translate("MainWindow", u"o", None))
        self.figure_line_como_marker.setItemText(1, QCoreApplication.translate("MainWindow", u"p", None))
        self.figure_line_como_marker.setItemText(2, QCoreApplication.translate("MainWindow", u"s", None))
        self.figure_line_como_marker.setItemText(3, QCoreApplication.translate("MainWindow", u"^", None))
        self.figure_line_como_marker.setItemText(4, QCoreApplication.translate("MainWindow", u"*", None))
        self.figure_line_como_marker.setItemText(5, QCoreApplication.translate("MainWindow", u"x", None))
        self.figure_line_como_marker.setItemText(6, QCoreApplication.translate("MainWindow", u"+", None))

        self.figure_line_edit_x_2.setText("")
        self.figure_line_edit_x_2.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Default if none", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_4), QCoreApplication.translate("MainWindow", u"Line", None))
        self.figure_bar_btn_color_2.setText(QCoreApplication.translate("MainWindow", u"Select Color", None))
        self.figure_bar_btn_trend_2.setText(QCoreApplication.translate("MainWindow", u"Select Color", None))
        self.figure_bar_label_trend.setText(QCoreApplication.translate("MainWindow", u"Trend Line", None))
        self.figure_bar_label_y.setText(QCoreApplication.translate("MainWindow", u"Y-Title", None))
        self.figure_bar_edit_x_2.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Default if none", None))
        self.figure_bar_edit_y_2.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Default if none", None))
        self.figure_bar_label_axis_title.setText(QCoreApplication.translate("MainWindow", u"Axis Title Size", None))
        self.figure_bar_label_error.setText(QCoreApplication.translate("MainWindow", u"Error Bar", None))
        self.figure_bar_label_x.setText(QCoreApplication.translate("MainWindow", u"X-Title", None))
        self.figure_bar_label_axis_tick.setText(QCoreApplication.translate("MainWindow", u"Axis Tick Size", None))
        self.figure_bar_label_color.setText(QCoreApplication.translate("MainWindow", u"Color", None))
        self.figure_bar_label_bar.setText(QCoreApplication.translate("MainWindow", u"Bar Value", None))
        self.figure_bar_label_y_range.setText(QCoreApplication.translate("MainWindow", u"Y-Range", None))
        self.figure_bar_radio_error_2.setText(QCoreApplication.translate("MainWindow", u"Yes", None))
        self.radioButton_6.setText(QCoreApplication.translate("MainWindow", u"No", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_5), QCoreApplication.translate("MainWindow", u"Bar", None))
        self.figure_scatter_label_shape_size.setText(QCoreApplication.translate("MainWindow", u"Shape size", None))
        self.figure_scatter_label_legend.setText(QCoreApplication.translate("MainWindow", u"Legend Size", None))
        self.figure_scatter_label_range.setText(QCoreApplication.translate("MainWindow", u"Value Range", None))
        self.figure_scatter_label_color.setText(QCoreApplication.translate("MainWindow", u"Color Type", None))
        self.figure_scatter_label_shape.setText(QCoreApplication.translate("MainWindow", u"Shape", None))
        self.figure_scatter_como_color_2.setItemText(0, QCoreApplication.translate("MainWindow", u"bwr", None))
        self.figure_scatter_como_color_2.setItemText(1, QCoreApplication.translate("MainWindow", u"binary", None))
        self.figure_scatter_como_color_2.setItemText(2, QCoreApplication.translate("MainWindow", u"Blues", None))
        self.figure_scatter_como_color_2.setItemText(3, QCoreApplication.translate("MainWindow", u"Greys", None))
        self.figure_scatter_como_color_2.setItemText(4, QCoreApplication.translate("MainWindow", u"inferno", None))
        self.figure_scatter_como_color_2.setItemText(5, QCoreApplication.translate("MainWindow", u"Oranges", None))
        self.figure_scatter_como_color_2.setItemText(6, QCoreApplication.translate("MainWindow", u"RaBu", None))
        self.figure_scatter_como_color_2.setItemText(7, QCoreApplication.translate("MainWindow", u"viridis", None))

        self.figure_scatter_btn_shape_2.setText(QCoreApplication.translate("MainWindow", u"Select Shape", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_6), QCoreApplication.translate("MainWindow", u"Scatter", None))
        self.btn_figure_run.setText(QCoreApplication.translate("MainWindow", u"RUN!", None))
        self.vmd_btn_start.setText(QCoreApplication.translate("MainWindow", u"Start VMD", None))
        self.vmd_btn_stop.setText(QCoreApplication.translate("MainWindow", u"Stop VMD", None))
        self.vmd_label.setText(QCoreApplication.translate("MainWindow", u"TextLabel", None))
        self.label_2.setText(QCoreApplication.translate("MainWindow", u"By\uff1aNanoBioMembrane Lab", None))
        self.version.setText(QCoreApplication.translate("MainWindow", u"<html><head/><body><p>v1.0</p></body></html>", None))
    # retranslateUi


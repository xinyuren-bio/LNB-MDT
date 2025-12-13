# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'Analysis_Figure.ui'
##
## Created by: Qt User Interface Compiler version 6.7.1
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
from PySide6.QtWidgets import (QApplication, QComboBox, QDoubleSpinBox, QFrame,
    QGridLayout, QHBoxLayout, QLabel, QLineEdit,
    QPushButton, QRadioButton, QSizePolicy, QTabWidget,
    QVBoxLayout, QWidget)

class Ui_Form(object):
    def setupUi(self, Form):
        if not Form.objectName():
            Form.setObjectName(u"Form")
        Form.resize(625, 405)
        self.horizontalLayout = QHBoxLayout(Form)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.frame = QFrame(Form)
        self.frame.setObjectName(u"frame")
        self.frame.setFrameShape(QFrame.StyledPanel)
        self.frame.setFrameShadow(QFrame.Raised)
        self.verticalLayout = QVBoxLayout(self.frame)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.tabWidget = QTabWidget(self.frame)
        self.tabWidget.setObjectName(u"tabWidget")
        self.tabWidget.setCursor(QCursor(Qt.PointingHandCursor))
        self.tabWidget.setStyleSheet(u"QTabBar::tab {\n"
"    background: lightgray;\n"
"    border: 2px solid #C4C4C3;\n"
"    border-bottom-color: #C4C4C3; /* same as the pane color */\n"
"    border-top-left-radius: 5px;\n"
"    border-top-right-radius: 5px;\n"
"    min-width: 16ex;\n"
"    padding: 2px;\n"
"	font: 18pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"	color:balck;\n"
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
        self.figure_line_btn_color_2.setCursor(QCursor(Qt.PointingHandCursor))
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
"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"color:white;\n"
"")

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
"font: 14pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"font-weight:bold;\n"
"border: 1px solid white;\n"
"color:white;")

        self.gridLayout_9.addWidget(self.figure_line_como_marker, 7, 2, 1, 1)

        self.figure_line_edit_x_2 = QLineEdit(self.tab_4)
        self.figure_line_edit_x_2.setObjectName(u"figure_line_edit_x_2")
        self.figure_line_edit_x_2.setMinimumSize(QSize(0, 30))
        self.figure_line_edit_x_2.setMaximumSize(QSize(16777215, 16777215))
        self.figure_line_edit_x_2.setStyleSheet(u"background-color: rgb(33, 37, 43);\n"
"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"color:white;")

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
"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"color:white;")

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
"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
"\n"
"color:white;")

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
"background-color:rgb(33, 37, 43);\n"
"\n"
"color:white;")

        self.gridLayout_12.addWidget(self.figure_scatter_como_color_2, 2, 1, 1, 1)

        self.figure_scatter_btn_shape_2 = QPushButton(self.tab_6)
        self.figure_scatter_btn_shape_2.setObjectName(u"figure_scatter_btn_shape_2")
        self.figure_scatter_btn_shape_2.setMinimumSize(QSize(0, 30))
        self.figure_scatter_btn_shape_2.setStyleSheet(u"background-color: rgb(189,147,249);\n"
"color:white;\n"
"font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";")

        self.gridLayout_12.addWidget(self.figure_scatter_btn_shape_2, 3, 1, 1, 1)

        self.tabWidget.addTab(self.tab_6, "")

        self.verticalLayout.addWidget(self.tabWidget)

        self.btn_figure_run = QPushButton(self.frame)
        self.btn_figure_run.setObjectName(u"btn_figure_run")
        self.btn_figure_run.setMinimumSize(QSize(150, 30))
        font = QFont()
        font.setFamilies([u"\u534e\u6587\u7ec6\u9ed1"])
        font.setPointSize(16)
        font.setBold(False)
        font.setItalic(False)
        self.btn_figure_run.setFont(font)
        self.btn_figure_run.setCursor(QCursor(Qt.PointingHandCursor))
        self.btn_figure_run.setLayoutDirection(Qt.LeftToRight)
        # 内联样式已移除，使用外部样式表统一控制按钮样式
        # self.btn_figure_run.setStyleSheet(u"background-color: rgb(189,147,249);\n"
# "color:white;\n"
# "font: 16pt \"\u534e\u6587\u7ec6\u9ed1\";\n"
# "")

        self.verticalLayout.addWidget(self.btn_figure_run)


        self.horizontalLayout.addWidget(self.frame)

        self.figure_color_extra_box = QFrame(Form)
        self.figure_color_extra_box.setObjectName(u"figure_color_extra_box")
        self.figure_color_extra_box.setMaximumSize(QSize(0, 16777215))
        self.figure_color_extra_box.setFrameShape(QFrame.StyledPanel)
        self.figure_color_extra_box.setFrameShadow(QFrame.Raised)

        self.horizontalLayout.addWidget(self.figure_color_extra_box)

        self.figure_shape_extra_box = QFrame(Form)
        self.figure_shape_extra_box.setObjectName(u"figure_shape_extra_box")
        self.figure_shape_extra_box.setMaximumSize(QSize(0, 16777215))
        self.figure_shape_extra_box.setFrameShape(QFrame.StyledPanel)
        self.figure_shape_extra_box.setFrameShadow(QFrame.Raised)

        self.horizontalLayout.addWidget(self.figure_shape_extra_box)


        self.retranslateUi(Form)

        self.tabWidget.setCurrentIndex(0)


        QMetaObject.connectSlotsByName(Form)
    # setupUi

    def retranslateUi(self, Form):
        Form.setWindowTitle(QCoreApplication.translate("Form", u"Form", None))
        self.figure_line_label_legend.setText(QCoreApplication.translate("Form", u"Legend Size(0=None)", None))
        self.figure_line_btn_color_2.setText(QCoreApplication.translate("Form", u"Select Color", None))
        self.figure_line_label_x.setText(QCoreApplication.translate("Form", u"X-Title", None))
        self.figure_line_label_color.setText(QCoreApplication.translate("Form", u"Color", None))
        self.figure_line_label_y_range.setText(QCoreApplication.translate("Form", u"Y-Range", None))
        self.figure_line_label_x_range.setText(QCoreApplication.translate("Form", u"X-Range", None))
        self.figure_line_label_y.setText(QCoreApplication.translate("Form", u"Y-Title", None))
        self.figure_line_label_axis_tick.setText(QCoreApplication.translate("Form", u"Axis Tick Size", None))
        self.figure_line_label_marker.setText(QCoreApplication.translate("Form", u"Marker Size", None))
        self.figure_line_label_axis_title.setText(QCoreApplication.translate("Form", u"Axis Title Size", None))
        self.figure_line_edit_y_2.setPlaceholderText(QCoreApplication.translate("Form", u"Default if none", None))
        self.figure_line_como_marker.setItemText(0, QCoreApplication.translate("Form", u"o", None))
        self.figure_line_como_marker.setItemText(1, QCoreApplication.translate("Form", u"p", None))
        self.figure_line_como_marker.setItemText(2, QCoreApplication.translate("Form", u"s", None))
        self.figure_line_como_marker.setItemText(3, QCoreApplication.translate("Form", u"^", None))
        self.figure_line_como_marker.setItemText(4, QCoreApplication.translate("Form", u"*", None))
        self.figure_line_como_marker.setItemText(5, QCoreApplication.translate("Form", u"x", None))
        self.figure_line_como_marker.setItemText(6, QCoreApplication.translate("Form", u"+", None))

        self.figure_line_edit_x_2.setText("")
        self.figure_line_edit_x_2.setPlaceholderText(QCoreApplication.translate("Form", u"Default if none", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_4), QCoreApplication.translate("Form", u"Line", None))
        self.figure_bar_btn_color_2.setText(QCoreApplication.translate("Form", u"Select Color", None))
        self.figure_bar_btn_trend_2.setText(QCoreApplication.translate("Form", u"Select Color", None))
        self.figure_bar_label_trend.setText(QCoreApplication.translate("Form", u"Trend Line", None))
        self.figure_bar_label_y.setText(QCoreApplication.translate("Form", u"Y-Title", None))
        self.figure_bar_edit_x_2.setPlaceholderText(QCoreApplication.translate("Form", u"Default if none", None))
        self.figure_bar_edit_y_2.setPlaceholderText(QCoreApplication.translate("Form", u"Default if none", None))
        self.figure_bar_label_axis_title.setText(QCoreApplication.translate("Form", u"Axis Title Size", None))
        self.figure_bar_label_error.setText(QCoreApplication.translate("Form", u"Error Bar", None))
        self.figure_bar_label_x.setText(QCoreApplication.translate("Form", u"X-Title", None))
        self.figure_bar_label_axis_tick.setText(QCoreApplication.translate("Form", u"Axis Tick Size", None))
        self.figure_bar_label_color.setText(QCoreApplication.translate("Form", u"Color", None))
        self.figure_bar_label_bar.setText(QCoreApplication.translate("Form", u"Bar Value", None))
        self.figure_bar_label_y_range.setText(QCoreApplication.translate("Form", u"Y-Range", None))
        self.figure_bar_radio_error_2.setText(QCoreApplication.translate("Form", u"Yes", None))
        self.radioButton_6.setText(QCoreApplication.translate("Form", u"No", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_5), QCoreApplication.translate("Form", u"Bar", None))
        self.figure_scatter_label_shape_size.setText(QCoreApplication.translate("Form", u"Shape size", None))
        self.figure_scatter_label_legend.setText(QCoreApplication.translate("Form", u"Legend Size", None))
        self.figure_scatter_label_range.setText(QCoreApplication.translate("Form", u"Value Range", None))
        self.figure_scatter_label_color.setText(QCoreApplication.translate("Form", u"Color Type", None))
        self.figure_scatter_label_shape.setText(QCoreApplication.translate("Form", u"Shape", None))
        self.figure_scatter_como_color_2.setItemText(0, QCoreApplication.translate("Form", u"bwr", None))
        self.figure_scatter_como_color_2.setItemText(1, QCoreApplication.translate("Form", u"binary", None))
        self.figure_scatter_como_color_2.setItemText(2, QCoreApplication.translate("Form", u"Blues", None))
        self.figure_scatter_como_color_2.setItemText(3, QCoreApplication.translate("Form", u"Greys", None))
        self.figure_scatter_como_color_2.setItemText(4, QCoreApplication.translate("Form", u"inferno", None))
        self.figure_scatter_como_color_2.setItemText(5, QCoreApplication.translate("Form", u"Oranges", None))
        self.figure_scatter_como_color_2.setItemText(6, QCoreApplication.translate("Form", u"RaBu", None))
        self.figure_scatter_como_color_2.setItemText(7, QCoreApplication.translate("Form", u"viridis", None))

        self.figure_scatter_btn_shape_2.setText(QCoreApplication.translate("Form", u"Select Shape", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_6), QCoreApplication.translate("Form", u"Scatter", None))
        self.btn_figure_run.setText(QCoreApplication.translate("Form", u"RUN!", None))
    # retranslateUi


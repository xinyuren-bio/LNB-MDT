from collections import defaultdict
from functools import partial

from PySide6.QtGui import QColor

from figure.figure import *
from .Tools import *
# from exception import *

global window


TYPE_ID = {'Height(nm)': 0
               , 'SZ': 0
               , 'Mean Curvature(nm-1)': 0
               , 'Area(nm^2)': 0
               , 'Anisotropy': 1
               , 'Gyration(nm)': 1
               , 'Cluster': 1
               , 'PCA': 1
               , 'RadialDistribution': 2
               , 'merge': 2
               , 'Bubble Height(nm)': 1
               , 'Bubble SZ': 1
               , 'Bubble Mean Curvature(nm-1)': 1
               , 'Bubble Area(nm^2)': 1
           }

FIGURE_TYPE = {0: ['Bar', 'Line', 'Scatter', 'Map']  # Lipids
    , 1: ['Bar', 'Line']  # bubble
    , 2: ['Line']}  # other


class FigureGetInfo:
    """
    获取FigureWidget内部的信息
    """
    FIGURE_METHOD: [str] = ['Line', 'Bar', 'Scatter']

    def __init__(self, ui):
        self.ui = ui
        # path
        try:
            self.path_figure = (self.ui.
                                figure_edit_path.text())
            self.figureMethod = self.FIGURE_METHOD[self.ui.tabWidget.currentIndex()]  # 0:line 1:bar 2:scatter

            self.LineInfo = defaultdict(lambda: None)
            self.BarInfo = defaultdict(lambda: None)
            self.ScaInfo = defaultdict(lambda: None)

            self.ColorInfo = defaultdict(lambda: None)
            self.ShapeInfo = []

            self.description, self.results = read_excel(self.path_figure)
            self.lipids_type = self.results['Resname'].unique() if 'Resname' in self.results.columns else None
            self.residues = None
        except:
            self.ui.figure_extra = 0

    @classmethod
    def _get_xlsx_path(cls, ui) -> str:
        return ui.figure_edit_path.text()

    def getLine(self):

        self.LineInfo.update({
            'axis_scale': self.ui.figure_line_spin_axis_scale_2.value()
            , 'axis_text': self.ui.figure_line_spin_axis_text_size_2.value()
            , 'grid_size': self.ui.figure_line_spin_grid_size_2.value()
            , 'x_title': self.ui.figure_line_edit_x_2.text() or "Frames"
            , 'y_title': self.ui.figure_line_edit_y_2.text() or self.description
            , 'x_min': self.ui.figure_line_spin_x_min_2.value()
            , 'x_max': self.ui.figure_line_spin_x_max_2.value()
            , 'y_min': self.ui.figure_line_spin_y_min_2.value()
            , 'y_max': self.ui.figure_line_spin_y_max_2.value()
            , 'marker_size': self.ui.figure_line_spin_marker_size_2.value()
            , 'marker_shape': self.ui.figure_line_como_marker.currentText()
            , 'color': self.ui.FigureInfo.ColorInfo
        })

        return self.LineInfo

    def getBar(self):

        self.BarInfo.update({
            'axis_scale': self.ui.figure_bar_spin_axis_scale_2.value()
            , 'axis_text': self.ui.figure_bar_spin_axis_text_size_2.value()
            , 'x_title': self.ui.figure_bar_edit_x_2.text()
            , 'y_title': self.ui.figure_bar_edit_y_2.text() or self.description
            , 'y_min': self.ui.figure_bar_spin_y_min_2.value()
            , 'y_max': self.ui.figure_bar_spin_y_max_2.value()
            , 'trend_size': self.ui.figure_bar_spin_trend_2.value()
            , 'up_bar_value': self.ui.figure_bar_spin_bar_2.value()
            , 'error_deci': self.ui.figure_bar_radio_error_2.isChecked()
            , 'color': self.ui.FigureInfo.ColorInfo
        })

        return self.BarInfo

    def getScatter(self):
        self.ScaInfo.update({
            'grid_size': self.ui.figure_scatter_spin_grid_size_2.value()
            , 'bar_min': self.ui.figure_scatter_color_min_2.value()
            , 'bar_max': self.ui.figure_scatter_color_max_2.value()
            , 'bar_color': self.ui.figure_scatter_como_color_2.currentText()
            , 'shape_size': self.ui.figure_scatter_spin_shape_size_2.value()
            , 'shape': {}
        })
        for group in self.ShapeInfo:
            self.ScaInfo['shape'][group.title()] = next(
                (radio.text() for radio in group.findChildren(QRadioButton) if radio.isChecked()), None
            )
        return self.ScaInfo


class FigurePage:


    @staticmethod
    def _ensure_figure_info(ui):
        """
        用来检测UI界面是否导入了结果的文件路径，以及结果文件路径有没有发生更改
        """
        if ui.FigureInfo is None or ui.FigureInfo.path_figure != ui.figure_edit_path.text():
                ui.FigureInfo = FigureGetInfo(ui)
        return ui.FigureInfo

    @classmethod
    def figureBtnMakeFigure(cls, ui):
        try:
            cls._ensure_figure_info(ui)
        except:
            create_warn_dialog("Please select the result file first.\nError in the function:figureBtnMakeFigure")
            return
        # 没有写入if的原因是，如果文件路径没有改变，则无需再次读取文件。
        ui.FigureInfo.figureMethod = FigureGetInfo.FIGURE_METHOD[ui.tabWidget.currentIndex()]  # 更新图形类型

        # 使用策略模式进行图表绘制类型的选择
        drawing_strategies = {
            'Line': lambda: FigureLine(ui.FigureInfo.description
                                       , ui.FigureInfo.results
                                       , ui.FigureInfo.getLine()).plot()
            , 'Bar': lambda: FigureBar(ui.FigureInfo.description
                                       , ui.FigureInfo.results
                                       , ui.FigureInfo.getBar()).plot()
            , 'Scatter': lambda: FigureScatter(ui.FigureInfo.description
                                               , ui.FigureInfo.results
                                               , ui.FigureInfo.getScatter()).plot()
        }
        draw_action = drawing_strategies[ui.FigureInfo.figureMethod]
        if draw_action:
            draw_action()

    @classmethod
    def figureBtnColor(cls, ui):
        # 确保 FigureInfo 已经初始化
        try:
            cls._ensure_figure_info(ui)
        except:
            create_warn_dialog("Please check your file type!\nError in the function:figureBtnColor")
            return
        color_strategies = {
            0: lambda: cls.lipids_colors(ui)
            , 1: lambda: cls.single_color(ui)
            , 2: lambda: cls.type_color(ui)
        }
        color_action = color_strategies[cls.TYPE_ID[ui.FigureInfo.description]]
        if color_action:
            color_action()

    @staticmethod
    def _ensure_color_lipids_sel(ui):
        if not getattr(ui, 'FigureColorLayout', None):
            ui.figure_color_extra_box.setStyleSheet(u"font: 15pt \"\u534e\u6587\u7ec6\u9ed1\";")
            ui.FigureColorLayout = QVBoxLayout()
            ui.figure_color_extra_box.setLayout(ui.FigureColorLayout)
        return None

    @staticmethod
    def _make_residues_btn(ui, callback):
        """用来创建选择Lipids颜色的时候的功能"""
        if ui.FigureInfo.path_figure != ui.figure_edit_path.text() or ui.FigureInfo.residues is None:
            residues = ui.FigureInfo.lipids_type
            ui.FigureInfo.residues = residues
            while ui.figure_color_extra_box.layout().count():
                item = ui.figure_color_extra_box.layout().takeAt(0)
                if item.widget():
                    item.widget().deleteLater()
                elif item.layout():
                    item.layout().deleteLater()
            for index, residue in enumerate(ui.FigureInfo.residues):
                btn = UIItemsMake.make_btn(
                    residue
                    , background_color='white'
                    , font_color='black'
                )
                btn.clicked.connect(partial(callback, index, btn, ui))
                ui.FigureColorLayout.addWidget(btn)

    @classmethod
    def lipids_colors(cls, ui):
        """
        类方法
        首先设定extra编号，用来决定是否显示右侧的框
        确保颜色布局建立
        确保读取了残基类别
        """
        ui.figure_extra = 1
        cls._ensure_color_lipids_sel(ui)  # 确保颜色选择布局已经创建

        cls._make_residues_btn(ui, cls.residues_color_sel)

        if ui.FigureInfo.lipids_type is None:
            raise ExcelFormatError("请严格遵守结果文件的格式！")
        # 得到残基的信息
        ui.figure_color_extra_box.show()

    @classmethod
    def single_color(cls, ui):
        """设置bubble颜色"""
        ui.figure_extra = 0
        if ui.FigureInfo.figureMethod == 'Line':
            ui.FigureInfo.LineInfo['bubble_color'] = cls.single_color_sel()
        elif ui.FigureInfo.figureMethod == 'Bar':
            ui.FigureInfo.BarInfo['bubble_color'] = cls.single_color_sel()

    @classmethod
    def type_color(cls, ui):
        """设置其他颜色，例如径向分布函数"""
        pass

    @classmethod
    def figureBtnShape(cls, ui):
        try:
            cls._ensure_figure_info(ui)
        except:
            create_warn_dialog("Please check your file type.\nError in the function:figureBtnShape")
            return
        residues_column = ui.FigureInfo.results.iloc[1:, 1]
        unique_residues = set(residues_column)
        ui.FigureInfo.residues = list(unique_residues)
        if not getattr(ui, 'FigureShapeWidget'):
            ui.FigureShapeLayout = QVBoxLayout(ui.figure_shape_extra_box)
            label_shape = QLabel('Shape')
            ui.FigureShapeLayout.addWidget(label_shape)
            scrollArea = QScrollArea()
            scrollArea.setWidgetResizable(True)
            ui.FigureShapeWidget = QWidget()
            containerLayout = QVBoxLayout()
            for sp in ui.FigureInfo.residues:
                groupBox = UIItemsMake.make_group_box(sp)
                groupLayout = QVBoxLayout(groupBox)
                for sh in ['o', 'p', 's', '^', '*', 'x', '+']:
                    radio = QRadioButton(sh)
                    groupLayout.addWidget(radio)
                ui.FigureInfo.ShapeInfo.append(groupBox)
                containerLayout.addWidget(groupBox)
            ui.FigureShapeWidget.setLayout(containerLayout)
            scrollArea.setWidget(ui.FigureShapeWidget)
            ui.FigureShapeLayout.addWidget(scrollArea)

    @classmethod
    def residues_color_sel(cls, index, btn, ui):
        new_color = cls._open_and_get_qcolor()
        ui.FigureInfo.ColorInfo[ui.FigureInfo.residues[index]] = new_color
        btn.setStyleSheet(f'background-color: rgb({new_color[0]*255}, {new_color[1]*255}, {new_color[2]*255});'
                          f'border-radius: {UISettings.BTN_BORDER_RADIUS};'
                          f'width: {UISettings.BTN_WIDTH};'
                          f'height: {UISettings.BTN_HEIGHT};')

    @classmethod
    def single_color_sel(cls):
        """用于选择趋势线的数值"""
        return cls._open_and_get_qcolor()

    @staticmethod
    def _open_and_get_qcolor():
        """用于打开PySide6的QColor按钮，并且返回所选择的color"""
        color = QColor()
        new_color = QColorDialog.getColor(color)
        if new_color.isValid():
            return (new_color.red() / 255
                    , new_color.green() / 255
                    , new_color.blue() / 255)
        return None
# ///////////////////////////////////////////////////////////////
#
# BY: RenXinYu
# PROJECT MADE WITH: Qt Designer and PySide6
# V: 1.0
#
#
#
# ///////////////////////////////////////////////////////////////
# MAIN FILE
# ///////////////////////////////////////////////////////////////
from main import *
from .Tools import *

# WITH ACCESS TO MAIN WINDOW WIDGETS
# ///////////////////////////////////////////////////////////////

global window


class AppFunctions(MainWindow):
    # def setThemeHack(self):
    #     Settings.BTN_LEFT_BOX_COLOR = "background-color: #495474;"
    #     Settings.BTN_RIGHT_BOX_COLOR = "background-color: #495474;"
    #     Settings.MENU_SELECTED_STYLESHEET = MENU_SELECTED_STYLESHEET = """
    #     border-left: 22px solid qlineargradient(spread:pad, x1:0.034, y1:0, x2:0.216, y2:0, stop:0.499 rgba(255, 121, 198, 255), stop:0.5 rgba(85, 170, 255, 0));
    #     background-color: #566388;
    #     """

        # SET MANUAL STYLES
        # self.lineEdit.setStyleSheet("background-color: #6272a4;font-size:15pt;")
        # self.ui.pushButton.setStyleSheet("background-color: #6272a4;")
        # self.ui.plainTextEdit.setStyleSheet("background-color: #6272a4;")
        # self.ui.tableWidget.setStyleSheet(
        #     "QScrollBar:vertical { background: #6272a4; } QScrollBar:horizontal { background: #6272a4; }")
        # self.ui.scrollArea.setStyleSheet(
        #     "QScrollBar:vertical { background: #6272a4; } QScrollBar:horizontal { background: #6272a4; }")
        # self.ui.comboBox.setStyleSheet("background-color: #6272a4;")
        # self.ui.horizontalScrollBar.setStyleSheet("background-color: #6272a4;")
        # self.ui.verticalScrollBar.setStyleSheet("background-color: #6272a4;")
        # self.ui.commandLinkButton.setStyleSheet("color: #ff79c6;")

    def btnLanguageClick(self):
        if self.btn_language.text() == "EN":
            font_style = "font: 16pt '华文细黑';"
            self.btn_language.setText('中')
            # Generation
            self.label_gene_label.setText('生成')
            self.label_gene_label.setStyleSheet(font_style)
            self.edit_gene_path.setPlaceholderText('请选择保存路径')
            self.btn_gene_path.setText('保存文件')
            self.label_gene_x.setText('盒子X长度(nm)')
            self.label_gene_x.setStyleSheet(font_style)
            self.label_gene_y.setText('盒子Y长度(nm)')
            self.label_gene_y.setStyleSheet(font_style)
            self.label_gene_z.setText('盒子Z长度(nm)')
            self.label_gene_z.setStyleSheet(font_style)
            self.label_gene_gas.setText('气体种类')
            self.label_gene_gas.setStyleSheet(font_style)
            self.label_gene_gden.setText('气体密度')
            self.label_gene_gden.setStyleSheet(font_style)
            self.label_gene_a.setText('单个磷脂面积')
            self.label_gene_a.setStyleSheet(font_style)
            self.label_gene_r.setText('脂质纳气泡半径')
            self.label_gene_r.setStyleSheet(font_style)
            self.label_gene_sol.setText('溶剂种类')
            self.label_gene_sol.setStyleSheet(font_style)
            self.label_gene_salt.setText('盐浓度')
            self.label_gene_salt.setStyleSheet(font_style)
            self.btn_gene_lipid.setText('磷脂选择')
            self.btn_gene_run.setText('开始生成！')

            # Figure
            self.btn_figure_run.setText('开始绘图！')
            self.figure_btn_path.setText('选择Excel文件')
            self.figure_edit_path.setPlaceholderText('请导入Excel文件')
            ## Line
            self.figure_line_label_axis_tick.setText('坐标刻度大小')
            self.figure_line_label_axis_tick.setStyleSheet(font_style)
            self.figure_line_label_axis_title.setText('标题大小')
            self.figure_line_label_axis_title.setStyleSheet(font_style)
            self.figure_line_label_legend.setText('图例大小')
            self.figure_line_label_legend.setStyleSheet(font_style)
            self.figure_line_label_x.setText('X轴标题')
            self.figure_line_label_x.setStyleSheet(font_style)
            self.figure_line_label_y.setText('Y轴标题')
            self.figure_line_label_y.setStyleSheet(font_style)
            self.figure_line_label_x_range.setText('X轴显示范围')
            self.figure_line_label_x_range.setStyleSheet(font_style)
            self.figure_line_label_y_range.setText('Y轴显示范围')
            self.figure_line_label_y_range.setStyleSheet(font_style)
            self.figure_line_label_marker.setText('标记点大小')
            self.figure_line_label_marker.setStyleSheet(font_style)
            self.figure_line_label_color.setText('颜色设置')
            self.figure_line_label_color.setStyleSheet(font_style)
            self.figure_line_btn_color_2.setText('选择颜色')
            ## Bar
            self.figure_bar_label_axis_tick.setText('坐标刻度大小')
            self.figure_bar_label_axis_tick.setStyleSheet(font_style)
            self.figure_bar_label_axis_title.setText('标题大小')
            self.figure_bar_label_axis_title.setStyleSheet(font_style)
            self.figure_bar_label_x.setText('X轴标题')
            self.figure_bar_label_x.setStyleSheet(font_style)
            self.figure_bar_label_y.setText('Y轴标题')
            self.figure_bar_label_y.setStyleSheet(font_style)
            self.figure_bar_label_y_range.setText('Y轴显示范围')
            self.figure_bar_label_y_range.setStyleSheet(font_style)
            self.figure_bar_label_trend.setText('趋势线设置')
            self.figure_bar_label_trend.setStyleSheet(font_style)
            self.figure_bar_label_bar.setText('柱图上数值设置')
            self.figure_bar_label_bar.setStyleSheet(font_style)
            self.figure_bar_label_color.setText('柱图颜色设置')
            self.figure_bar_label_color.setStyleSheet(font_style)
            self.figure_bar_label_error.setText('误差棒设置')
            self.figure_bar_label_error.setStyleSheet(font_style)
            self.figure_bar_btn_trend_2.setText('选择颜色')
            self.figure_bar_btn_color_2.setText('选择颜色')
            ## Scatter
            self.figure_scatter_label_legend.setText('图例大小')
            self.figure_scatter_label_legend.setStyleSheet(font_style)
            self.figure_scatter_label_range.setText('数值显示范围')
            self.figure_scatter_label_range.setStyleSheet(font_style)
            self.figure_scatter_label_color.setText('颜色类型设置')
            self.figure_scatter_label_color.setStyleSheet(font_style)
            self.figure_scatter_label_shape.setText('形状设置')
            self.figure_scatter_label_shape.setStyleSheet(font_style)
            self.figure_scatter_label_shape_size.setText('形状大小')
            self.figure_scatter_label_shape_size.setStyleSheet(font_style)

            # Analysis
            self.editStructure.setPlaceholderText('支持的结构文件格式：GRO PDB')
            self.editTrajectory.setPlaceholderText('支持的轨迹文件格式：XTC TPR')
            self.editResult.setPlaceholderText('请选择结果的保存路径')
            self.btnSructure.setText('结构文件')
            self.btnResult.setText('结果路径')
            self.btnTrajectory.setText('轨迹文件')
            self.analysis_label_frame.setText('帧数')
            self.analysis_label_frame.setStyleSheet(font_style)
            self.analysis_label_first.setText('起始帧')
            self.analysis_label_first.setStyleSheet(font_style)
            self.analysis_label_last.setText('结束帧')
            self.analysis_label_last.setStyleSheet(font_style)
            self.analysis_label_step.setText('间隔步数')
            self.analysis_label_step.setStyleSheet(font_style)
            self.analysis_label_k.setText('K值设置')
            self.analysis_label_k.setStyleSheet(font_style)
            self.analysis_label_method.setText('分析方法')
            self.analysis_label_method.setStyleSheet(font_style)

        elif self.btn_language.text() == "中":
            font_style = "font: 16pt '华文细黑';"
            self.btn_language.setText('EN')
            # Generation
            self.label_gene_label.setText('Generation')
            self.edit_gene_path.setPlaceholderText('Please select path to save gro and topol file')
            self.btn_gene_path.setText('Save Path')
            self.label_gene_x.setText('Box X Length(nm)')
            self.label_gene_x.setStyleSheet(font_style)
            self.label_gene_y.setText('Box Y Length(nm)')
            self.label_gene_y.setStyleSheet(font_style)
            self.label_gene_z.setText('Box Z Length(nm)')
            self.label_gene_z.setStyleSheet(font_style)
            self.label_gene_gas.setText('Gas Type')
            self.label_gene_gas.setStyleSheet(font_style)
            self.label_gene_gden.setText('Gas Density')
            self.label_gene_gden.setStyleSheet(font_style)
            self.label_gene_a.setText('Area per Lipid')
            self.label_gene_a.setStyleSheet(font_style)
            self.label_gene_r.setText('Radius of LNB')
            self.label_gene_r.setStyleSheet(font_style)
            self.label_gene_sol.setText('Solvent Type')
            self.label_gene_sol.setStyleSheet(font_style)
            self.label_gene_salt.setText('Salt Concentration')
            self.label_gene_salt.setStyleSheet(font_style)
            self.btn_gene_lipid.setText('Select Lipids')
            self.btn_gene_run.setText('Run！')
            # Figure
            self.btn_figure_run.setText('Run！')
            self.figure_btn_path.setText('Select Excel File')

            self.figure_edit_path.setPlaceholderText('Please import Excel file')

            ## Line
            self.figure_line_label_axis_tick.setText('Axis Tick Size')
            self.figure_line_label_axis_tick.setStyleSheet(font_style)
            self.figure_line_label_axis_title.setText('Axis Title Size')
            self.figure_line_label_axis_title.setStyleSheet(font_style)
            self.figure_line_label_legend.setText('Legend Size')
            self.figure_line_label_legend.setStyleSheet(font_style)
            self.figure_line_label_x.setText('X-Title')
            self.figure_line_label_x.setStyleSheet(font_style)
            self.figure_line_label_y.setText('Y-Title')
            self.figure_line_label_y.setStyleSheet(font_style)
            self.figure_line_label_x_range.setText('X-Range')
            self.figure_line_label_x_range.setStyleSheet(font_style)
            self.figure_line_label_y_range.setText('Y-Range')
            self.figure_line_label_y_range.setStyleSheet(font_style)
            self.figure_line_label_marker.setText('Marker Size')
            self.figure_line_label_marker.setStyleSheet(font_style)
            self.figure_line_label_color.setText('Color')
            self.figure_line_label_color.setStyleSheet(font_style)
            self.figure_line_btn_color_2.setText('Select Color')
            self.figure_line_btn_color_2.setText('Select Color')
            ## Bar
            self.figure_bar_label_axis_tick.setText('Axis Tick Size')
            self.figure_bar_label_axis_tick.setStyleSheet(font_style)
            self.figure_bar_label_axis_title.setText('Axis Title Size')
            self.figure_bar_label_axis_title.setStyleSheet(font_style)
            self.figure_bar_label_x.setText('X-Title')
            self.figure_bar_label_x.setStyleSheet(font_style)
            self.figure_bar_label_y.setText('Y-Title')
            self.figure_bar_label_y.setStyleSheet(font_style)
            self.figure_bar_label_y_range.setText('Y-Range')
            self.figure_bar_label_y_range.setStyleSheet(font_style)
            self.figure_bar_label_trend.setText('Trend Line')
            self.figure_bar_label_trend.setStyleSheet(font_style)
            self.figure_bar_label_bar.setText('Bar Value')
            self.figure_bar_label_bar.setStyleSheet(font_style)
            self.figure_bar_label_color.setText('Color')
            self.figure_bar_label_color.setStyleSheet(font_style)
            self.figure_bar_label_error.setText('Error Bar')
            self.figure_bar_label_error.setStyleSheet(font_style)
            self.figure_bar_btn_trend_2.setText('Select Color')

            self.figure_bar_btn_color_2.setText('Select Color')
            ## Scatter
            self.figure_scatter_label_legend.setText('Legend Size')
            self.figure_scatter_label_legend.setStyleSheet(font_style)
            self.figure_scatter_label_range.setText('Value Range')
            self.figure_scatter_label_range.setStyleSheet(font_style)
            self.figure_scatter_label_color.setText('Color Type')
            self.figure_scatter_label_color.setStyleSheet(font_style)
            self.figure_scatter_label_shape.setText('Shape')
            self.figure_scatter_label_shape.setStyleSheet(font_style)
            self.figure_scatter_label_shape_size.setText('Shape Size')
            self.figure_scatter_label_shape_size.setStyleSheet(font_style)
            self.figure_scatter_btn_shape_2.setText('Select Shape')

            # Analysis
            self.editStructure.setPlaceholderText('Supported file formats: GRO PDB')
            self.editTrajectory.setPlaceholderText('Supported file formats: XTC TPR')
            self.editResult.setPlaceholderText('Please select the result save path')
            self.btnSructure.setText('Structure')

            self.btnResult.setText('Result')
            self.btnTrajectory.setText('Trajectory')

            self.analysis_label_frame.setText('Frame')
            self.analysis_label_frame.setStyleSheet(font_style)
            self.analysis_label_first.setText('First')
            self.analysis_label_first.setStyleSheet(font_style)
            self.analysis_label_last.setText('Last')
            self.analysis_label_last.setStyleSheet(font_style)
            self.analysis_label_step.setText('Step')
            self.analysis_label_step.setStyleSheet(font_style)
            self.analysis_label_k.setText('K-Neighboors')
            self.analysis_label_k.setStyleSheet(font_style)
            self.analysis_label_method.setText('Analysis Method')
            self.analysis_label_method.setStyleSheet(font_style)


# ///////////////////////////导入文件////////////////////////////////////////
class BtnGetPath:
    __slots__ = ()
    TYPE = {
        'analysis_gro':{
            'Qtype':QFileDialog.getOpenFileName
            ,'caption':"Select Structure File"
            ,'dir':""
            ,'filter':"GRO Files (*.gro);;PDB Files (*.pdb)"
        }
        ,'analysis_xtc':{
            'Qtype': QFileDialog.getOpenFileName
            , 'caption': "Select Trajectory File"
            , 'dir': ""
            , 'filter': "XTC Files (*.xtc);;tpr Files (*.tpr)"
        }
        ,'analysis_result': {
            'Qtype': QFileDialog.getSaveFileName,
            'caption': "Save CSV File",
            'dir': QDir.currentPath() + "/untitled.csv",  # 默认保存路径和文件名
            'filter': "CSV Files (*.csv)"  # 文件类型过滤器
}
        ,'gene_gro':{
            'Qtype': QFileDialog.getSaveFileName
            , 'caption': "Save gro file"
            , 'dir': QDir.currentPath() + "/untitled.gro"
            , 'filter': "Gro Files (*gro)"
        }
        ,'figure_xlsx':{
            'Qtype': QFileDialog.getOpenFileName
            , 'caption': "Select csv file"
            , 'dir': ''
            , 'filter': "CSV Files (*.csv)"
        }
        ,'multiple_files': {
            'Qtype': QFileDialog.getOpenFileNames,  # 用于选择多个文件
            'caption': "Select Files",
            'dir': "",
            'filter': "All Files (*)"
        }
        , 'data_save':{
            'Qtype': QFileDialog.getSaveFileName
            , 'caption': "Save Data"
            , 'dir': QDir.currentPath() + "/untitled.xlsx"
            , 'filter': "Excel Files (*.xlsx)"
        }
    }

    @classmethod
    def open_file(cls, Qtype, caption, dir, filter):
        path, _ = Qtype(widgets
                        , caption=caption
                        , dir=dir
                        , filter=filter
                        , options=QFileDialog.DontUseNativeDialog)
        return path

    @classmethod
    def run(cls, label, name, edit=None):
        path = cls.open_file(**cls.TYPE[name])
        print('path:', path)
        if isinstance(path, list):
            label.setText(f'{len(path)} files were selected!')
            edit.clear()
            edit.setText('\n'.join(path))
        else:
            label.setText(path) if path else label.setText('Import file failed!!')


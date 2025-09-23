from collections import defaultdict
from functools import partial

from PySide6.QtGui import QColor
from PySide6.QtWidgets import QDoubleSpinBox, QLineEdit, QRadioButton, QComboBox

from figure.figure import *
from .Tools import *
from _exception import ExcelFormatError

global window


TYPE_ID = {'Height(nm)': 0
               , 'Sz Order Parameter (chain: sn1 and sn2)': 0  # 添加Sz Order Parameter支持
               , 'Sz Order Parameter (chain: sn1)': 0  # 添加Sz Order Parameter支持
               , 'Sz Order Parameter (chain: sn2)': 0  # 添加Sz Order Parameter支持
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
               , 'DensityTime': 2
               , 'Density With Time': 2  # 新的Density With Time类型
               , 'Density With Radius': 3  # 新的Density With Radius类型
           }

FIGURE_TYPE = {0: ['Bar', 'Line', 'Scatter', 'Map']  # Lipids
    , 1: ['Bar', 'Line']  # bubble
    , 2: ['Line']  # other
    , 3: ['Line', 'Heatmap']}  # Density


class BaseParameterReader:
    """参数读取器基类"""
    
    def __init__(self, ui_instance):
        self.ui = ui_instance
        # 文件信息管理
        self.path_figure = ui_instance.figure_edit_path.text()
        self.description, self.results, self.time_unit = read_excel(self.path_figure)
        self.lipids_type = self.results['Resname'].unique() if 'Resname' in self.results.columns else None
        self.residues = None
        
        # 参数存储
        self.LineInfo = defaultdict(lambda: None)
        self.BarInfo = defaultdict(lambda: None)
        self.ScaInfo = defaultdict(lambda: None)
        self.MapInfo = defaultdict(lambda: None)
        self.ColorInfo = defaultdict(lambda: None)
        self.ShapeInfo = []
        
        # 图表方法
        self.figureMethod = None
    
    def _get_current_tab_widget(self):
        """获取当前活动的TabWidget"""
        if hasattr(self.ui, 'tab_manager') and hasattr(self.ui.tab_manager, 'current_tab_widget'):
            return self.ui.tab_manager.current_tab_widget
        return None
    
    def _get_current_tab(self):
        """获取当前活动的标签页"""
        current_tab_widget = self._get_current_tab_widget()
        if current_tab_widget:
            return current_tab_widget.currentWidget()
        return None
    
    def _get_file_type_from_description(self):
        """根据描述判断文件类型"""
        if 'Bubble' in self.description:
            return 'bubble'
        elif 'Density With Time' in self.description:
            return 'density_time'
        elif 'Density With Radius' in self.description:
            return 'density_radius'
        else:
            return 'lipids'
    
    def _get_time_unit_label(self):
        """根据CSV文件中的TIME_UNIT注释获取时间单位标签"""
        try:
            if self.time_unit == 'ns':
                return "Time (ns)"
            elif self.time_unit == 'ps':
                return "Time (ps)"
            else:
                return "Frames"
        except Exception:
            return "Frames"


class LipidsParameterReader(BaseParameterReader):
    """脂质数据参数读取器"""
    
    def get_line(self):
        """读取Line图参数"""
        current_tab = self._get_current_tab()
        if current_tab:
            # 从当前标签页中查找控件
            axis_scale_spin = current_tab.findChild(QDoubleSpinBox, "lipids_line_spin_axis_tick")
            axis_text_spin = current_tab.findChild(QDoubleSpinBox, "lipids_line_spin_axis_title")
            grid_size_spin = current_tab.findChild(QDoubleSpinBox, "lipids_line_spin_legend")
            x_title_edit = current_tab.findChild(QLineEdit, "lipids_line_edit_x")
            y_title_edit = current_tab.findChild(QLineEdit, "lipids_line_edit_y")
            x_min_spin = current_tab.findChild(QDoubleSpinBox, "lipids_line_spin_x_min")
            x_max_spin = current_tab.findChild(QDoubleSpinBox, "lipids_line_spin_x_max")
            y_min_spin = current_tab.findChild(QDoubleSpinBox, "lipids_line_spin_y_min")
            y_max_spin = current_tab.findChild(QDoubleSpinBox, "lipids_line_spin_y_max")
            marker_size_spin = current_tab.findChild(QDoubleSpinBox, "lipids_line_spin_marker")
            marker_shape_combo = None  # file_detection.py中没有创建这个控件
            
            # 读取参数值
            line_info = {
                'axis_scale': axis_scale_spin.value() if axis_scale_spin else 12.0,
                'axis_text': axis_text_spin.value() if axis_text_spin else 16.0,
                'grid_size': grid_size_spin.value() if grid_size_spin else 1.0,
                'x_title': x_title_edit.text() if x_title_edit and x_title_edit.text() else self._get_time_unit_label(),
                'y_title': y_title_edit.text() if y_title_edit and y_title_edit.text() else self.description,
                'x_min': x_min_spin.value() if x_min_spin else 0.0,
                'x_max': x_max_spin.value() if x_max_spin else 0.0,
                'y_min': y_min_spin.value() if y_min_spin else 0.0,
                'y_max': y_max_spin.value() if y_max_spin else 1.0,
                'marker_size': marker_size_spin.value() if marker_size_spin else 0,
                'marker_shape': marker_shape_combo.currentText() if marker_shape_combo else 'o',
                'color': self.ColorInfo
            }
        else:
            line_info = {}
        return line_info
    
    def get_bar(self):
        """读取Bar图参数"""
        current_tab = self._get_current_tab()
        if current_tab:
            # 从当前标签页中查找控件
            axis_tick_spin = current_tab.findChild(QDoubleSpinBox, "lipids_bar_spin_axis_tick")
            axis_title_spin = current_tab.findChild(QDoubleSpinBox, "lipids_bar_spin_axis_title")
            x_title_edit = current_tab.findChild(QLineEdit, "lipids_bar_edit_x")
            y_title_edit = current_tab.findChild(QLineEdit, "lipids_bar_edit_y")
            y_min_spin = current_tab.findChild(QDoubleSpinBox, "lipids_bar_spin_y_min")
            y_max_spin = current_tab.findChild(QDoubleSpinBox, "lipids_bar_spin_y_max")
            error_yes_radio = current_tab.findChild(QRadioButton, "lipids_bar_radio_error_yes")
            error_no_radio = current_tab.findChild(QRadioButton, "lipids_bar_radio_error_no")
            
            # 读取参数值
            bar_info = {
                'axis_scale': axis_tick_spin.value() if axis_tick_spin else 12.0,
                'axis_text': axis_title_spin.value() if axis_title_spin else 16.0,
                'x_title': x_title_edit.text() if x_title_edit and x_title_edit.text() else self._get_time_unit_label(),
                'y_title': y_title_edit.text() if y_title_edit and y_title_edit.text() else self.description,
                'y_min': y_min_spin.value() if y_min_spin else 0.0,
                'y_max': y_max_spin.value() if y_max_spin else 1.0,
                'trend_size': 0.0,  # 暂时设为默认值
                'up_bar_value': 0.0,  # 暂时设为默认值
                'error_deci': error_yes_radio.isChecked() if error_yes_radio else False,
                'color': self.ColorInfo
            }
        else:
            bar_info = {}
        return bar_info
    
    def get_scatter(self):
        """读取Scatter图参数"""
        current_tab = self._get_current_tab()
        if current_tab:
            # 从当前标签页中查找控件
            grid_size_spin = current_tab.findChild(QDoubleSpinBox, "lipids_scatter_spin_axis_tick")
            color_min_spin = current_tab.findChild(QDoubleSpinBox, "lipids_scatter_spin_range_min")
            color_max_spin = current_tab.findChild(QDoubleSpinBox, "lipids_scatter_spin_range_max")
            color_map_combo = None  # file_detection.py中没有创建这个控件
            shape_size_spin = current_tab.findChild(QDoubleSpinBox, "lipids_scatter_spin_shape_size")
            
            # 读取参数值
            scatter_info = {
                'grid_size': grid_size_spin.value() if grid_size_spin else 1.0,
                'bar_min': color_min_spin.value() if color_min_spin else 0.0,
                'bar_max': color_max_spin.value() if color_max_spin else 1.0,
                'bar_color': color_map_combo.currentText() if color_map_combo else 'viridis',
                'shape_size': shape_size_spin.value() if shape_size_spin else 50.0,
                'shape': {}
            }
            
            # 读取形状信息
            for group in self.ShapeInfo:
                scatter_info['shape'][group.title()] = next(
                    (radio.text() for radio in group.findChildren(QRadioButton) if radio.isChecked()), None
                )
        else:
            scatter_info = {}
        return scatter_info
    
    def get_map(self):
        """读取Map图参数"""
        current_tab = self._get_current_tab()
        if current_tab:
            # 从当前标签页中查找控件
            color_bar_yes_radio = current_tab.findChild(QRadioButton, "lipids_map_radio_color_bar_yes")
            color_bar_no_radio = current_tab.findChild(QRadioButton, "lipids_map_radio_color_bar_no")
            x_title_edit = current_tab.findChild(QLineEdit, "lipids_map_edit_x")
            y_title_edit = current_tab.findChild(QLineEdit, "lipids_map_edit_y")
            value_min_spin = current_tab.findChild(QDoubleSpinBox, "lipids_map_spin_value_min")
            value_max_spin = current_tab.findChild(QDoubleSpinBox, "lipids_map_spin_value_max")
            color_map_combo = current_tab.findChild(QComboBox, "lipids_map_combo_color_map")
            
            # 读取参数值
            map_info = {
                'color_bar': color_bar_yes_radio.isChecked() if color_bar_yes_radio else True,
                'x_title': x_title_edit.text() if x_title_edit and x_title_edit.text() else "X Position (nm)",
                'y_title': y_title_edit.text() if y_title_edit and y_title_edit.text() else "Y Position (nm)",
                'value_min': value_min_spin.value() if value_min_spin else 0.0,
                'value_max': value_max_spin.value() if value_max_spin else 1.0,
                'color_map': color_map_combo.currentText() if color_map_combo else "viridis"
            }
        else:
            map_info = {}
        return map_info
    


class BubbleParameterReader(BaseParameterReader):
    """气泡数据参数读取器"""
    
    def get_line(self):
        """读取Line图参数"""
        current_tab = self._get_current_tab()
        if current_tab:
            # 从当前标签页中查找控件
            axis_scale_spin = current_tab.findChild(QDoubleSpinBox, "bubble_line_spin_axis_tick")
            axis_text_spin = current_tab.findChild(QDoubleSpinBox, "bubble_line_spin_axis_title")
            grid_size_spin = current_tab.findChild(QDoubleSpinBox, "bubble_line_spin_legend")
            x_title_edit = current_tab.findChild(QLineEdit, "bubble_line_edit_x")
            y_title_edit = current_tab.findChild(QLineEdit, "bubble_line_edit_y")
            x_min_spin = current_tab.findChild(QDoubleSpinBox, "bubble_line_spin_x_min")
            x_max_spin = current_tab.findChild(QDoubleSpinBox, "bubble_line_spin_x_max")
            y_min_spin = current_tab.findChild(QDoubleSpinBox, "bubble_line_spin_y_min")
            y_max_spin = current_tab.findChild(QDoubleSpinBox, "bubble_line_spin_y_max")
            marker_size_spin = current_tab.findChild(QDoubleSpinBox, "bubble_line_spin_marker")
            marker_shape_combo = None  # file_detection.py中没有创建这个控件
            
            # 读取参数值
            line_info = {
                'axis_scale': axis_scale_spin.value() if axis_scale_spin else 12.0,
                'axis_text': axis_text_spin.value() if axis_text_spin else 16.0,
                'grid_size': grid_size_spin.value() if grid_size_spin else 1.0,
                'x_title': x_title_edit.text() if x_title_edit and x_title_edit.text() else self._get_time_unit_label(),
                'y_title': y_title_edit.text() if y_title_edit and y_title_edit.text() else self.description,
                'x_min': x_min_spin.value() if x_min_spin else 0.0,
                'x_max': x_max_spin.value() if x_max_spin else 100.0,
                'y_min': y_min_spin.value() if y_min_spin else 0.0,
                'y_max': y_max_spin.value() if y_max_spin else 1.0,
                'marker_size': marker_size_spin.value() if marker_size_spin else 0,
                'marker_shape': marker_shape_combo.currentText() if marker_shape_combo else 'o',
                'color': self.ColorInfo
            }
        else:
            line_info = {}
        return line_info
    
    def get_bar(self):
        """读取Bar图参数"""
        current_tab = self._get_current_tab()
        if current_tab:
            # 从当前标签页中查找控件
            axis_tick_spin = current_tab.findChild(QDoubleSpinBox, "bubble_bar_spin_axis_tick")
            axis_title_spin = current_tab.findChild(QDoubleSpinBox, "bubble_bar_spin_axis_title")
            x_title_edit = current_tab.findChild(QLineEdit, "bubble_bar_edit_x")
            y_title_edit = current_tab.findChild(QLineEdit, "bubble_bar_edit_y")
            y_min_spin = current_tab.findChild(QDoubleSpinBox, "bubble_bar_spin_y_min")
            y_max_spin = current_tab.findChild(QDoubleSpinBox, "bubble_bar_spin_y_max")
            error_yes_radio = current_tab.findChild(QRadioButton, "bubble_bar_radio_error_yes")
            error_no_radio = current_tab.findChild(QRadioButton, "bubble_bar_radio_error_no")
            
            # 读取参数值
            bar_info = {
                'axis_scale': axis_tick_spin.value() if axis_tick_spin else 12.0,
                'axis_text': axis_title_spin.value() if axis_title_spin else 16.0,
                'x_title': x_title_edit.text() if x_title_edit and x_title_edit.text() else self._get_time_unit_label(),
                'y_title': y_title_edit.text() if y_title_edit and y_title_edit.text() else self.description,
                'y_min': y_min_spin.value() if y_min_spin else 0.0,
                'y_max': y_max_spin.value() if y_max_spin else 1.0,
                'trend_size': 0.0,  # 暂时设为默认值
                'up_bar_value': 0.0,  # 暂时设为默认值
                'error_deci': error_yes_radio.isChecked() if error_yes_radio else False,
                'color': self.ColorInfo
            }
        else:
            bar_info = {}
        return bar_info
    
    def get_heatmap(self):
        """读取Heatmap图参数"""
        current_tab = self._get_current_tab()
        if current_tab:
            # 从当前标签页中查找控件
            x_title_edit = current_tab.findChild(QLineEdit, "density_heatmap_edit_x")
            y_title_edit = current_tab.findChild(QLineEdit, "density_heatmap_edit_y")
            color_map_combo = current_tab.findChild(QComboBox, "density_heatmap_combo_color_map")
            
            # 读取参数值
            heatmap_info = {
                'x_title': x_title_edit.text() if x_title_edit and x_title_edit.text() else self._get_time_unit_label(),
                'y_title': y_title_edit.text() if y_title_edit and y_title_edit.text() else self.description,
                'color_map': color_map_combo.currentText() if color_map_combo else "viridis"
            }
        else:
            heatmap_info = {}
        return heatmap_info
    
    def get_heatmap(self):
        """读取Heatmap图参数 - BubbleParameterReader版本"""
        current_tab = self._get_current_tab()
        if current_tab:
            # 从当前标签页中查找控件
            x_title_edit = current_tab.findChild(QLineEdit, "density_heatmap_edit_x")
            y_title_edit = current_tab.findChild(QLineEdit, "density_heatmap_edit_y")
            color_map_combo = current_tab.findChild(QComboBox, "density_heatmap_combo_color_map")
            
            # 读取参数值
            heatmap_info = {
                'x_title': x_title_edit.text() if x_title_edit and x_title_edit.text() else self._get_time_unit_label(),
                'y_title': y_title_edit.text() if y_title_edit and y_title_edit.text() else self.description,
                'cmap': color_map_combo.currentText() if color_map_combo else "viridis"
            }
        else:
            heatmap_info = {}
        return heatmap_info


def create_parameter_reader(ui):
    """
    根据文件类型创建相应的参数读取器
    
    Args:
        ui: UI实例
        
    Returns:
        LipidsParameterReader或BubbleParameterReader实例
    """
    try:
        # 优先使用已保存的文件类型检测结果
        if hasattr(ui, 'detected_file_type') and ui.detected_file_type:
            file_type = ui.detected_file_type
        else:
            # 如果没有保存的结果，则进行文件检测
            from .file_detection import FileTypeDetector
            file_type = FileTypeDetector.detect_file_type(ui.figure_edit_path.text())
            
        # 如果文件检测器返回unknown，则回退到description判断
        if file_type in ['unknown_csv_format', 'unknown_csv_type']:
            file_type = _get_file_type_from_description(ui)
            
        if file_type == 'lipids':
            return LipidsParameterReader(ui)
        elif file_type == 'bubble':
            return BubbleParameterReader(ui)
        elif file_type == 'density_time':
            return BubbleParameterReader(ui)  # 使用bubble类型的参数读取器
        elif file_type == 'density_radius':
            return BubbleParameterReader(ui)  # 使用bubble类型的参数读取器
        else:
            return None
    except Exception as e:
        print(f"创建参数读取器失败: {e}")
        ui.figure_extra = 0
        return None

def _get_file_type_from_description(ui):
    """根据描述判断文件类型"""
    try:
        description, _, _ = read_excel(ui.figure_edit_path.text())
        if 'Bubble' in description:
            return 'bubble'
        elif 'Density With Time' in description or 'DensityTime' in description:
            return 'density_time'
        elif 'Density With Radius' in description or 'DensityRadius' in description or 'Multi-Radius Density Analysis' in description:
            return 'density_radius'
        else:
            return 'lipids'
    except Exception as e:
        print(f"Exception in _get_file_type_from_description: {e}")
        return 'lipids'

def _get_current_chart_type(ui):
    """从当前TabWidget获取图表类型"""
    try:
        # 检查是否有tab_manager
        if hasattr(ui, 'tab_manager') and hasattr(ui.tab_manager, 'current_tab_widget'):
            current_tab_widget = ui.tab_manager.current_tab_widget
            if current_tab_widget:
                current_index = current_tab_widget.currentIndex()
                # 根据文件类型和索引确定图表类型
                file_type = _get_file_type_from_description(ui)
                if file_type == 'lipids':
                    chart_types = ['Line', 'Bar', 'Scatter', 'Map']
                elif file_type == 'bubble':
                    chart_types = ['Line', 'Bar']
                elif file_type == 'density_time':
                    chart_types = ['Line', 'Heatmap']
                elif file_type == 'density_radius':
                    chart_types = ['Line', 'Heatmap']
                else:
                    chart_types = ['Line', 'Bar']
                
                if 0 <= current_index < len(chart_types):
                    return chart_types[current_index]
        
        # 如果无法获取，返回默认值
        print("无法获取当前图表类型，使用默认的Line")
        return 'Line'
    except Exception as e:
        print(f"获取当前图表类型失败: {e}")
        return 'Line'

def get_xlsx_path(ui) -> str:
    """获取Excel文件路径"""
    return ui.figure_edit_path.text()


class FigurePage:

    @staticmethod
    def _get_time_unit_label(figure_info):
        """根据CSV文件中的TIME_UNIT注释获取时间单位标签"""
        try:
            # 从FigureInfo中获取时间单位
            time_unit = getattr(figure_info, 'time_unit', 'frame')
            
            # 根据时间单位返回相应的标签
            if time_unit == 'ns':
                return "Time (ns)"
            elif time_unit == 'ps':
                return "Time (ps)"
            else:
                return "Frames"
                
        except Exception:
            return "Frames"

    @staticmethod
    def _ensure_figure_info(ui):
        """
        用来检测UI界面是否导入了结果的文件路径，以及结果文件路径有没有发生更改
        """
        if ui.FigureInfo is None or ui.FigureInfo.path_figure != ui.figure_edit_path.text():
            ui.FigureInfo = create_parameter_reader(ui)
        return ui.FigureInfo

    @classmethod
    def figureBtnMakeFigure(cls, ui):
        try:
            cls._ensure_figure_info(ui)
        except Exception as e:
            create_warn_dialog("Please select the result file first.\nError in the function:figureBtnMakeFigure")
            return
        
        # 从当前TabWidget获取图表类型
        chart_type = _get_current_chart_type(ui)

        # 优先使用已保存的文件类型检测结果
        if hasattr(ui, 'detected_file_type') and ui.detected_file_type:
            file_type = ui.detected_file_type
            print(f"DEBUG: Using cached file_type: {file_type}")
        else:
            # 如果没有保存的结果，则进行文件检测
            print(f"DEBUG: Detecting file type from path: {ui.FigureInfo.path_figure}")
            from .file_detection import FileTypeDetector
            file_type = FileTypeDetector.detect_file_type(ui.FigureInfo.path_figure)
            print(f"DEBUG: Detected file_type: {file_type}")
            
        # 如果文件检测器返回unknown，则回退到description判断
        if file_type in ['unknown_csv_format', 'unknown_csv_type']:
            print(f"DEBUG: File type unknown, using description fallback")
            print(f"DEBUG: Description: {ui.FigureInfo.description}")
            file_type = cls._get_file_type_from_data(ui.FigureInfo.description)
            print(f"DEBUG: Fallback file_type: {file_type}")
        
        if file_type == 'lipids':
            from figure.figure import LipidsFigure
            chart_settings = cls._get_chart_settings(ui, chart_type)
            lipids_figure = LipidsFigure(ui.FigureInfo.description, ui.FigureInfo.results, chart_settings)
            
            if chart_type == 'Line':
                lipids_figure.plot_line()
            elif chart_type == 'Bar':
                lipids_figure.plot_bar()
            elif chart_type == 'Scatter':
                lipids_figure.plot_scatter()
            elif chart_type == 'Map':
                lipids_figure.plot_map()
                
        elif file_type == 'bubble':
            from figure.figure import BubbleFigure
            chart_settings = cls._get_chart_settings(ui, chart_type)
            bubble_figure = BubbleFigure(ui.FigureInfo.description, ui.FigureInfo.results, chart_settings)
            
            if chart_type == 'Line':
                bubble_figure.plot_line()
            elif chart_type == 'Bar':
                bubble_figure.plot_bar()
                
        elif file_type == 'density':
            from figure.density_figure import DensityFigure
            chart_settings = cls._get_chart_settings(ui, chart_type)
            density_figure = DensityFigure(ui.FigureInfo.description, ui.FigureInfo.results, chart_settings)
            
            if chart_type == 'Line':
                density_figure.plot_line()
            elif chart_type == 'Heatmap':
                density_figure.plot_heatmap()
                
        elif file_type == 'density_time':
            from figure.density_figure import DensityFigure
            chart_settings = cls._get_chart_settings(ui, chart_type)
            density_figure = DensityFigure(ui.FigureInfo.description, ui.FigureInfo.results, chart_settings)
            
            if chart_type == 'Line':
                density_figure.plot_line()
            elif chart_type == 'Heatmap':
                density_figure.plot_heatmap()
                
        elif file_type == 'density_radius':
            from figure.density_figure import DensityFigure
            chart_settings = cls._get_chart_settings(ui, chart_type)
            density_figure = DensityFigure(ui.FigureInfo.description, ui.FigureInfo.results, chart_settings)
            
            if chart_type == 'Line':
                density_figure.plot_line()
            elif chart_type == 'Heatmap':
                density_figure.plot_heatmap()
    
    @staticmethod
    def _get_file_type_from_data(description):
        """根据描述判断文件类型"""
        if 'Bubble' in description:
            return 'bubble'
        elif 'DensityTime' in description:
            return 'density_time'
        elif 'DensityRadius' in description:
            return 'density_radius'
        elif 'Density' in description:
            return 'density'
        else:
            return 'lipids'
    
    @staticmethod
    def _get_chart_settings(ui, chart_type):
        """根据图表类型获取相应的设置"""
        if ui.FigureInfo:
            if chart_type == 'Line':
                return ui.FigureInfo.get_line()
            elif chart_type == 'Bar':
                return ui.FigureInfo.get_bar()
            elif chart_type == 'Scatter' and hasattr(ui.FigureInfo, 'get_scatter'):
                return ui.FigureInfo.get_scatter()
            elif chart_type == 'Map' and hasattr(ui.FigureInfo, 'get_map'):
                return ui.FigureInfo.get_map()
            elif chart_type == 'Heatmap' and hasattr(ui.FigureInfo, 'get_heatmap'):
                return ui.FigureInfo.get_heatmap()
        return {}

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
            if hasattr(ui, 'figure_color_extra_box'):
                ui.figure_color_extra_box.setStyleSheet(u"font: 15pt \"\u534e\u6587\u7ec6\u9ed1\";")
                # 修复：设置合适的尺寸，让颜色选择框可见
                ui.figure_color_extra_box.setMinimumSize(200, 300)
                ui.figure_color_extra_box.setMaximumSize(300, 500)
                ui.FigureColorLayout = QVBoxLayout()
                ui.figure_color_extra_box.setLayout(ui.FigureColorLayout)
            else:
                raise AttributeError("figure_color_extra_box不存在")
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
                    # , background_color='white'
                    # , font_color='black'
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
        try:
            ui.figure_extra = 1
            
            cls._ensure_color_lipids_sel(ui)  # 确保颜色选择布局已经创建

            cls._make_residues_btn(ui, cls.residues_color_sel)

            if ui.FigureInfo.lipids_type is None:
                raise ExcelFormatError("请严格遵守结果文件的格式！")
            
            # 得到残基的信息
            if hasattr(ui, 'figure_color_extra_box'):
                ui.figure_color_extra_box.show()
                ui.figure_color_extra_box.raise_()  # 确保在最前面
            else:
                pass  # 警告：figure_color_extra_box不存在
                
        except Exception as e:
            raise

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
        btn.setStyleSheet(
            f'background-color: rgb({new_color[0]*255}, {new_color[1]*255}, {new_color[2]*255});'
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
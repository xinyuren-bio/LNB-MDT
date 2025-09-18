# -*- coding: utf-8 -*-
"""
按钮功能绑定辅助工具

这个模块提供了一个通用的按钮绑定系统，用于为动态创建的TabWidget中的按钮绑定功能。

使用方法：
1. 创建ButtonBindingHelper实例
2. 注册按钮功能映射
3. 调用bind_buttons方法绑定按钮

示例：
    helper = ButtonBindingHelper(ui_instance)
    helper.register_button_function("color_btn", handle_color_selection)
    helper.bind_buttons(tab_widget, file_type, file_path)
"""

from PySide6.QtWidgets import QPushButton, QMessageBox
from typing import Dict, Callable, Any, Optional


class ButtonBindingHelper:
    """按钮功能绑定辅助类"""
    
    def __init__(self, ui_instance):
        """
        初始化按钮绑定辅助类
        
        Args:
            ui_instance: 主UI实例
        """
        self.ui_instance = ui_instance
        self.button_functions: Dict[str, Callable] = {}
        self.file_type: Optional[str] = None
        
    def register_button_function(self, button_name: str, function: Callable):
        """
        注册按钮功能
        
        Args:
            button_name: 按钮名称（用于查找按钮）
            function: 按钮点击时调用的函数
        """
        self.button_functions[button_name] = function
        
    def register_multiple_functions(self, function_map: Dict[str, Callable]):
        """
        批量注册按钮功能
        
        Args:
            function_map: 按钮名称到函数的映射字典
        """
        self.button_functions.update(function_map)
        
    def bind_buttons(self, tab_widget, file_type: str, file_path: str):
        """
        为TabWidget绑定按钮功能
        
        Args:
            tab_widget: TabWidget实例
            file_type: 文件类型
            file_path: 文件路径
        """
        self.file_type = file_type
        
        # 绑定tab切换事件
        tab_widget.currentChanged.connect(
            lambda index: self._bind_current_tab_buttons(tab_widget, file_path)
        )
        
        # 立即绑定当前活跃tab的按钮
        self._bind_current_tab_buttons(tab_widget, file_path)
        
    def _bind_current_tab_buttons(self, tab_widget, file_path: str):
        """绑定当前活跃tab的按钮"""
        current_index = tab_widget.currentIndex()
        current_tab = tab_widget.widget(current_index)
        
        if current_tab:
            self._bind_tab_buttons(current_tab, current_index, file_path)
            
    def _bind_tab_buttons(self, tab_widget, tab_index: int, file_path: str):
        """
        绑定特定tab的按钮
        
        Args:
            tab_widget: tab widget实例
            tab_index: tab索引
            file_path: 文件路径
        """
        # 根据文件类型和tab索引确定按钮映射
        button_mapping = self._get_button_mapping(tab_index)
        
        for button_name, function_name in button_mapping.items():
            button = tab_widget.findChild(QPushButton, button_name)
            if button and function_name in self.button_functions:
                # 断开之前的连接（避免重复绑定）
                try:
                    button.clicked.disconnect()
                except:
                    pass
                # 绑定新功能
                button.clicked.connect(
                    lambda checked=False, func=self.button_functions[function_name]: 
                    self._safe_call_function(func, file_path)
                )
                
    def _get_button_mapping(self, tab_index: int) -> Dict[str, str]:
        """
        根据tab索引获取按钮映射
        
        Args:
            tab_index: tab索引
            
        Returns:
            按钮名称到函数名称的映射
        """
        if self.file_type == 'lipids':
            return self._get_lipids_button_mapping(tab_index)
        elif self.file_type == 'bubble':
            return self._get_bubble_button_mapping(tab_index)
        else:
            return {}
            
    def _get_lipids_button_mapping(self, tab_index: int) -> Dict[str, str]:
        """获取Lipids类型的按钮映射"""
        mappings = {
            0: {  # Line tab
                "lipids_line_btn_color": "handle_color_selection",
            },
            1: {  # Bar tab
                "lipids_bar_btn_color": "handle_color_selection",
                "lipids_bar_btn_trend": "handle_trend_selection",
            },
            2: {  # Scatter tab
                "lipids_scatter_btn_shape": "handle_shape_selection",
            }
        }
        return mappings.get(tab_index, {})
        
    def _get_bubble_button_mapping(self, tab_index: int) -> Dict[str, str]:
        """获取Bubble类型的按钮映射"""
        mappings = {
            0: {  # Line tab
                "bubble_line_btn_color": "handle_color_selection",
            },
            1: {  # Bar tab
                "bubble_bar_btn_color": "handle_color_selection",
                "bubble_bar_btn_trend": "handle_trend_selection",
            }
        }
        return mappings.get(tab_index, {})
        
    def _safe_call_function(self, function: Callable, file_path: str):
        """
        安全调用函数，包含错误处理
        
        Args:
            function: 要调用的函数
            file_path: 文件路径
        """
        try:
            function(file_path)
        except Exception as e:
            QMessageBox.critical(None, "功能执行错误", f"执行功能时出错: {str(e)}")
            
    def add_custom_button_mapping(self, file_type: str, tab_index: int, mapping: Dict[str, str]):
        """
        添加自定义按钮映射
        
        Args:
            file_type: 文件类型
            tab_index: tab索引
            mapping: 按钮名称到函数名称的映射
        """
        # 这个方法可以用于扩展新的文件类型或自定义按钮映射
        # 具体实现可以根据需要扩展
        pass


class DefaultButtonHandlers:
    """默认按钮处理器"""
    
    def __init__(self, ui_instance):
        self.ui_instance = ui_instance
        
    def handle_color_selection(self, file_path: str):
        """处理颜色选择功能"""
        try:
            from .Fuctions_Figure import FigurePage
            
            # 确保FigureInfo已初始化
            if not hasattr(self.ui_instance, 'FigureInfo') or self.ui_instance.FigureInfo is None:
                QMessageBox.warning(None, "警告", "请先选择结果文件！")
                return
            
            # 根据文件类型调用相应的颜色选择方法
            # 这里需要根据实际的文件类型判断逻辑来实现
            FigurePage.lipids_colors(self.ui_instance)
                
        except Exception as e:
            QMessageBox.critical(None, "错误", f"颜色选择功能出错: {str(e)}")
    
    def handle_trend_selection(self, file_path: str):
        """处理趋势线选择功能"""
        try:
            from .Fuctions_Figure import FigurePage
            
            if not hasattr(self.ui_instance, 'FigureInfo') or self.ui_instance.FigureInfo is None:
                QMessageBox.warning(None, "警告", "请先选择结果文件！")
                return
            
            FigurePage.single_color(self.ui_instance)
            
        except Exception as e:
            QMessageBox.critical(None, "错误", f"趋势线选择功能出错: {str(e)}")
    
    def handle_shape_selection(self, file_path: str):
        """处理形状选择功能"""
        try:
            from .Fuctions_Figure import FigurePage
            
            if not hasattr(self.ui_instance, 'FigureInfo') or self.ui_instance.FigureInfo is None:
                QMessageBox.warning(None, "警告", "请先选择结果文件！")
                return
            
            FigurePage.figureBtnShape(self.ui_instance)
            
        except Exception as e:
            QMessageBox.critical(None, "错误", f"形状选择功能出错: {str(e)}")


def create_button_binding_helper(ui_instance) -> ButtonBindingHelper:
    """
    创建并配置按钮绑定辅助类
    
    Args:
        ui_instance: 主UI实例
        
    Returns:
        配置好的ButtonBindingHelper实例
    """
    helper = ButtonBindingHelper(ui_instance)
    handlers = DefaultButtonHandlers(ui_instance)
    
    # 注册默认按钮功能
    helper.register_multiple_functions({
        "handle_color_selection": handlers.handle_color_selection,
        "handle_trend_selection": handlers.handle_trend_selection,
        "handle_shape_selection": handlers.handle_shape_selection,
    })
    
    return helper

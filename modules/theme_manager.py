# -*- coding: utf-8 -*-
"""
主题管理模块

提供主题切换功能，支持白天/黑夜主题切换
"""

import os
from PySide6.QtWidgets import QPushButton, QHBoxLayout, QWidget, QApplication
from PySide6.QtCore import Qt, QSettings, QObject
from PySide6.QtGui import QIcon, QPixmap


class ThemeManager:
    """主题管理器"""
    
    def __init__(self, main_window):
        """
        初始化主题管理器
        
        Args:
            main_window: 主窗口实例
        """
        self.main_window = main_window
        self.current_theme = "light"  # 默认主题为白天
        self.theme_path = os.path.join(os.path.dirname(__file__), "..", "themes")
        
        # 加载保存的主题设置
        self.load_theme_setting()
        
    def load_theme_setting(self):
        """从设置中加载主题"""
        settings = QSettings("LNB-MDT", "Theme")
        saved_theme = settings.value("theme", "light")
        self.set_theme(saved_theme)
        
    def save_theme_setting(self):
        """保存主题设置"""
        settings = QSettings("LNB-MDT", "Theme")
        settings.setValue("theme", self.current_theme)
        
    def set_theme(self, theme_name):
        """
        设置主题
        
        Args:
            theme_name: 主题名称 ('dark' 或 'light')
        """
        if theme_name not in ["dark", "light"]:
            return False
            
        theme_file = os.path.join(self.theme_path, f"py_dracula_{theme_name}.qss")
        
        if not os.path.exists(theme_file):
            return False
            
        try:
            with open(theme_file, 'r', encoding='utf-8') as f:
                style_sheet = f.read()
            
            # 应用主题到UI的styleSheet组件
            if hasattr(self.main_window.ui, 'styleSheet'):
                self.main_window.ui.styleSheet.setStyleSheet(style_sheet)
            else:
                # 如果没有styleSheet组件，则应用到主窗口
                self.main_window.setStyleSheet(style_sheet)
            
            # 延迟应用按钮背景色，确保QSS已经加载完成
            from PySide6.QtCore import QTimer
            QTimer.singleShot(200, self._apply_button_backgrounds)
            
            # 主题切换完成，样式由QSS文件控制
            self.current_theme = theme_name
            
            # 保存设置
            self.save_theme_setting()
            
            # 通知UI组件更新图标
            self.notify_theme_change()
            
            return True
            
        except Exception as e:
            return False
    
    def _apply_button_backgrounds(self):
        """强制应用按钮背景色 - 确保背景色在macOS上显示（排除左侧和顶部按钮）"""
        try:
            from PySide6.QtWidgets import QPushButton
            from PySide6.QtGui import QPalette, QColor
            from PySide6.QtCore import QObject
            
            # 需要排除的按钮ID列表（左侧菜单和顶部控制按钮保持原始样式）
            excluded_object_names = {
                'btn_home', 'btn_generate', 'btn_figure', 'btn_analysis', 'btn_data_process',
                'toggleButton', 'btn_language', 'minimizeAppBtn', 'maximizeRestoreAppBtn', 
                'closeAppBtn', 'themeSwitchButton'
            }
            
            # 需要排除的父容器（包含这些按钮的容器）
            excluded_parents = set()
            # 查找 topMenu 和 rightButtons 容器中的按钮
            top_menu = self.main_window.findChild(QObject, "topMenu")
            right_buttons = self.main_window.findChild(QObject, "rightButtons")
            if top_menu:
                excluded_parents.add(top_menu)
            if right_buttons:
                excluded_parents.add(right_buttons)
            
            # 目标背景色和边框色
            target_bg = "#9faeda"
            
            for button in self.main_window.findChildren(QPushButton):
                # 跳过左侧菜单按钮和顶部控制按钮
                if button.objectName() in excluded_object_names:
                    continue
                # 跳过 topMenu 和 rightButtons 容器中的按钮
                parent = button.parent()
                if parent and parent in excluded_parents:
                    continue
                
                # 获取按钮的当前样式表
                current_style = button.styleSheet()
                
                # 确保按钮填充背景（macOS必需）
                button.setAutoFillBackground(True)
                
                # 使用QPalette设置背景色（macOS需要）
                palette = button.palette()
                palette.setColor(button.backgroundRole(), QColor(target_bg))
                button.setPalette(palette)
                
                # 重新应用样式表以保持border-radius等样式
                # 如果当前样式为空或不符合make_btn的样式，使用默认样式
                if current_style and ('!important' in current_style or 'border-radius' in current_style):
                    button.setStyleSheet(current_style)
                else:
                    # 使用make_btn的默认样式
                    default_style = (
                        f"background-color: {target_bg} !important;"
                        f"color: black !important;"
                        f"border: 2px solid {target_bg} !important;"
                        f"border-radius: 5px !important;"
                        f"font: 12pt \"华文细黑\";"
                        f"padding: 5px;"
                    )
                    button.setStyleSheet(default_style)
                
        except Exception as e:
            print(f"[Theme Manager] 应用按钮背景色时出错: {e}")
    
    def notify_theme_change(self):
        """
        通知UI组件主题已更改，需要更新图标
        """
        try:
            # 方法1：通过全局变量通知
            import main
            if hasattr(main, 'window') and main.window:
                # 查找AnalysisBtnClick实例
                if hasattr(main.window.ui, 'extraRightBox') and main.window.ui.extraRightBox:
                    if hasattr(main.window.ui, 'VLayoutRightMain'):
                        layout = main.window.ui.VLayoutRightMain
                        if layout:
                            for i in range(layout.count()):
                                item = layout.itemAt(i)
                                if item and item.widget():
                                    widget = item.widget()
                                    # 检查widget内部是否有按钮
                                    self.update_buttons_in_widget(widget)
            
            # 方法2：直接通过UI查找按钮并更新
            if hasattr(self.main_window.ui, 'btnBack') and hasattr(self.main_window.ui, 'btnRefresh'):
                self.update_direct_buttons()
            
        except Exception as e:
            pass  # 发送主题变更通知时出错
    
    def update_direct_buttons(self):
        """
        直接更新按钮图标
        """
        try:
            is_dark = self.is_dark_theme()
            
            if is_dark:
                back_icon_path = 'images/icons/houtui_black.png'
                refresh_icon_path = 'images/icons/shuaxin_black.png'
            else:
                back_icon_path = 'images/icons/houtui.png'
                refresh_icon_path = 'images/icons/shuaxin.png'
            
            # 更新按钮图标
            if hasattr(self.main_window.ui, 'btnBack'):
                self.main_window.ui.btnBack.setIcon(QIcon(back_icon_path))
            
            if hasattr(self.main_window.ui, 'btnRefresh'):
                self.main_window.ui.btnRefresh.setIcon(QIcon(refresh_icon_path))
                
        except Exception as e:
            pass  # 直接更新按钮图标时出错
    
    def update_buttons_in_widget(self, widget):
        """
        更新指定widget中的按钮图标
        """
        try:
            # 递归查找按钮
            def find_and_update_buttons(obj):
                if hasattr(obj, 'btnBack') and hasattr(obj, 'btnRefresh'):
                    # 如果找到了按钮，尝试调用update_button_icons方法
                    if hasattr(obj, 'update_button_icons'):
                        obj.update_button_icons()
                        return True
                
                # 递归查找子对象
                for child in obj.children():
                    if find_and_update_buttons(child):
                        return True
                return False
            
            find_and_update_buttons(widget)
        except Exception as e:
            pass  # 更新widget中的按钮图标时出错
    
    def toggle_theme(self):
        """切换主题"""
        new_theme = "light" if self.current_theme == "dark" else "dark"
        return self.set_theme(new_theme)
    
    def get_current_theme(self):
        """获取当前主题"""
        return self.current_theme
    
    def is_dark_theme(self):
        """判断是否为暗色主题"""
        return self.current_theme == "dark"
    


class ThemeSwitchButton(QPushButton):
    """主题切换按钮"""
    
    def __init__(self, theme_manager, parent=None):
        """
        初始化主题切换按钮
        
        Args:
            theme_manager: 主题管理器实例
            parent: 父组件
        """
        super().__init__(parent)
        self.theme_manager = theme_manager
        
        # 设置按钮属性
        self.setObjectName("themeSwitchButton")  # 设置对象名以便QSS样式控制
        self.setFixedSize(28, 28)  # 与EN按钮相同大小
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setToolTip("切换主题 (白天/黑夜)")
        
        # 设置初始图标
        self.update_icon()
        
        # 绑定点击事件
        self.clicked.connect(self.toggle_theme)
        
    def update_icon(self):
        """更新按钮图标"""
        if self.theme_manager.is_dark_theme():
            # 暗色主题时显示月亮图标（表示当前是暗色）
            self.setText("🌙")
            self.setToolTip("切换到白天模式")
        else:
            # 亮色主题时显示太阳图标（表示当前是亮色）
            self.setText("☀️")
            self.setToolTip("切换到黑夜模式")
            
        # 按钮样式由QSS文件控制
    
    def toggle_theme(self):
        """切换主题"""
        if self.theme_manager.toggle_theme():
            self.update_icon()


def create_theme_switch_widget(main_window):
    """
    创建主题切换组件
    
    Args:
        main_window: 主窗口实例
        
    Returns:
        tuple: (theme_manager, theme_widget)
    """
    # 创建主题管理器
    theme_manager = ThemeManager(main_window)
    
    # 创建主题切换按钮
    theme_button = ThemeSwitchButton(theme_manager)
    
    # 创建容器组件
    theme_widget = QWidget()
    layout = QHBoxLayout(theme_widget)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.addWidget(theme_button)
    layout.addStretch()  # 添加弹性空间，让按钮靠右
    
    return theme_manager, theme_widget


def add_theme_switch_to_ui(ui_instance, main_window):
    """
    将主题切换功能添加到UI中
    
    Args:
        ui_instance: UI实例
        main_window: 主窗口实例
        
    Returns:
        ThemeManager: 主题管理器实例
    """
    # 创建主题切换组件
    theme_manager, theme_widget = create_theme_switch_widget(main_window)
    
    # 将主题切换按钮添加到主窗口的顶部
    # 这里需要根据你的UI结构调整
    if hasattr(ui_instance, 'centralwidget'):
        # 如果有centralwidget，添加到其布局中
        if hasattr(ui_instance.centralwidget, 'layout'):
            layout = ui_instance.centralwidget.layout()
            if layout:
                # 在布局顶部添加主题切换组件
                layout.insertWidget(0, theme_widget)
    
    return theme_manager


def add_theme_button_to_top_menu(ui_instance, main_window):
    """
    将主题切换按钮添加到顶部按钮区域（EN按钮旁边）
    
    Args:
        ui_instance: UI实例
        main_window: 主窗口实例
        
    Returns:
        ThemeManager: 主题管理器实例
    """
    # 创建主题管理器
    theme_manager = ThemeManager(main_window)
    
    # 创建主题切换按钮
    theme_button = ThemeSwitchButton(theme_manager)
    
    # 添加到rightButtons的布局中
    if hasattr(ui_instance, 'rightButtons'):
        # 获取rightButtons的布局
        layout = ui_instance.rightButtons.layout()
        if layout:
            # 在EN按钮之前添加主题切换按钮
            layout.insertWidget(0, theme_button)
    else:
        pass  # 未找到rightButtons组件
    
    return theme_manager

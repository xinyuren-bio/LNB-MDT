# -*- coding: utf-8 -*-
"""
主题管理模块

提供主题切换功能，支持白天/黑夜主题切换
"""

import os
from PySide6.QtWidgets import QPushButton, QHBoxLayout, QWidget, QApplication
from PySide6.QtCore import Qt, QSettings
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
        
        # 调试信息
        print(f"主题文件路径: {self.theme_path}")
        print(f"主题文件是否存在: {os.path.exists(self.theme_path)}")
        
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
            print(f"主题文件不存在: {theme_file}")
            return False
            
        try:
            with open(theme_file, 'r', encoding='utf-8') as f:
                style_sheet = f.read()
                
            # 应用主题到UI的styleSheet组件
            if hasattr(self.main_window.ui, 'styleSheet'):
                self.main_window.ui.styleSheet.setStyleSheet(style_sheet)
                print(f"主题已切换到: {theme_name}")
            else:
                # 如果没有styleSheet组件，则应用到主窗口
                self.main_window.setStyleSheet(style_sheet)
                print(f"主题已切换到: {theme_name} (应用到主窗口)")
            
            # 强制更新analysis页面的组件样式
            self._update_analysis_components_style(theme_name)
                
            self.current_theme = theme_name
            
            # 保存设置
            self.save_theme_setting()
            
            return True
            
        except Exception as e:
            print(f"加载主题失败: {e}")
            return False
    
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
    
    def _update_analysis_components_style(self, theme_name):
        """更新analysis页面组件的样式"""
        try:
            if not hasattr(self.main_window.ui, 'page_analysis'):
                return
                
            # 定义主题相关的颜色
            if theme_name == "dark":
                bg_color = "rgb(33, 37, 43)"
                text_color = "white"
                border_color = "white"
            else:  # light theme
                bg_color = "rgb(248, 248, 242)"
                text_color = "rgb(40, 42, 54)"
                border_color = "rgb(68, 71, 90)"
            
            # 更新QSpinBox样式
            spinbox_style = f"""
                QSpinBox {{
                    font: 12pt "华文细黑";
                    color: {text_color};
                    border: 1px solid {border_color};
                    background-color: {bg_color};
                }}
                QSpinBox::up-button {{
                    background-color: {bg_color};
                    border: 1px solid {border_color};
                    width: 30px;
                    subcontrol-position: top right;
                }}
                QSpinBox::down-button {{
                    background-color: {bg_color};
                    border: 1px solid {border_color};
                    width: 30px;
                    subcontrol-position: bottom right;
                }}
            """
            
            # 更新QDoubleSpinBox样式
            doublespinbox_style = f"""
                QDoubleSpinBox {{
                    font: 16pt "华文细黑";
                    color: {text_color};
                    border: 1px solid {border_color};
                    background-color: {bg_color};
                }}
                QDoubleSpinBox::up-button {{
                    background-color: {bg_color};
                    border: 1px solid {border_color};
                    width: 30px;
                    subcontrol-position: top right;
                }}
                QDoubleSpinBox::down-button {{
                    background-color: {bg_color};
                    border: 1px solid {border_color};
                    width: 30px;
                    subcontrol-position: bottom right;
                }}
            """
            
            # 更新Label样式
            label_style = f"""
                QLabel {{
                    color: {text_color};
                    font: 12pt "华文细黑";
                }}
            """
            
            # 应用样式到analysis页面的组件
            analysis_page = self.main_window.ui.page_analysis
            
            # 更新QSpinBox组件
            for spinbox in analysis_page.findChildren("QSpinBox"):
                if hasattr(spinbox, 'setStyleSheet'):
                    spinbox.setStyleSheet(spinbox_style)
            
            # 更新QDoubleSpinBox组件
            for doublespinbox in analysis_page.findChildren("QDoubleSpinBox"):
                if hasattr(doublespinbox, 'setStyleSheet'):
                    doublespinbox.setStyleSheet(doublespinbox_style)
            
            # 更新Label组件
            for label in analysis_page.findChildren("QLabel"):
                if hasattr(label, 'setStyleSheet'):
                    label.setStyleSheet(label_style)
            
            # 更新QLineEdit组件
            placeholder_color = "rgb(139, 139, 139)" if theme_name == "dark" else "rgb(50, 50, 50)"
            lineedit_style = f"""
                QLineEdit {{
                    background-color: {bg_color};
                    border-radius: 5px;
                    border: 2px solid {bg_color};
                    padding-left: 10px;
                    color: {text_color};
                    selection-color: rgb(255, 255, 255);
                    selection-background-color: rgb(255, 121, 198);
                }}
                QLineEdit::placeholder {{
                    color: {placeholder_color};
                }}
                QLineEdit:focus {{
                    border: 2px solid {border_color};
                }}
            """
            for lineedit in analysis_page.findChildren("QLineEdit"):
                if hasattr(lineedit, 'setStyleSheet'):
                    lineedit.setStyleSheet(lineedit_style)
            
            # 特别处理analysis_label_first，因为它有内嵌的HTML内容和白色文字样式
            if hasattr(self.main_window.ui, 'analysis_label_first'):
                # 直接修改HTML内容中的颜色，使用更强的样式覆盖
                if theme_name == "dark":
                    html_content = '<html><head/><body><p><span style=" font-size:14pt; font-weight:600; color:white !important;">First</span></p></body></html>'
                else:  # light theme
                    html_content = '<html><head/><body><p><span style=" font-size:14pt; font-weight:600; color:rgb(40, 42, 54) !important;">First</span></p></body></html>'
                
                self.main_window.ui.analysis_label_first.setText(html_content)
                
                # 设置更强的样式表覆盖
                first_label_style = f"""
                    QLabel {{
                        color: {text_color} !important;
                        font: 14pt "华文细黑" !important;
                        font-family: Verdana !important;
                    }}
                    QLabel * {{
                        color: {text_color} !important;
                    }}
                """
                self.main_window.ui.analysis_label_first.setStyleSheet(first_label_style)
                print(f"已特别更新analysis_label_first样式和内容: {theme_name}")
                
                # 强制刷新组件
                self.main_window.ui.analysis_label_first.update()
                self.main_window.ui.analysis_label_first.repaint()
                
                # 延迟再次设置，确保生效
                from PySide6.QtCore import QTimer
                QTimer.singleShot(100, lambda: self.main_window.ui.analysis_label_first.setStyleSheet(first_label_style))
            
            # 更新所有页面的组件样式（不仅仅是analysis页面）
            self._update_all_pages_components_style(theme_name, text_color, bg_color, border_color)
                    
            print(f"已更新analysis页面组件样式: {theme_name}")
            
        except Exception as e:
            print(f"更新analysis组件样式失败: {e}")
    
    def _update_all_pages_components_style(self, theme_name, text_color, bg_color, border_color):
        """更新所有页面的组件样式"""
        try:
            # 更新QSpinBox样式
            spinbox_style = f"""
                QSpinBox {{
                    font: 12pt "华文细黑";
                    color: {text_color};
                    border: 1px solid {border_color};
                    background-color: {bg_color};
                }}
                QSpinBox::up-button {{
                    background-color: {bg_color};
                    border: 1px solid {border_color};
                    width: 30px;
                    subcontrol-position: top right;
                }}
                QSpinBox::down-button {{
                    background-color: {bg_color};
                    border: 1px solid {border_color};
                    width: 30px;
                    subcontrol-position: bottom right;
                }}
            """
            
            # 更新QDoubleSpinBox样式
            doublespinbox_style = f"""
                QDoubleSpinBox {{
                    font: 16pt "华文细黑";
                    color: {text_color};
                    border: 1px solid {border_color};
                    background-color: {bg_color};
                }}
                QDoubleSpinBox::up-button {{
                    background-color: {bg_color};
                    border: 1px solid {border_color};
                    width: 30px;
                    subcontrol-position: top right;
                }}
                QDoubleSpinBox::down-button {{
                    background-color: {bg_color};
                    border: 1px solid {border_color};
                    width: 30px;
                    subcontrol-position: bottom right;
                }}
            """
            
            # 更新Label样式
            label_style = f"""
                QLabel {{
                    color: {text_color};
                    font: 12pt "华文细黑";
                }}
            """
            
            # 更新所有页面的组件
            for page_name in ['page_home', 'page_analysis', 'page_generation', 'page_figure', 'page_data_process']:
                if hasattr(self.main_window.ui, page_name):
                    page = getattr(self.main_window.ui, page_name)
                    
                    # 更新QSpinBox组件
                    for spinbox in page.findChildren("QSpinBox"):
                        if hasattr(spinbox, 'setStyleSheet'):
                            spinbox.setStyleSheet(spinbox_style)
                    
                    # 更新QDoubleSpinBox组件
                    for doublespinbox in page.findChildren("QDoubleSpinBox"):
                        if hasattr(doublespinbox, 'setStyleSheet'):
                            doublespinbox.setStyleSheet(doublespinbox_style)
                    
                    # 更新Label组件
                    for label in page.findChildren("QLabel"):
                        if hasattr(label, 'setStyleSheet'):
                            label.setStyleSheet(label_style)
                    
                    # 更新QLineEdit组件
                    placeholder_color = "rgb(139, 139, 139)" if theme_name == "dark" else "rgb(50, 50, 50)"
                    lineedit_style = f"""
                        QLineEdit {{
                            background-color: {bg_color};
                            border-radius: 5px;
                            border: 2px solid {bg_color};
                            padding-left: 10px;
                            color: {text_color};
                            selection-color: rgb(255, 255, 255);
                            selection-background-color: rgb(255, 121, 198);
                        }}
                        QLineEdit::placeholder {{
                            color: {placeholder_color};
                        }}
                        QLineEdit:focus {{
                            border: 2px solid {border_color};
                        }}
                    """
                    for lineedit in page.findChildren("QLineEdit"):
                        if hasattr(lineedit, 'setStyleSheet'):
                            lineedit.setStyleSheet(lineedit_style)
            
            print(f"已更新所有页面组件样式: {theme_name}")
            
        except Exception as e:
            print(f"更新所有页面组件样式失败: {e}")
    


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
            
        # 设置按钮样式 - 与EN按钮保持一致
        self.setStyleSheet("""
            QPushButton {
                background-color: rgba(255, 255, 255, 0);
                border: none;
                border-radius: 5px;
                font-size: 14px;
                color: white;
                font: 14pt "华文细黑";
            }
            QPushButton:hover {
                # background-color: rgb(44, 49, 57);
                # border-style: solid;
                # border-radius: 4px;
            }
            QPushButton:pressed {
                # background-color: rgb(23, 26, 30);
                # border-style: solid;
                # border-radius: 4px;
            }
        """)
    
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
            print("主题切换按钮已添加到顶部按钮区域")
    else:
        print("未找到rightButtons组件")
    
    return theme_manager

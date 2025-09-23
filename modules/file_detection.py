# -*- coding: utf-8 -*-
"""
CSV文件类型检测和动态TabWidget管理模块

功能：
1. 检测CSV文件类型（lipids、bubble等）
2. 根据文件类型动态创建相应的TabWidget界面
3. 管理TabWidget的切换和布局

支持的文件类型：
- lipids: 显示Line、Bar、Scatter三个tab
- bubble: 显示Line、Bar两个tab
"""

import os
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                               QPushButton, QTextEdit, QFileDialog, QMessageBox, QTabWidget,
                               QGridLayout, QLineEdit, QDoubleSpinBox, QComboBox, QRadioButton)
from PySide6.QtCore import Qt, QSize
from PySide6.QtGui import QFont, QIcon, QCursor

class FileTypeDetector:
    """文件类型检测器"""
    
    @staticmethod
    def detect_file_type(file_path):
        """检测CSV文件类型（lipids或bubble）"""
        print(f"DEBUG: FileTypeDetector.detect_file_type called with: {file_path}")
        
        if not os.path.exists(file_path):
            print(f"DEBUG: File does not exist: {file_path}")
            return "unknown"
        
        # 获取文件扩展名
        _, ext = os.path.splitext(file_path.lower())
        print(f"DEBUG: File extension: {ext}")
        
        # 只支持CSV文件
        if ext != '.csv':
            print(f"DEBUG: Unsupported file type: {ext}")
            return "unsupported"
        
        # 读取CSV文件的前几行来检测类型
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                print(f"DEBUG: Read {len(lines)} lines from file")
                
                # 检查前几行（最多10行）的TYPE信息
                for i in range(min(10, len(lines))):
                    line = lines[i].strip()
                    print(f"DEBUG: Line {i}: {line}")
                    if line.startswith('# TYPE:'):
                        type_info = line.split('# TYPE:')[1].strip()
                        print(f"DEBUG: Found TYPE line: '{type_info}'")
                        
                        if type_info.lower() == 'lipids':
                            print(f"DEBUG: Returning 'lipids'")
                            return 'lipids'
                        elif type_info.lower() == 'bubble':
                            print(f"DEBUG: Returning 'bubble'")
                            return 'bubble'
                        elif type_info.lower() == 'density with time':
                            print(f"DEBUG: Returning 'density_time'")
                            return 'density_time'
                        elif type_info.lower() == 'density with radius':
                            print(f"DEBUG: Returning 'density_radius'")
                            return 'density_radius'
                        else:
                            print(f"DEBUG: Unknown CSV type: '{type_info}', returning 'unknown_csv_type'")
                            return 'unknown_csv_type'
                
                # 如果没有找到TYPE信息，返回未知格式
                print(f"DEBUG: No TYPE line found, returning 'unknown_csv_format'")
                return 'unknown_csv_format'
                    
        except Exception as e:
            print(f"DEBUG: Exception reading CSV file: {e}")
            pass  # 读取CSV文件时出错

class DynamicTabManager:
    """动态Tab Widget管理器"""
    
    def __init__(self, ui_instance):
        self.ui_instance = ui_instance
        self.current_tab_widget = None
        self.original_tab_widget = None
    
    def replace_tab_widget(self, file_type, file_path):
        """根据文件类型创建新的TabWidget并替换当前的TabWidget"""
        
        # 保存原始的TabWidget（如果还没有保存的话）
        if self.original_tab_widget is None:
            # 由于tabWidget_lipids已被删除，我们需要创建一个默认的TabWidget作为参考
            # 或者直接使用None，让后续代码处理
            self.original_tab_widget = None
        
        # 创建新的TabWidget
        new_tab_widget = QTabWidget()
        
        # 设置新TabWidget的样式（使用默认样式）
        new_tab_widget.setStyleSheet("""
            QTabBar::tab {
                background: lightgray;
                border: 2px solid #C4C4C3;
                border-bottom-color: #C4C4C3;
                border-top-left-radius: 5px;
                border-top-right-radius: 5px;
                min-width: 16ex;
                padding: 2px;
                font: 18pt "华文细黑";
                color: black;
            }
            QTabBar::tab:selected {
                background: lightblue;
            }
            QTabBar::tab:hover {
                background: pink;
            }
        """)
        new_tab_widget.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        
        # 根据文件类型创建不同的内容
        if file_type == 'lipids':
            # 对于lipids类型，创建美化的TabWidget
            return self._create_lipids_tab_widget(file_path)
        elif file_type == 'bubble':
            # 对于bubble类型，创建新的tabWidget，复制Line和Bar
            return self._create_bubble_tab_widget(file_path)
        elif file_type == 'density_time':
            # 对于density_time类型，创建类似bubble的TabWidget
            return self._create_density_time_tab_widget(file_path)
        elif file_type == 'density_radius':
            # 对于density_radius类型，创建包含Heatmap的TabWidget
            return self._create_density_radius_tab_widget(file_path)
        else:
            # 创建错误信息tab
            content_widget = QWidget()
            layout = QVBoxLayout(content_widget)
            self._create_error_tab(layout, file_path, file_type)
            
            # 添加内容到新的TabWidget
            tab_name = self._get_tab_name(file_type)
            new_tab_widget.addTab(content_widget, tab_name)
            
            # 替换父widget中的TabWidget
            self._replace_widget_in_layout(new_tab_widget)
            
            # 保存当前TabWidget的引用
            self.current_tab_widget = new_tab_widget
            
            return new_tab_widget
    
    def restore_original_tab_widget(self):
        """恢复原始的TabWidget"""
        if self.original_tab_widget is not None:
            self._replace_widget_in_layout(self.original_tab_widget)
            self.current_tab_widget = None
    
    def _replace_widget_in_layout(self, new_widget):
        """在布局中替换TabWidget，确保RUN按钮位置不变"""
        # 获取widget_2的布局
        parent_layout = self.ui_instance.widget_2.layout()
        
        # 找到当前TabWidget的位置
        current_widget = self.current_tab_widget if self.current_tab_widget else self.original_tab_widget
        
        if current_widget and parent_layout:
            # 获取RUN按钮的引用
            run_button = self.ui_instance.btn_figure_run
            
            # 从布局中移除当前TabWidget
            parent_layout.removeWidget(current_widget)
            
            # 隐藏当前TabWidget
            current_widget.hide()
            
            # 添加新的TabWidget到布局中（在TabWidget的位置，RUN按钮之前）
            # 找到RUN按钮在布局中的位置
            run_button_index = parent_layout.indexOf(run_button)
            if run_button_index >= 0:
                # 在RUN按钮之前插入新的TabWidget
                parent_layout.insertWidget(run_button_index, new_widget)
            else:
                # 如果找不到RUN按钮，直接添加到末尾
                parent_layout.addWidget(new_widget)
            
            # 显示新的TabWidget
            new_widget.show()
        elif parent_layout:
            # 如果没有当前TabWidget，直接添加到布局中
            # 获取RUN按钮的引用
            run_button = self.ui_instance.btn_figure_run
            
            # 找到RUN按钮在布局中的位置
            run_button_index = parent_layout.indexOf(run_button)
            if run_button_index >= 0:
                # 在RUN按钮之前插入新的TabWidget
                parent_layout.insertWidget(run_button_index, new_widget)
            else:
                # 如果找不到RUN按钮，直接添加到末尾
                parent_layout.addWidget(new_widget)
            
            # 显示新的TabWidget
            new_widget.show()
            
            # 确保布局更新
            parent_layout.update()
    
    def _get_tab_name(self, file_type):
        """获取tab显示名称"""
        name_mapping = {
            'lipids': '🧬 Lipids分析',
            'bubble': '🫧 Bubble分析',
            'density_time': '📊 Density Time分析',
            'density_radius': '📈 Density Radius分析',
            'unsupported': '❌ 不支持的文件类型',
            'unknown_csv_type': '❓ 未知CSV类型',
            'unknown_csv_format': '❓ 未知CSV格式',
            'invalid_csv_format': '❌ 无效CSV格式',
            'csv_read_error': '❌ CSV读取错误',
            'unknown': '❓ 未知文件'
        }
        return name_mapping.get(file_type, '❓ 未知文件')
    
    
    def _create_bubble_tab_widget(self, file_path):
        """创建Bubble类型的TabWidget，完全仿照原有UI设计"""
        # 创建新的TabWidget
        bubble_tab_widget = QTabWidget()
        
        # 设置样式（使用默认样式）
        bubble_tab_widget.setStyleSheet("""
            QTabBar::tab {
                background: lightgray;
                border: 2px solid #C4C4C3;
                border-bottom-color: #C4C4C3;
                border-top-left-radius: 5px;
                border-top-right-radius: 5px;
                min-width: 16ex;
                padding: 2px;
                font: 18pt "华文细黑";
                color: black;
            }
            QTabBar::tab:selected {
                background: lightblue;
            }
            QTabBar::tab:hover {
                background: pink;
            }
        """)
        bubble_tab_widget.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        
        # 创建Line tab - 完全仿照原有设计
        line_tab = self._create_bubble_line_tab()
        bubble_tab_widget.addTab(line_tab, "Line")
        
        # 创建Bar tab - 完全仿照原有设计
        bar_tab = self._create_bubble_bar_tab()
        bubble_tab_widget.addTab(bar_tab, "Bar")
        
        # 替换当前的TabWidget
        self._replace_widget_in_layout(bubble_tab_widget)
        self.current_tab_widget = bubble_tab_widget
        
        # 设置按钮绑定
        self.setup_button_bindings_for_new_tab(bubble_tab_widget, 'bubble', file_path)
        
        return bubble_tab_widget
    
    def _create_error_tab(self, layout, file_path, file_type):
        """创建错误信息tab"""
        error_label = QLabel(f"❌ 文件处理错误")
        error_label.setWordWrap(True)
        error_label.setStyleSheet("""
            QLabel {
                font-size: 14px; 
                padding: 10px; 
                background-color: #ffebee; 
                border-radius: 5px;
                # color: #C62828;
                font-weight: bold;
            }
        """)
        layout.addWidget(error_label)
        
        # 显示错误信息
        error_text = QTextEdit()
        error_text.setReadOnly(True)
        error_text.setStyleSheet("""
            QTextEdit {
                background-color: #f5f5f5;
                border: 1px solid #ddd;
                border-radius: 5px;
                padding: 10px;
                font-size: 14px;
            }
        """)
        
        error_messages = {
            'unsupported': f"不支持的文件类型。\n文件: {file_path}\n\n只支持CSV文件。",
            'unknown_csv_type': f"未知的CSV类型。\n文件: {file_path}\n\n支持的CSV类型: lipids, bubble, density_time, density_radius",
            'unknown_csv_format': f"未知的CSV格式。\n文件: {file_path}\n\nCSV文件第4行应包含: # TYPE:lipids 或 # TYPE:bubble 或 # TYPE:Density With Time 或 # TYPE:Density With Radius",
            'invalid_csv_format': f"无效的CSV格式。\n文件: {file_path}\n\nCSV文件至少需要4行。",
            'csv_read_error': f"CSV文件读取错误。\n文件: {file_path}\n\n请检查文件是否存在且可读。"
        }
        
        error_text.setText(error_messages.get(file_type, f"未知错误: {file_type}"))
        layout.addWidget(error_text)
    
    
    def _create_bubble_line_tab(self):
        """创建Bubble Line tab - 完全仿照原有UI设计"""
        # 创建Line tab widget
        line_tab = QWidget()
        line_tab.setObjectName("bubble_line_tab")
        
        # 创建GridLayout - 仿照原有的gridLayout_9
        grid_layout = QGridLayout(line_tab)
        grid_layout.setObjectName("bubble_line_gridLayout")
        
        # 创建所有控件 - 完全仿照原有设计
        
        # Axis Tick Size
        axis_tick_label = QLabel(line_tab)
        axis_tick_label.setObjectName("bubble_line_label_axis_tick")
        axis_tick_label.setStyleSheet("font: 14pt \"华文细黑\";")
        axis_tick_label.setText("Axis Tick Size")
        grid_layout.addWidget(axis_tick_label, 0, 0, 1, 1)
        
        axis_tick_spin = QDoubleSpinBox(line_tab)
        axis_tick_spin.setObjectName("bubble_line_spin_axis_tick")
        axis_tick_spin.setMinimum(0.0)
        axis_tick_spin.setMaximum(100.0)
        axis_tick_spin.setValue(12.0)
        grid_layout.addWidget(axis_tick_spin, 0, 1, 1, 1)
        
        # Axis Title Size
        axis_title_label = QLabel(line_tab)
        axis_title_label.setObjectName("bubble_line_label_axis_title")
        axis_title_label.setMinimumSize(QSize(50, 0))
        axis_title_label.setStyleSheet("font: 14pt \"华文细黑\";")
        axis_title_label.setText("Axis Title Size")
        grid_layout.addWidget(axis_title_label, 1, 0, 1, 1)
        
        axis_title_spin = QDoubleSpinBox(line_tab)
        axis_title_spin.setObjectName("bubble_line_spin_axis_title")
        axis_title_spin.setMinimum(0.0)
        axis_title_spin.setMaximum(100.0)
        axis_title_spin.setValue(16.0)
        grid_layout.addWidget(axis_title_spin, 1, 1, 1, 1)
        
        # Legend Size
        legend_label = QLabel(line_tab)
        legend_label.setObjectName("bubble_line_label_legend")
        legend_label.setStyleSheet("font: 14pt \"华文细黑\";")
        legend_label.setText("Legend Size(0=None)")
        grid_layout.addWidget(legend_label, 2, 0, 1, 1)
        
        legend_spin = QDoubleSpinBox(line_tab)
        legend_spin.setObjectName("bubble_line_spin_legend")
        legend_spin.setMinimum(0.0)
        legend_spin.setMaximum(100.0)
        legend_spin.setValue(14.0)
        grid_layout.addWidget(legend_spin, 2, 1, 1, 1)
        
        # X-Title
        x_title_label = QLabel(line_tab)
        x_title_label.setObjectName("bubble_line_label_x")
        x_title_label.setStyleSheet("font: 14pt \"华文细黑\";")
        x_title_label.setText("X-Title")
        grid_layout.addWidget(x_title_label, 3, 0, 1, 1)
        
        x_title_edit = QLineEdit(line_tab)
        x_title_edit.setObjectName("bubble_line_edit_x")
        x_title_edit.setPlaceholderText("Default if none")
        grid_layout.addWidget(x_title_edit, 3, 1, 1, 1)
        
        # Y-Title
        y_title_label = QLabel(line_tab)
        y_title_label.setObjectName("bubble_line_label_y")
        y_title_label.setStyleSheet("font: 14pt \"华文细黑\";")
        y_title_label.setText("Y-Title")
        grid_layout.addWidget(y_title_label, 4, 0, 1, 1)
        
        y_title_edit = QLineEdit(line_tab)
        y_title_edit.setObjectName("bubble_line_edit_y")
        y_title_edit.setPlaceholderText("Default if none")
        grid_layout.addWidget(y_title_edit, 4, 1, 1, 1)
        
        # X-Range
        x_range_label = QLabel(line_tab)
        x_range_label.setObjectName("bubble_line_label_x_range")
        x_range_label.setStyleSheet("font: 14pt \"华文细黑\";")
        x_range_label.setText("X-Range")
        grid_layout.addWidget(x_range_label, 5, 0, 1, 1)
        
        x_min_spin = QDoubleSpinBox(line_tab)
        x_min_spin.setObjectName("bubble_line_spin_x_min")
        x_min_spin.setMinimum(-1000000000.0)
        x_min_spin.setMaximum(1000000000.0)
        grid_layout.addWidget(x_min_spin, 5, 1, 1, 1)
        
        x_max_spin = QDoubleSpinBox(line_tab)
        x_max_spin.setObjectName("bubble_line_spin_x_max")
        x_max_spin.setMinimum(-1000000000.0)
        x_max_spin.setMaximum(1000000000.0)
        grid_layout.addWidget(x_max_spin, 5, 2, 1, 1)
        
        # Y-Range
        y_range_label = QLabel(line_tab)
        y_range_label.setObjectName("bubble_line_label_y_range")
        y_range_label.setStyleSheet("font: 14pt \"华文细黑\";")
        y_range_label.setText("Y-Range")
        grid_layout.addWidget(y_range_label, 6, 0, 1, 1)
        
        y_min_spin = QDoubleSpinBox(line_tab)
        y_min_spin.setObjectName("bubble_line_spin_y_min")
        y_min_spin.setMinimum(-1000000000.0)
        y_min_spin.setMaximum(1000000000.0)
        grid_layout.addWidget(y_min_spin, 6, 1, 1, 1)
        
        y_max_spin = QDoubleSpinBox(line_tab)
        y_max_spin.setObjectName("bubble_line_spin_y_max")
        y_max_spin.setMinimum(-1000000000.0)
        y_max_spin.setMaximum(1000000000.0)
        grid_layout.addWidget(y_max_spin, 6, 2, 1, 1)
        
        # Marker Size
        marker_label = QLabel(line_tab)
        marker_label.setObjectName("bubble_line_label_marker")
        marker_label.setStyleSheet("font: 14pt \"华文细黑\";")
        marker_label.setText("Marker Size")
        grid_layout.addWidget(marker_label, 7, 0, 1, 1)
        
        marker_spin = QDoubleSpinBox(line_tab)
        marker_spin.setObjectName("bubble_line_spin_marker")
        marker_spin.setMinimum(0.0)
        marker_spin.setMaximum(100.0)
        marker_spin.setValue(0)
        grid_layout.addWidget(marker_spin, 7, 1, 1, 1)
        
        # Color
        color_label = QLabel(line_tab)
        color_label.setObjectName("bubble_line_label_color")
        color_label.setStyleSheet("font: 14pt \"华文细黑\";")
        color_label.setText("Color")
        grid_layout.addWidget(color_label, 8, 0, 1, 1)
        
        color_btn = QPushButton(line_tab)
        color_btn.setObjectName("bubble_line_btn_color")
        color_btn.setEnabled(True)
        color_btn.setMinimumSize(QSize(0, 38))
        color_btn.setMaximumSize(QSize(16000000, 16777215))
        color_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        color_btn.setStyleSheet("background-color: rgb(189,147,249);\ncolor:white;\nfont: 16pt \"华文细黑\";")
        color_btn.setText("Select Color")
        grid_layout.addWidget(color_btn, 8, 1, 1, 1)
        
        return line_tab
    
    def _create_bubble_bar_tab(self):
        """创建Bubble Bar tab - 完全仿照原有UI设计"""
        # 创建Bar tab widget
        bar_tab = QWidget()
        bar_tab.setObjectName("bubble_bar_tab")
        
        # 创建GridLayout - 仿照原有的gridLayout_11
        grid_layout = QGridLayout(bar_tab)
        grid_layout.setObjectName("bubble_bar_gridLayout")
        
        # 创建所有控件 - 完全仿照原有设计
        
        # Axis Tick Size
        axis_tick_label = QLabel(bar_tab)
        axis_tick_label.setObjectName("bubble_bar_label_axis_tick")
        axis_tick_label.setStyleSheet("font: 16pt \"华文细黑\";")
        axis_tick_label.setText("Axis Tick Size")
        grid_layout.addWidget(axis_tick_label, 0, 0, 1, 1)
        
        axis_tick_spin = QDoubleSpinBox(bar_tab)
        axis_tick_spin.setObjectName("bubble_bar_spin_axis_tick")
        axis_tick_spin.setMinimum(0.0)
        axis_tick_spin.setMaximum(100.0)
        axis_tick_spin.setValue(12.0)
        grid_layout.addWidget(axis_tick_spin, 0, 1, 1, 1)
        
        # Axis Title Size
        axis_title_label = QLabel(bar_tab)
        axis_title_label.setObjectName("bubble_bar_label_axis_title")
        axis_title_label.setStyleSheet("font: 16pt \"华文细黑\";")
        axis_title_label.setText("Axis Title Size")
        grid_layout.addWidget(axis_title_label, 1, 0, 1, 1)
        
        axis_title_spin = QDoubleSpinBox(bar_tab)
        axis_title_spin.setObjectName("bubble_bar_spin_axis_title")
        axis_title_spin.setMinimum(0.0)
        axis_title_spin.setMaximum(100.0)
        axis_title_spin.setValue(16.0)
        grid_layout.addWidget(axis_title_spin, 1, 1, 1, 1)
        
        # X-Title
        x_title_label = QLabel(bar_tab)
        x_title_label.setObjectName("bubble_bar_label_x")
        x_title_label.setStyleSheet("font: 16pt \"华文细黑\";")
        x_title_label.setText("X-Title")
        grid_layout.addWidget(x_title_label, 2, 0, 1, 1)
        
        x_title_edit = QLineEdit(bar_tab)
        x_title_edit.setObjectName("bubble_bar_edit_x")
        x_title_edit.setPlaceholderText("Default if none")
        grid_layout.addWidget(x_title_edit, 2, 1, 1, 1)
        
        # Y-Title
        y_title_label = QLabel(bar_tab)
        y_title_label.setObjectName("bubble_bar_label_y")
        y_title_label.setStyleSheet("font: 16pt \"华文细黑\";")
        y_title_label.setText("Y-Title")
        grid_layout.addWidget(y_title_label, 3, 0, 1, 1)
        
        y_title_edit = QLineEdit(bar_tab)
        y_title_edit.setObjectName("bubble_bar_edit_y")
        y_title_edit.setPlaceholderText("Default if none")
        grid_layout.addWidget(y_title_edit, 3, 1, 1, 1)
        
        # Y-Range
        y_range_label = QLabel(bar_tab)
        y_range_label.setObjectName("bubble_bar_label_y_range")
        y_range_label.setStyleSheet("font: 16pt \"华文细黑\";")
        y_range_label.setText("Y-Range")
        grid_layout.addWidget(y_range_label, 4, 0, 1, 1)
        
        y_min_spin = QDoubleSpinBox(bar_tab)
        y_min_spin.setObjectName("bubble_bar_spin_y_min")
        y_min_spin.setMinimum(-1000000000.0)
        y_min_spin.setMaximum(1000000000.0)
        grid_layout.addWidget(y_min_spin, 4, 1, 1, 1)
        
        y_max_spin = QDoubleSpinBox(bar_tab)
        y_max_spin.setObjectName("bubble_bar_spin_y_max")
        y_max_spin.setMinimum(-1000000000.0)
        y_max_spin.setMaximum(1000000000.0)
        grid_layout.addWidget(y_max_spin, 4, 2, 1, 1)
        
        # Trend Line
        trend_label = QLabel(bar_tab)
        trend_label.setObjectName("bubble_bar_label_trend")
        trend_label.setStyleSheet("font: 16pt \"华文细黑\";")
        trend_label.setText("Trend Line")
        grid_layout.addWidget(trend_label, 5, 0, 1, 1)
        
        trend_btn = QPushButton(bar_tab)
        trend_btn.setObjectName("bubble_bar_btn_trend")
        trend_btn.setEnabled(True)
        trend_btn.setMinimumSize(QSize(0, 38))
        trend_btn.setMaximumSize(QSize(16000000, 16777215))
        trend_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        trend_btn.setStyleSheet("background-color: rgb(189,147,249);\ncolor:white;\nfont: 16pt \"华文细黑\";")
        trend_btn.setText("Select Color")
        grid_layout.addWidget(trend_btn, 5, 2, 1, 1)
        
        # Bar Value
        bar_value_label = QLabel(bar_tab)
        bar_value_label.setObjectName("bubble_bar_label_bar")
        bar_value_label.setStyleSheet("font: 16pt \"华文细黑\";")
        bar_value_label.setText("Bar Value")
        grid_layout.addWidget(bar_value_label, 6, 0, 1, 1)
        
        bar_value_spin = QDoubleSpinBox(bar_tab)
        bar_value_spin.setObjectName("bubble_bar_spin_bar")
        bar_value_spin.setMinimum(0.0)
        bar_value_spin.setMaximum(100.0)
        bar_value_spin.setValue(0.8)
        grid_layout.addWidget(bar_value_spin, 6, 1, 1, 1)
        
        # Color
        color_label = QLabel(bar_tab)
        color_label.setObjectName("bubble_bar_label_color")
        color_label.setStyleSheet("font: 16pt \"华文细黑\";")
        color_label.setText("Color")
        grid_layout.addWidget(color_label, 7, 0, 1, 1)
        
        color_btn = QPushButton(bar_tab)
        color_btn.setObjectName("bubble_bar_btn_color")
        color_btn.setEnabled(True)
        color_btn.setMinimumSize(QSize(0, 38))
        color_btn.setMaximumSize(QSize(16000000, 16777215))
        color_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        color_btn.setStyleSheet("background-color: rgb(189,147,249);\ncolor:white;\nfont: 16pt \"华文细黑\";")
        color_btn.setText("Select Color")
        grid_layout.addWidget(color_btn, 7, 1, 1, 1)
        
        # Error Bar
        error_label = QLabel(bar_tab)
        error_label.setObjectName("bubble_bar_label_error")
        error_label.setStyleSheet("font: 16pt \"华文细黑\";")
        error_label.setText("Error Bar")
        grid_layout.addWidget(error_label, 8, 0, 1, 1)
        
        error_yes_radio = QRadioButton(bar_tab)
        error_yes_radio.setObjectName("bubble_bar_radio_error_yes")
        error_yes_radio.setText("Yes")
        grid_layout.addWidget(error_yes_radio, 8, 1, 1, 1)
        
        error_no_radio = QRadioButton(bar_tab)
        error_no_radio.setObjectName("bubble_bar_radio_error_no")
        error_no_radio.setText("No")
        error_no_radio.setChecked(True)
        grid_layout.addWidget(error_no_radio, 8, 2, 1, 1)
        
        return bar_tab
    
    def _create_lipids_map_tab(self):
        """创建Lipids Map tab - 包含Map图相关控件"""
        # 创建Map tab widget
        map_tab = QWidget()
        map_tab.setObjectName("lipids_map_tab")
        
        # 创建GridLayout
        grid_layout = QGridLayout(map_tab)
        grid_layout.setObjectName("lipids_map_gridLayout")
        
        # Color Bar
        color_bar_label = QLabel(map_tab)
        color_bar_label.setObjectName("lipids_map_label_color_bar")
        color_bar_label.setStyleSheet("font: 16pt \"华文细黑\";")
        color_bar_label.setText("Color Bar")
        grid_layout.addWidget(color_bar_label, 0, 0, 1, 1)
        
        color_bar_yes_radio = QRadioButton(map_tab)
        color_bar_yes_radio.setObjectName("lipids_map_radio_color_bar_yes")
        color_bar_yes_radio.setText("Yes")
        color_bar_yes_radio.setChecked(True)  # 默认选中
        grid_layout.addWidget(color_bar_yes_radio, 0, 1, 1, 1)
        
        color_bar_no_radio = QRadioButton(map_tab)
        color_bar_no_radio.setObjectName("lipids_map_radio_color_bar_no")
        color_bar_no_radio.setText("No")
        grid_layout.addWidget(color_bar_no_radio, 0, 2, 1, 1)
        
        # X-Title
        x_title_label = QLabel(map_tab)
        x_title_label.setObjectName("lipids_map_label_x")
        x_title_label.setStyleSheet("font: 16pt \"华文细黑\";")
        x_title_label.setText("X-Title")
        grid_layout.addWidget(x_title_label, 1, 0, 1, 1)
        
        x_title_edit = QLineEdit(map_tab)
        x_title_edit.setObjectName("lipids_map_edit_x")
        x_title_edit.setPlaceholderText("Default if none")
        grid_layout.addWidget(x_title_edit, 1, 1, 1, 2)
        
        # Y-Title
        y_title_label = QLabel(map_tab)
        y_title_label.setObjectName("lipids_map_label_y")
        y_title_label.setStyleSheet("font: 16pt \"华文细黑\";")
        y_title_label.setText("Y-Title")
        grid_layout.addWidget(y_title_label, 2, 0, 1, 1)
        
        y_title_edit = QLineEdit(map_tab)
        y_title_edit.setObjectName("lipids_map_edit_y")
        y_title_edit.setPlaceholderText("Default if none")
        grid_layout.addWidget(y_title_edit, 2, 1, 1, 2)
        
        # Value Range Min
        value_min_label = QLabel(map_tab)
        value_min_label.setObjectName("lipids_map_label_value_min")
        value_min_label.setStyleSheet("font: 16pt \"华文细黑\";")
        value_min_label.setText("Value Rangre")
        grid_layout.addWidget(value_min_label, 3, 0, 1, 1)
        
        value_min_spin = QDoubleSpinBox(map_tab)
        value_min_spin.setObjectName("lipids_map_spin_value_min")
        value_min_spin.setMinimum(-1000000000.0)
        value_min_spin.setMaximum(1000000000.0)
        value_min_spin.setValue(0.0)  # 设置默认值
        grid_layout.addWidget(value_min_spin, 3, 1, 1, 1)
        
        # # Value Range Max
        # value_max_label = QLabel(map_tab)
        # value_max_label.setObjectName("lipids_map_label_value_max")
        # value_max_label.setStyleSheet("font: 16pt \"华文细黑\";")
        # value_max_label.setText("Value Range Max")
        # grid_layout.addWidget(value_max_label, 3, 2, 1, 1)
        
        value_max_spin = QDoubleSpinBox(map_tab)
        value_max_spin.setObjectName("lipids_map_spin_value_max")
        value_max_spin.setMinimum(-1000000000.0)
        value_max_spin.setMaximum(1000000000.0)
        value_max_spin.setValue(1.0)  # 设置默认值
        grid_layout.addWidget(value_max_spin, 3, 2, 1, 1)
        
        # Color Map
        color_map_label = QLabel(map_tab)
        color_map_label.setObjectName("lipids_map_label_color_map")
        color_map_label.setStyleSheet("font: 16pt \"华文细黑\";")
        color_map_label.setText("Color Map")
        grid_layout.addWidget(color_map_label, 4, 0, 1, 1)
        
        color_map_combo = QComboBox(map_tab)
        color_map_combo.setObjectName("lipids_map_combo_color_map")
        color_map_combo.addItems(["viridis", "plasma", "inferno", "magma", "jet", "hot", "cool", "spring", "summer", "autumn", "winter"])
        color_map_combo.setCurrentText("viridis")  # 默认选择viridis
        grid_layout.addWidget(color_map_combo, 4, 1, 1, 2)
        
        return map_tab
    
    def _create_lipids_tab_widget(self, file_path):
        """创建Lipids类型的TabWidget，包含Line、Bar、Scatter、Map四个tab"""
        # 创建新的TabWidget
        lipids_tab_widget = QTabWidget()
        
        # 设置样式（使用默认样式）
        lipids_tab_widget.setStyleSheet("""
            QTabBar::tab {
                background: lightgray;
                border: 2px solid #C4C4C3;
                border-bottom-color: #C4C4C3;
                border-top-left-radius: 5px;
                border-top-right-radius: 5px;
                min-width: 16ex;
                padding: 2px;
                font: 18pt "华文细黑";
                color: black;
            }
            QTabBar::tab:selected {
                background: lightblue;
            }
            QTabBar::tab:hover {
                background: pink;
            }
        """)
        lipids_tab_widget.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        
        # 创建Line tab - 使用美化的设计
        line_tab = self._create_lipids_line_tab()
        lipids_tab_widget.addTab(line_tab, "Line")
        
        # 创建Bar tab - 使用美化的设计
        bar_tab = self._create_lipids_bar_tab()
        lipids_tab_widget.addTab(bar_tab, "Bar")
        
        # 创建Scatter tab - 使用美化的设计
        scatter_tab = self._create_lipids_scatter_tab()
        lipids_tab_widget.addTab(scatter_tab, "Scatter")
        
        # 创建Map tab - 新增的Map图标签页
        map_tab = self._create_lipids_map_tab()
        lipids_tab_widget.addTab(map_tab, "Map")
        
        # 替换当前的TabWidget
        self._replace_widget_in_layout(lipids_tab_widget)
        self.current_tab_widget = lipids_tab_widget
        
        # 设置按钮绑定
        self.setup_button_bindings_for_new_tab(lipids_tab_widget, 'lipids', file_path)
        
        return lipids_tab_widget
    
    def _create_lipids_line_tab(self):
        """创建Lipids Line tab - 完全仿照Bubble的设计"""
        # 创建Line tab widget
        line_tab = QWidget()
        line_tab.setObjectName("lipids_line_tab")
        
        # 创建GridLayout
        grid_layout = QGridLayout(line_tab)
        grid_layout.setObjectName("lipids_line_gridLayout")
        
        # 创建所有控件 - 完全仿照Bubble的设计
        
        # Axis Tick Size
        axis_tick_label = QLabel(line_tab)
        axis_tick_label.setObjectName("lipids_line_label_axis_tick")
        axis_tick_label.setStyleSheet("font: 14pt \"华文细黑\";")
        axis_tick_label.setText("Axis Tick Size")
        grid_layout.addWidget(axis_tick_label, 0, 0, 1, 1)
        
        axis_tick_spin = QDoubleSpinBox(line_tab)
        axis_tick_spin.setObjectName("lipids_line_spin_axis_tick")
        axis_tick_spin.setMinimum(0.0)
        axis_tick_spin.setMaximum(100.0)
        axis_tick_spin.setValue(12.0)
        grid_layout.addWidget(axis_tick_spin, 0, 1, 1, 1)
        
        # Axis Title Size
        axis_title_label = QLabel(line_tab)
        axis_title_label.setObjectName("lipids_line_label_axis_title")
        axis_title_label.setMinimumSize(QSize(50, 0))
        axis_title_label.setStyleSheet("font: 14pt \"华文细黑\";")
        axis_title_label.setText("Axis Title Size")
        grid_layout.addWidget(axis_title_label, 1, 0, 1, 1)
        
        axis_title_spin = QDoubleSpinBox(line_tab)
        axis_title_spin.setObjectName("lipids_line_spin_axis_title")
        axis_title_spin.setMinimum(0.0)
        axis_title_spin.setMaximum(100.0)
        axis_title_spin.setValue(16.0)
        grid_layout.addWidget(axis_title_spin, 1, 1, 1, 1)
        
        # Legend Size
        legend_label = QLabel(line_tab)
        legend_label.setObjectName("lipids_line_label_legend")
        legend_label.setStyleSheet("font: 14pt \"华文细黑\";")
        legend_label.setText("Legend Size(0=None)")
        grid_layout.addWidget(legend_label, 2, 0, 1, 1)
        
        legend_spin = QDoubleSpinBox(line_tab)
        legend_spin.setObjectName("lipids_line_spin_legend")
        legend_spin.setMinimum(0.0)
        legend_spin.setMaximum(100.0)
        legend_spin.setValue(14.0)
        grid_layout.addWidget(legend_spin, 2, 1, 1, 1)
        
        # X-Title
        x_title_label = QLabel(line_tab)
        x_title_label.setObjectName("lipids_line_label_x")
        x_title_label.setStyleSheet("font: 14pt \"华文细黑\";")
        x_title_label.setText("X-Title")
        grid_layout.addWidget(x_title_label, 3, 0, 1, 1)
        
        x_title_edit = QLineEdit(line_tab)
        x_title_edit.setObjectName("lipids_line_edit_x")
        x_title_edit.setPlaceholderText("Default if none")
        grid_layout.addWidget(x_title_edit, 3, 1, 1, 1)
        
        # Y-Title
        y_title_label = QLabel(line_tab)
        y_title_label.setObjectName("lipids_line_label_y")
        y_title_label.setStyleSheet("font: 14pt \"华文细黑\";")
        y_title_label.setText("Y-Title")
        grid_layout.addWidget(y_title_label, 4, 0, 1, 1)
        
        y_title_edit = QLineEdit(line_tab)
        y_title_edit.setObjectName("lipids_line_edit_y")
        y_title_edit.setPlaceholderText("Default if none")
        grid_layout.addWidget(y_title_edit, 4, 1, 1, 1)
        
        # X-Range
        x_range_label = QLabel(line_tab)
        x_range_label.setObjectName("lipids_line_label_x_range")
        x_range_label.setStyleSheet("font: 14pt \"华文细黑\";")
        x_range_label.setText("X-Range")
        grid_layout.addWidget(x_range_label, 5, 0, 1, 1)
        
        x_min_spin = QDoubleSpinBox(line_tab)
        x_min_spin.setObjectName("lipids_line_spin_x_min")
        x_min_spin.setMinimum(-1000000000.0)
        x_min_spin.setMaximum(1000000000.0)
        grid_layout.addWidget(x_min_spin, 5, 1, 1, 1)
        
        x_max_spin = QDoubleSpinBox(line_tab)
        x_max_spin.setObjectName("lipids_line_spin_x_max")
        x_max_spin.setMinimum(-1000000000.0)
        x_max_spin.setMaximum(1000000000.0)
        grid_layout.addWidget(x_max_spin, 5, 2, 1, 1)
        
        # Y-Range
        y_range_label = QLabel(line_tab)
        y_range_label.setObjectName("lipids_line_label_y_range")
        y_range_label.setStyleSheet("font: 14pt \"华文细黑\";")
        y_range_label.setText("Y-Range")
        grid_layout.addWidget(y_range_label, 6, 0, 1, 1)
        
        y_min_spin = QDoubleSpinBox(line_tab)
        y_min_spin.setObjectName("lipids_line_spin_y_min")
        y_min_spin.setMinimum(-1000000000.0)
        y_min_spin.setMaximum(1000000000.0)
        grid_layout.addWidget(y_min_spin, 6, 1, 1, 1)
        
        y_max_spin = QDoubleSpinBox(line_tab)
        y_max_spin.setObjectName("lipids_line_spin_y_max")
        y_max_spin.setMinimum(-1000000000.0)
        y_max_spin.setMaximum(1000000000.0)
        grid_layout.addWidget(y_max_spin, 6, 2, 1, 1)
        
        # Marker Size
        marker_label = QLabel(line_tab)
        marker_label.setObjectName("lipids_line_label_marker")
        marker_label.setStyleSheet("font: 14pt \"华文细黑\";")
        marker_label.setText("Marker Size")
        grid_layout.addWidget(marker_label, 7, 0, 1, 1)
        
        marker_spin = QDoubleSpinBox(line_tab)
        marker_spin.setObjectName("lipids_line_spin_marker")
        marker_spin.setMinimum(0.0)
        marker_spin.setMaximum(100.0)
        marker_spin.setValue(0)
        grid_layout.addWidget(marker_spin, 7, 1, 1, 1)
        
        # Color
        color_label = QLabel(line_tab)
        color_label.setObjectName("lipids_line_label_color")
        color_label.setStyleSheet("font: 14pt \"华文细黑\";")
        color_label.setText("Color")
        grid_layout.addWidget(color_label, 8, 0, 1, 1)
        
        color_btn = QPushButton(line_tab)
        color_btn.setObjectName("lipids_line_btn_color")
        color_btn.setEnabled(True)
        color_btn.setMinimumSize(QSize(0, 38))
        color_btn.setMaximumSize(QSize(16000000, 16777215))
        color_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        color_btn.setStyleSheet("background-color: rgb(189,147,249);\ncolor:white;\nfont: 16pt \"华文细黑\";")
        color_btn.setText("Select Color")
        grid_layout.addWidget(color_btn, 8, 1, 1, 1)
        
        return line_tab
    
    def _create_lipids_bar_tab(self):
        """创建Lipids Bar tab - 完全仿照Bubble的设计"""
        # 创建Bar tab widget
        bar_tab = QWidget()
        bar_tab.setObjectName("lipids_bar_tab")
        
        # 创建GridLayout
        grid_layout = QGridLayout(bar_tab)
        grid_layout.setObjectName("lipids_bar_gridLayout")
        
        # 创建所有控件 - 完全仿照Bubble的设计
        
        # Axis Tick Size
        axis_tick_label = QLabel(bar_tab)
        axis_tick_label.setObjectName("lipids_bar_label_axis_tick")
        axis_tick_label.setStyleSheet("font: 16pt \"华文细黑\";")
        axis_tick_label.setText("Axis Tick Size")
        grid_layout.addWidget(axis_tick_label, 0, 0, 1, 1)
        
        axis_tick_spin = QDoubleSpinBox(bar_tab)
        axis_tick_spin.setObjectName("lipids_bar_spin_axis_tick")
        axis_tick_spin.setMinimum(0.0)
        axis_tick_spin.setMaximum(100.0)
        axis_tick_spin.setValue(12.0)
        grid_layout.addWidget(axis_tick_spin, 0, 1, 1, 1)
        
        # Axis Title Size
        axis_title_label = QLabel(bar_tab)
        axis_title_label.setObjectName("lipids_bar_label_axis_title")
        axis_title_label.setStyleSheet("font: 16pt \"华文细黑\";")
        axis_title_label.setText("Axis Title Size")
        grid_layout.addWidget(axis_title_label, 1, 0, 1, 1)
        
        axis_title_spin = QDoubleSpinBox(bar_tab)
        axis_title_spin.setObjectName("lipids_bar_spin_axis_title")
        axis_title_spin.setMinimum(0.0)
        axis_title_spin.setMaximum(100.0)
        axis_title_spin.setValue(16.0)
        grid_layout.addWidget(axis_title_spin, 1, 1, 1, 1)
        
        # X-Title
        x_title_label = QLabel(bar_tab)
        x_title_label.setObjectName("lipids_bar_label_x")
        x_title_label.setStyleSheet("font: 16pt \"华文细黑\";")
        x_title_label.setText("X-Title")
        grid_layout.addWidget(x_title_label, 2, 0, 1, 1)
        
        x_title_edit = QLineEdit(bar_tab)
        x_title_edit.setObjectName("lipids_bar_edit_x")
        x_title_edit.setPlaceholderText("Default if none")
        grid_layout.addWidget(x_title_edit, 2, 1, 1, 1)
        
        # Y-Title
        y_title_label = QLabel(bar_tab)
        y_title_label.setObjectName("lipids_bar_label_y")
        y_title_label.setStyleSheet("font: 16pt \"华文细黑\";")
        y_title_label.setText("Y-Title")
        grid_layout.addWidget(y_title_label, 3, 0, 1, 1)
        
        y_title_edit = QLineEdit(bar_tab)
        y_title_edit.setObjectName("lipids_bar_edit_y")
        y_title_edit.setPlaceholderText("Default if none")
        grid_layout.addWidget(y_title_edit, 3, 1, 1, 1)
        
        # Y-Range
        y_range_label = QLabel(bar_tab)
        y_range_label.setObjectName("lipids_bar_label_y_range")
        y_range_label.setStyleSheet("font: 16pt \"华文细黑\";")
        y_range_label.setText("Y-Range")
        grid_layout.addWidget(y_range_label, 4, 0, 1, 1)
        
        y_min_spin = QDoubleSpinBox(bar_tab)
        y_min_spin.setObjectName("lipids_bar_spin_y_min")
        y_min_spin.setMinimum(-1000000000.0)
        y_min_spin.setMaximum(1000000000.0)
        y_min_spin.setValue(0.0)  # 设置默认值
        grid_layout.addWidget(y_min_spin, 4, 1, 1, 1)
        
        y_max_spin = QDoubleSpinBox(bar_tab)
        y_max_spin.setObjectName("lipids_bar_spin_y_max")
        y_max_spin.setMinimum(-1000000000.0)
        y_max_spin.setMaximum(1000000000.0)
        y_max_spin.setValue(1.0)  # 设置默认值
        grid_layout.addWidget(y_max_spin, 4, 2, 1, 1)
        
        # Trend Line
        trend_label = QLabel(bar_tab)
        trend_label.setObjectName("lipids_bar_label_trend")
        trend_label.setStyleSheet("font: 16pt \"华文细黑\";")
        trend_label.setText("Trend Line")
        grid_layout.addWidget(trend_label, 5, 0, 1, 1)
        
        trend_btn = QPushButton(bar_tab)
        trend_btn.setObjectName("lipids_bar_btn_trend")
        trend_btn.setEnabled(True)
        trend_btn.setMinimumSize(QSize(0, 38))
        trend_btn.setMaximumSize(QSize(16000000, 16777215))
        trend_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        trend_btn.setStyleSheet("background-color: rgb(189,147,249);\ncolor:white;\nfont: 16pt \"华文细黑\";")
        trend_btn.setText("Select Color")
        grid_layout.addWidget(trend_btn, 5, 2, 1, 1)
        
        # Bar Value
        bar_value_label = QLabel(bar_tab)
        bar_value_label.setObjectName("lipids_bar_label_bar")
        bar_value_label.setStyleSheet("font: 16pt \"华文细黑\";")
        bar_value_label.setText("Bar Value")
        grid_layout.addWidget(bar_value_label, 6, 0, 1, 1)
        
        bar_value_spin = QDoubleSpinBox(bar_tab)
        bar_value_spin.setObjectName("lipids_bar_spin_bar")
        bar_value_spin.setMinimum(0.0)
        bar_value_spin.setMaximum(100.0)
        bar_value_spin.setValue(0.8)
        grid_layout.addWidget(bar_value_spin, 6, 1, 1, 1)
        
        # Color
        color_label = QLabel(bar_tab)
        color_label.setObjectName("lipids_bar_label_color")
        color_label.setStyleSheet("font: 16pt \"华文细黑\";")
        color_label.setText("Color")
        grid_layout.addWidget(color_label, 7, 0, 1, 1)
        
        color_btn = QPushButton(bar_tab)
        color_btn.setObjectName("lipids_bar_btn_color")
        color_btn.setEnabled(True)
        color_btn.setMinimumSize(QSize(0, 38))
        color_btn.setMaximumSize(QSize(16000000, 16777215))
        color_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        color_btn.setStyleSheet("background-color: rgb(189,147,249);\ncolor:white;\nfont: 16pt \"华文细黑\";")
        color_btn.setText("Select Color")
        grid_layout.addWidget(color_btn, 7, 1, 1, 1)
        
        # Error Bar
        error_label = QLabel(bar_tab)
        error_label.setObjectName("lipids_bar_label_error")
        error_label.setStyleSheet("font: 16pt \"华文细黑\";")
        error_label.setText("Error Bar")
        grid_layout.addWidget(error_label, 8, 0, 1, 1)
        
        error_yes_radio = QRadioButton(bar_tab)
        error_yes_radio.setObjectName("lipids_bar_radio_error_yes")
        error_yes_radio.setText("Yes")
        grid_layout.addWidget(error_yes_radio, 8, 1, 1, 1)
        
        error_no_radio = QRadioButton(bar_tab)
        error_no_radio.setObjectName("lipids_bar_radio_error_no")
        error_no_radio.setText("No")
        error_no_radio.setChecked(True)
        grid_layout.addWidget(error_no_radio, 8, 2, 1, 1)
        
        return bar_tab
    
    def _create_lipids_scatter_tab(self):
        """创建Lipids Scatter tab - 完全仿照Bubble的设计"""
        # 创建Scatter tab widget
        scatter_tab = QWidget()
        scatter_tab.setObjectName("lipids_scatter_tab")
        
        # 创建GridLayout
        grid_layout = QGridLayout(scatter_tab)
        grid_layout.setObjectName("lipids_scatter_gridLayout")
        
        # 创建所有控件 - 完全仿照Bubble的设计
        
        # Axis Tick Size
        axis_tick_label = QLabel(scatter_tab)
        axis_tick_label.setObjectName("lipids_scatter_label_axis_tick")
        axis_tick_label.setStyleSheet("font: 16pt \"华文细黑\";")
        axis_tick_label.setText("Axis Tick Size")
        grid_layout.addWidget(axis_tick_label, 0, 0, 1, 1)
        
        axis_tick_spin = QDoubleSpinBox(scatter_tab)
        axis_tick_spin.setObjectName("lipids_scatter_spin_axis_tick")
        axis_tick_spin.setMinimum(0.0)
        axis_tick_spin.setMaximum(100.0)
        axis_tick_spin.setValue(12.0)
        grid_layout.addWidget(axis_tick_spin, 0, 1, 1, 1)
        
        # Axis Title Size
        axis_title_label = QLabel(scatter_tab)
        axis_title_label.setObjectName("lipids_scatter_label_axis_title")
        axis_title_label.setStyleSheet("font: 16pt \"华文细黑\";")
        axis_title_label.setText("Axis Title Size")
        grid_layout.addWidget(axis_title_label, 1, 0, 1, 1)
        
        axis_title_spin = QDoubleSpinBox(scatter_tab)
        axis_title_spin.setObjectName("lipids_scatter_spin_axis_title")
        axis_title_spin.setMinimum(0.0)
        axis_title_spin.setMaximum(100.0)
        axis_title_spin.setValue(16.0)
        grid_layout.addWidget(axis_title_spin, 1, 1, 1, 1)
        
        # X-Title
        x_title_label = QLabel(scatter_tab)
        x_title_label.setObjectName("lipids_scatter_label_x")
        x_title_label.setStyleSheet("font: 16pt \"华文细黑\";")
        x_title_label.setText("X-Title")
        grid_layout.addWidget(x_title_label, 2, 0, 1, 1)
        
        x_title_edit = QLineEdit(scatter_tab)
        x_title_edit.setObjectName("lipids_scatter_edit_x")
        x_title_edit.setPlaceholderText("Default if none")
        grid_layout.addWidget(x_title_edit, 2, 1, 1, 1)
        
        # Y-Title
        y_title_label = QLabel(scatter_tab)
        y_title_label.setObjectName("lipids_scatter_label_y")
        y_title_label.setStyleSheet("font: 16pt \"华文细黑\";")
        y_title_label.setText("Y-Title")
        grid_layout.addWidget(y_title_label, 3, 0, 1, 1)
        
        y_title_edit = QLineEdit(scatter_tab)
        y_title_edit.setObjectName("lipids_scatter_edit_y")
        y_title_edit.setPlaceholderText("Default if none")
        grid_layout.addWidget(y_title_edit, 3, 1, 1, 1)
        
        # Shape Size
        shape_size_label = QLabel(scatter_tab)
        shape_size_label.setObjectName("lipids_scatter_label_shape_size")
        shape_size_label.setStyleSheet("font: 16pt \"华文细黑\";")
        shape_size_label.setText("Shape Size")
        grid_layout.addWidget(shape_size_label, 4, 0, 1, 1)
        
        shape_size_spin = QDoubleSpinBox(scatter_tab)
        shape_size_spin.setObjectName("lipids_scatter_spin_shape_size")
        shape_size_spin.setMinimum(0.0)
        shape_size_spin.setMaximum(100.0)
        shape_size_spin.setValue(50.0)
        grid_layout.addWidget(shape_size_spin, 4, 1, 1, 1)
        
        # Legend Size
        legend_label = QLabel(scatter_tab)
        legend_label.setObjectName("lipids_scatter_label_legend")
        legend_label.setStyleSheet("font: 16pt \"华文细黑\";")
        legend_label.setText("Legend Size")
        grid_layout.addWidget(legend_label, 5, 0, 1, 1)
        
        legend_spin = QDoubleSpinBox(scatter_tab)
        legend_spin.setObjectName("lipids_scatter_spin_legend")
        legend_spin.setMinimum(0.0)
        legend_spin.setMaximum(100.0)
        legend_spin.setValue(14.0)
        grid_layout.addWidget(legend_spin, 5, 1, 1, 1)
        
        # Value Range
        range_label = QLabel(scatter_tab)
        range_label.setObjectName("lipids_scatter_label_range")
        range_label.setStyleSheet("font: 16pt \"华文细黑\";")
        range_label.setText("Value Range")
        grid_layout.addWidget(range_label, 6, 0, 1, 1)
        
        range_min_spin = QDoubleSpinBox(scatter_tab)
        range_min_spin.setObjectName("lipids_scatter_spin_range_min")
        range_min_spin.setMinimum(-1000000000.0)
        range_min_spin.setMaximum(1000000000.0)
        grid_layout.addWidget(range_min_spin, 6, 1, 1, 1)
        
        range_max_spin = QDoubleSpinBox(scatter_tab)
        range_max_spin.setObjectName("lipids_scatter_spin_range_max")
        range_max_spin.setMinimum(-1000000000.0)
        range_max_spin.setMaximum(1000000000.0)
        grid_layout.addWidget(range_max_spin, 6, 2, 1, 1)
        
        # Shape Selection
        shape_label = QLabel(scatter_tab)
        shape_label.setObjectName("lipids_scatter_label_shape")
        shape_label.setStyleSheet("font: 16pt \"华文细黑\";")
        shape_label.setText("Shape")
        grid_layout.addWidget(shape_label, 7, 0, 1, 1)
        
        shape_btn = QPushButton(scatter_tab)
        shape_btn.setObjectName("lipids_scatter_btn_shape")
        shape_btn.setEnabled(True)
        shape_btn.setMinimumSize(QSize(0, 38))
        shape_btn.setMaximumSize(QSize(16000000, 16777215))
        shape_btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        shape_btn.setStyleSheet("background-color: rgb(189,147,249);\ncolor:white;\nfont: 16pt \"华文细黑\";")
        shape_btn.setText("Select Shape")
        grid_layout.addWidget(shape_btn, 7, 1, 1, 1)
        
        return scatter_tab
    
    def bind_dynamic_tab_buttons(self, tab_widget, file_path):
        """
        为动态创建的TabWidget中的按钮绑定功能
        
        Args:
            tab_widget: 动态创建的TabWidget实例
            file_path: 当前文件路径
        """
        # 获取当前活跃的tab索引
        current_index = tab_widget.currentIndex()
        current_tab = tab_widget.widget(current_index)
        
        # 根据文件类型和当前tab类型绑定相应的按钮
        if hasattr(self, 'current_file_type'):
            if self.current_file_type == 'lipids':
                self._bind_lipids_buttons(tab_widget, current_index, file_path)
            elif self.current_file_type == 'bubble':
                self._bind_bubble_buttons(tab_widget, current_index, file_path)
        else:
            pass  # current_file_type属性不存在
    
    def _bind_lipids_buttons(self, tab_widget, tab_index, file_path):
        """绑定Lipids类型TabWidget的按钮功能"""
        current_tab = tab_widget.widget(tab_index)
        
        if tab_index == 0:  # Line tab
            self._bind_line_buttons(current_tab, file_path)
        elif tab_index == 1:  # Bar tab
            self._bind_bar_buttons(current_tab, file_path)
        elif tab_index == 2:  # Scatter tab
            self._bind_scatter_buttons(current_tab, file_path)
    
    def _bind_bubble_buttons(self, tab_widget, tab_index, file_path):
        """绑定Bubble类型TabWidget的按钮功能"""
        current_tab = tab_widget.widget(tab_index)
        
        if tab_index == 0:  # Line tab
            self._bind_line_buttons(current_tab, file_path)
        elif tab_index == 1:  # Bar tab
            self._bind_bar_buttons(current_tab, file_path)
    
    def _bind_line_buttons(self, line_tab, file_path):
        """绑定Line tab中的按钮功能"""
        # 查找颜色选择按钮
        color_btn = line_tab.findChild(QPushButton, "lipids_line_btn_color")
        if color_btn:
            color_btn.clicked.connect(lambda: self._handle_color_selection(line_tab, "line"))
        else:
            pass  # 未找到lipids_line_btn_color按钮
        
        # 查找bubble颜色按钮（如果存在）
        bubble_color_btn = line_tab.findChild(QPushButton, "bubble_line_btn_color")
        if bubble_color_btn:
            bubble_color_btn.clicked.connect(lambda: self._handle_color_selection(line_tab, "line"))
        else:
            pass  # 未找到bubble_line_btn_color按钮
        
        # 调试：列出所有按钮
        all_buttons = line_tab.findChildren(QPushButton)
        pass  # Line tab中找到的所有按钮
    
    def _bind_bar_buttons(self, bar_tab, file_path):
        """绑定Bar tab中的按钮功能"""
        # 查找颜色选择按钮
        color_btn = bar_tab.findChild(QPushButton, "lipids_bar_btn_color")
        if color_btn:
            color_btn.clicked.connect(lambda: self._handle_color_selection(bar_tab, "bar"))
        
        # 查找bubble颜色按钮（如果存在）
        bubble_color_btn = bar_tab.findChild(QPushButton, "bubble_bar_btn_color")
        if bubble_color_btn:
            bubble_color_btn.clicked.connect(lambda: self._handle_color_selection(bar_tab, "bar"))
        
        # 查找趋势线按钮（如果存在）
        trend_btn = bar_tab.findChild(QPushButton, "lipids_bar_btn_trend")
        if trend_btn:
            trend_btn.clicked.connect(lambda: self._handle_trend_selection(bar_tab))
        
        bubble_trend_btn = bar_tab.findChild(QPushButton, "bubble_bar_btn_trend")
        if bubble_trend_btn:
            bubble_trend_btn.clicked.connect(lambda: self._handle_trend_selection(bar_tab))
    
    def _bind_scatter_buttons(self, scatter_tab, file_path):
        """绑定Scatter tab中的按钮功能"""
        # 查找形状选择按钮
        shape_btn = scatter_tab.findChild(QPushButton, "lipids_scatter_btn_shape")
        if shape_btn:
            shape_btn.clicked.connect(lambda: self._handle_shape_selection(scatter_tab))
    
    def _handle_color_selection(self, tab_widget, tab_type):
        """处理颜色选择功能"""
        try:
            # 导入颜色选择功能
            from .Fuctions_Figure import FigurePage
            
            # 获取主UI实例
            main_ui = self.ui_instance
            
            # 确保FigureInfo已初始化
            if not hasattr(main_ui, 'FigureInfo') or main_ui.FigureInfo is None:
                QMessageBox.warning(None, "警告", "请先选择结果文件！")
                return
            
            # 根据文件类型调用相应的颜色选择方法
            if self.current_file_type == 'lipids':
                FigurePage.lipids_colors(main_ui)
            elif self.current_file_type == 'bubble':
                FigurePage.single_color(main_ui)
            else:
                pass  # 未知的文件类型
                
        except Exception as e:
            QMessageBox.critical(None, "错误", f"颜色选择功能出错: {str(e)}")
    
    def _handle_trend_selection(self, tab_widget):
        """处理趋势线选择功能"""
        try:
            from .Fuctions_Figure import FigurePage
            
            main_ui = self.ui_instance
            
            if not hasattr(main_ui, 'FigureInfo') or main_ui.FigureInfo is None:
                QMessageBox.warning(None, "警告", "请先选择结果文件！")
                return
            
            # 调用趋势线颜色选择
            FigurePage.single_color(main_ui)
            
        except Exception as e:
            QMessageBox.critical(None, "错误", f"趋势线选择功能出错: {str(e)}")
    
    def _handle_shape_selection(self, tab_widget):
        """处理形状选择功能"""
        try:
            from .Fuctions_Figure import FigurePage
            
            main_ui = self.ui_instance
            
            if not hasattr(main_ui, 'FigureInfo') or main_ui.FigureInfo is None:
                QMessageBox.warning(None, "警告", "请先选择结果文件！")
                return
            
            # 调用形状选择功能
            FigurePage.figureBtnShape(main_ui)
            
        except Exception as e:
            QMessageBox.critical(None, "错误", f"形状选择功能出错: {str(e)}")
    
    def setup_button_bindings_for_new_tab(self, tab_widget, file_type, file_path):
        """
        为新创建的TabWidget设置按钮绑定
        
        Args:
            tab_widget: 新创建的TabWidget
            file_type: 文件类型 ('lipids' 或 'bubble')
            file_path: 文件路径
        """
        # 保存当前文件类型
        self.current_file_type = file_type
        
        # 绑定tab切换事件
        tab_widget.currentChanged.connect(
            lambda index: self.bind_dynamic_tab_buttons(tab_widget, file_path)
        )
        
        # 标签页切换时不需要额外处理，图表类型会从TabWidget直接获取
        
        # 立即绑定当前活跃tab的按钮
        self.bind_dynamic_tab_buttons(tab_widget, file_path)
    
    def _create_density_time_tab_widget(self, file_path):
        """创建Density Time类型的TabWidget，包含Line和Bar两个tab"""
        # 创建新的TabWidget
        density_time_tab_widget = QTabWidget()
        
        # 设置样式（使用默认样式）
        density_time_tab_widget.setStyleSheet("""
            QTabBar::tab {
                background: lightgray;
                border: 2px solid #C4C4C3;
                border-bottom-color: #C4C4C3;
                border-top-left-radius: 5px;
                border-top-right-radius: 5px;
                min-width: 16ex;
                padding: 2px;
                font: 18pt "华文细黑";
                color: black;
            }
            QTabBar::tab:selected {
                background: lightblue;
            }
            QTabBar::tab:hover {
                background: pink;
            }
        """)
        density_time_tab_widget.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        
        # 创建Line tab
        line_tab = self._create_bubble_line_tab()
        density_time_tab_widget.addTab(line_tab, "📈 Density-Time Line")
        
        # 创建Heatmap tab
        heatmap_tab = self._create_density_heatmap_tab()
        density_time_tab_widget.addTab(heatmap_tab, "🗺️ Density-Time Heatmap")
        
        # 替换当前的TabWidget
        self._replace_widget_in_layout(density_time_tab_widget)
        self.current_tab_widget = density_time_tab_widget
        
        # 设置按钮绑定
        self.setup_button_bindings_for_new_tab(density_time_tab_widget, 'density_time', file_path)
        
        return density_time_tab_widget
    
    def _create_density_radius_tab_widget(self, file_path):
        """创建Density Radius类型的TabWidget，包含Line和Heatmap两个tab"""
        # 创建新的TabWidget
        density_radius_tab_widget = QTabWidget()
        
        # 设置样式（使用默认样式）
        density_radius_tab_widget.setStyleSheet("""
            QTabBar::tab {
                background: lightgray;
                border: 2px solid #C4C4C3;
                border-bottom-color: #C4C4C3;
                border-top-left-radius: 5px;
                border-top-right-radius: 5px;
                min-width: 16ex;
                padding: 2px;
                font: 18pt "华文细黑";
                color: black;
            }
            QTabBar::tab:selected {
                background: lightblue;
            }
            QTabBar::tab:hover {
                background: pink;
            }
        """)
        density_radius_tab_widget.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        
        # 创建Line tab
        line_tab = self._create_bubble_line_tab()
        density_radius_tab_widget.addTab(line_tab, "📈 Density-Radius Line")
        
        # 创建Heatmap tab
        heatmap_tab = self._create_density_heatmap_tab()
        density_radius_tab_widget.addTab(heatmap_tab, "🗺️ Density-Radius Heatmap")
        
        # 替换当前的TabWidget
        self._replace_widget_in_layout(density_radius_tab_widget)
        self.current_tab_widget = density_radius_tab_widget
        
        # 设置按钮绑定
        self.setup_button_bindings_for_new_tab(density_radius_tab_widget, 'density_radius', file_path)
        
        return density_radius_tab_widget
    
    def _create_density_heatmap_tab(self):
        """创建Density Heatmap tab"""
        # 创建Heatmap tab widget
        heatmap_tab = QWidget()
        heatmap_tab.setObjectName("density_heatmap_tab")
        
        # 创建GridLayout
        grid_layout = QGridLayout(heatmap_tab)
        grid_layout.setObjectName("density_heatmap_gridLayout")
        
        # X-Title
        x_title_label = QLabel(heatmap_tab)
        x_title_label.setObjectName("density_heatmap_label_x")
        x_title_label.setStyleSheet("font: 16pt \"华文细黑\";")
        x_title_label.setText("X-Title")
        grid_layout.addWidget(x_title_label, 0, 0, 1, 1)
        
        x_title_edit = QLineEdit(heatmap_tab)
        x_title_edit.setObjectName("density_heatmap_edit_x")
        x_title_edit.setPlaceholderText("Default if none")
        grid_layout.addWidget(x_title_edit, 0, 1, 1, 1)
        
        # Y-Title
        y_title_label = QLabel(heatmap_tab)
        y_title_label.setObjectName("density_heatmap_label_y")
        y_title_label.setStyleSheet("font: 16pt \"华文细黑\";")
        y_title_label.setText("Y-Title")
        grid_layout.addWidget(y_title_label, 1, 0, 1, 1)
        
        y_title_edit = QLineEdit(heatmap_tab)
        y_title_edit.setObjectName("density_heatmap_edit_y")
        y_title_edit.setPlaceholderText("Default if none")
        grid_layout.addWidget(y_title_edit, 1, 1, 1, 1)
        
        # Color Map
        color_map_label = QLabel(heatmap_tab)
        color_map_label.setObjectName("density_heatmap_label_color_map")
        color_map_label.setStyleSheet("font: 16pt \"华文细黑\";")
        color_map_label.setText("Color Map")
        grid_layout.addWidget(color_map_label, 2, 0, 1, 1)
        
        color_map_combo = QComboBox(heatmap_tab)
        color_map_combo.setObjectName("density_heatmap_combo_color_map")
        color_map_combo.addItems(["viridis", "plasma", "inferno", "magma", "jet", "hot", "cool", "spring", "summer", "autumn", "winter"])
        color_map_combo.setCurrentText("viridis")
        grid_layout.addWidget(color_map_combo, 2, 1, 1, 1)
        
        return heatmap_tab
    

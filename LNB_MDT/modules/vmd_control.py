import sys
import os
import time
import pandas as pd
from socket import socket, AF_INET, SOCK_STREAM
import subprocess
import tracemalloc
import threading
from PySide6.QtWidgets import QApplication, QWidget, QGridLayout, QLabel, QPushButton, QTableWidget, QTableWidgetItem, \
    QSizePolicy, QFrame, QHBoxLayout, QVBoxLayout
from PySide6.QtCore import Qt, QThread, Signal, QMutex, QWaitCondition, QTimer, QElapsedTimer
from PySide6.QtGui import QDragEnterEvent, QDropEvent

# 性能监控类（可选，默认关闭）
class PerformanceMonitor:
    """性能监控工具，记录每个步骤的耗时和系统状态"""
    
    def __init__(self, enabled=False):
        self.enabled = enabled
        self.timer = QElapsedTimer()
        self.steps = []
        self.memory_tracking = False
        
    def start(self, label="Total"):
        """开始计时"""
        if not self.enabled:
            return
        self.timer.start()
        self.steps = []
        self.steps.append(("START", 0))
        
    def checkpoint(self, label):
        """记录检查点"""
        if not self.enabled:
            return
        elapsed = self.timer.elapsed()
        self.steps.append((label, elapsed))
    
    def finish(self):
        """结束计时"""
        if not self.enabled:
            return
        elapsed = self.timer.elapsed()
        
    def get_thread_info(self):
        """获取线程信息"""
        if not self.enabled:
            return ""
        try:
            thread_count = threading.active_count()
            return f"Threads: {thread_count}"
        except:
            return ""

# VMD 命令类
class VMDCommands:
    @staticmethod
    def gotoFrame(frame):
        return f"animate goto {frame}"

    @staticmethod
    def highlightResid(resids):
        """
        高亮显示选中的resid，使用licorice风格和Beta着色
        注意：这不会影响背景representation的coloring method
        """
        resid_str = " ".join(map(str, resids))
        # 使用rep 1来高亮，不影响rep 0的背景
        return f"mol delrep 1 0; mol selection resid {resid_str}; mol representation Licorice; mol color Beta; mol addrep 0"
    
    @staticmethod
    def setBetaValues(resid_value_dict):
        """
        设置每个resid的beta值
        resid_value_dict: {resid: value, ...}
        """
        commands = []
        # 首先清除所有beta值
        commands.append("set sel [atomselect top \"all\"]; $sel set beta 0; $sel delete")
        
        for resid, value in resid_value_dict.items():
            # 为每个resid设置beta值
            cmd = f"set sel [atomselect top \"resid {resid}\"]; $sel set beta {value}; $sel delete"
            commands.append(cmd)
        
        return "; ".join(commands)
    
    @staticmethod
    def setupInitialDisplay():
        """
        设置VMD初始显示：白色背景、红色Points作为背景
        """
        commands = [
            "color Display Background white",  # 设置背景为白色
            "display projection orthographic",  # 设置为正交投影
            "display depthcue off",  # 关闭深度提示
            "mol delrep 0 0",  # 删除默认representation
            "mol representation Points",  # 使用Points表示
            "mol color ColorID 1",  # 红色
            "mol selection all",  # 选择所有atoms
            "mol material Opaque",  # 使用不透明材质
            "mol addrep 0"  # 添加representation
        ]
        return "; ".join(commands)
    
    @staticmethod
    def setupBetaColoring(min_value, max_value):
        """
        设置VMD使用beta值进行着色
        min_value, max_value: 用于设置颜色范围
        颜色映射：蓝色(小值) → 白色(中值) → 红色(大值)
        """
        commands = [
            "mol delrep 0 0",  # 删除当前representation
            "mol representation Points",  # 使用Points表示
            "mol color Beta",  # 使用beta着色
            "mol selection all",  # 选择所有atoms
            "mol material Opaque",  # 使用不透明材质
            "mol addrep 0",  # 添加representation
            "color scale method RWB",  # 设置颜色刻度为Red-White-Blue
            f"mol scaleminmax 0 0 {min_value} {max_value}"  # 设置beta值的范围：min对应蓝色，max对应红色
        ]
        return "; ".join(commands)
    
    @staticmethod
    def showBetaInfo(column_name, min_val, max_val, resid_values):
        """
        在VMD控制台显示beta值信息
        """
        info_commands = [
            f'puts "=== Beta值信息 - {column_name} ==="',
            f'puts "数值范围: {min_val:.6f} 到 {max_val:.6f}"',
            f'puts "已设置 {len(resid_values)} 个resid的beta值"',
            'puts "颜色映射: 蓝色(最小值) → 白色(中值) → 红色(最大值)"',
            'puts "=========================================="'
        ]
        return "; ".join(info_commands)
    
    @staticmethod
    def setColumnBetaValues(column_data, resid_column):
        """
        设置某一列的beta值，直接使用原始数值
        column_data: 某一列的数值数据
        resid_column: resid列的数据
        返回resid_value_dict
        """
        # 直接使用原始数值，不进行归一化
        resid_values = {}
        for resid, value in zip(resid_column, column_data):
            if pd.notna(value) and isinstance(value, (int, float)):
                resid_values[resid] = value
        
        return resid_values

# 后台命令发送线程类
class CommandSenderThread(QThread):
    """在后台线程中发送VMD命令，避免阻塞UI"""
    
    def __init__(self, socket, command):
        super().__init__()
        self.socket = socket
        self.command = command
        
    def run(self):
        """在线程中执行命令发送"""
        try:
            cmd_bytes = (self.command + "\n").encode('utf-8')
            
            # 设置超时，避免无限等待
            try:
                original_timeout = self.socket.gettimeout()
            except (AttributeError, OSError):
                original_timeout = None
            
            self.socket.settimeout(5.0)  # 5秒超时，给大命令足够时间
            
            try:
                self.socket.send(cmd_bytes)
            finally:
                if original_timeout is not None:
                    self.socket.settimeout(original_timeout)
                else:
                    self.socket.settimeout(None)
                    
        except (BrokenPipeError, ConnectionResetError, socket.error) as e:
            print(f"[ERROR] Failed to send command to VMD: {e}")
        except Exception as e:
            print(f"[ERROR] Unexpected error sending command: {e}")

# VMD TCP 控制类
class VMDTcp:
    def __init__(self, rctl_path, vmd_path):
        self.rctl = rctl_path
        self.vmd_path = vmd_path
        self.HOST = 'localhost'
        self.PORT = 5050
        self.ADDR = (self.HOST, self.PORT)
        self.tcpClientSocket = None
        self.vmd_process = None
        # 后台线程相关设置
        self._use_background_thread = True  # 启用后台线程
        self._large_command_threshold = 10000  # 10KB，超过此大小使用后台线程
        self._command_threads = []  # 保存线程引用，避免被垃圾回收

    def attemptConnection(self):
        max_attempts = 5
        for attempt in range(max_attempts):
            try:
                self.tcpClientSocket = socket(AF_INET, SOCK_STREAM)
                self.tcpClientSocket.connect(self.ADDR)
                return 0
            except ConnectionRefusedError:
                time.sleep(1)
                if attempt == max_attempts - 1:
                    return -1

    def start(self):
        if not os.path.exists(self.rctl):
            raise FileNotFoundError(f"remote_ctl.tcl not found at {self.rctl}")
        if not os.path.exists(self.vmd_path):
            raise FileNotFoundError(f"VMD executable not found at {self.vmd_path}")
        self.vmd_process = subprocess.Popen([self.vmd_path, "-e", self.rctl])
        return self.attemptConnection()

    def send_command(self, cmd, use_background=None):
        """
        向VMD发送命令，对于大命令（>10KB），自动使用后台线程发送，避免阻塞UI
        
        Args:
            cmd: 要发送的命令字符串
            use_background: 是否强制使用后台线程（None=自动判断）
        
        Returns:
            "sent"、"sent_async" 或错误信息
        """
        if not self.tcpClientSocket:
            return "error: no socket"
        
        cmd_bytes = (cmd + "\n").encode('utf-8')
        cmd_size = len(cmd_bytes)
        
        # 决定是否使用后台线程
        if use_background is None:
            use_background = self._use_background_thread and cmd_size > self._large_command_threshold
        
        if use_background:
            # 创建后台线程发送命令
            thread = CommandSenderThread(self.tcpClientSocket, cmd)
            self._command_threads.append(thread)  # 保存引用
            thread.start()
            return "sent_async"
        
        # 小命令直接在主线程发送
        try:
            try:
                original_timeout = self.tcpClientSocket.gettimeout()
            except (AttributeError, OSError):
                original_timeout = None
            
            self.tcpClientSocket.settimeout(1.0)
            
            try:
                self.tcpClientSocket.send(cmd_bytes)
                return "sent"
            finally:
                if original_timeout is not None:
                    self.tcpClientSocket.settimeout(original_timeout)
                else:
                    self.tcpClientSocket.settimeout(None)
            
        except (BrokenPipeError, ConnectionResetError, socket.error) as e:
            return f"error: {e}"
        except Exception as e:
            return f"error: unexpected - {e}"

    def stop(self):
        if self.tcpClientSocket:
            self.send_command("quit")
            self.tcpClientSocket.close()
        if self.vmd_process:
            self.vmd_process.terminate()

# CSV 读取函数
def read_excel_vmd(file_path):
    try:
        comments = []
        frame_info = None
        gro_file = None
        xtc_file = None
        
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    comment = line.strip()[2:]
                    comments.append(comment)
                    # 提取FRAME_INFO
                    if comment.startswith('FRAME_INFO:'):
                        frame_info = comment.split(':', 1)[1].split(',')
                    # 提取GRO_FILE
                    elif comment.startswith('GRO_FILE:'):
                        gro_file = comment.split(':', 1)[1]
                    # 提取XTC_FILE
                    elif comment.startswith('XTC_FILE:'):
                        xtc_file = comment.split(':', 1)[1]
                else:
                    break

        df = pd.read_csv(file_path, skiprows=len(comments), header=0)
        
        valid_comments = comments[1] if len(comments) > 1 else ""
        return valid_comments, df, frame_info, gro_file, xtc_file
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return None, None, None

class VMDControlPanel(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.vmd_process = None
        self.vmd_tcp = None
        self.setupUi()
        self.initConnections()
        self.initState()

    def setupUi(self):
        page_vmd = self.parent()
        if not page_vmd:
            print("Error: Parent widget (page_vmd) not found.")
            return

        self.startButton = page_vmd.findChild(QPushButton, "vmd_btn_start")
        self.stopButton = page_vmd.findChild(QPushButton, "vmd_btn_stop")
        self.table = page_vmd.findChild(QTableWidget, "vmd_tablewidget")
        self.infoLabel = page_vmd.findChild(QLabel, "vmd_label")

        scroll_area = page_vmd.findChild(QWidget, "qt_scrollarea")
        if scroll_area:
            scroll_area.setAcceptDrops(False)
            viewport = scroll_area.findChild(QWidget, "qt_scrollarea_viewport")
            if viewport:
                viewport.setAcceptDrops(False)

        if not all([self.startButton, self.stopButton, self.table, self.infoLabel]):
            print("Error: Some UI elements not found. Check object names in Qt Designer.")
            print("Available children in page_vmd:", [child.objectName() for child in page_vmd.findChildren(QWidget)])
            return

        self.table.setSelectionMode(QTableWidget.ExtendedSelection)
        self.table.setVisible(True)
        self.table.setRowCount(0)
        self.table.setColumnCount(0)

        self.stopButton.setEnabled(False)
        self.infoLabel.setText("Click 'Start VMD' to launch VMD, then drag and drop a CSV file")

        self.setAcceptDrops(True)  # 启用拖放

    def initConnections(self):
        self.startButton.clicked.connect(self.pushStartVMD)
        self.stopButton.clicked.connect(self.pushStopVMD)
        self.table.selectionModel().selectionChanged.connect(self.onSelectionChanged)

    def initState(self):
        self.vmd_running = False
        self.vmd_process = None
        self.vmd_tcp = None
        self.df = None
        self.valid_comments = None
        self.frame_info = None
        self._processing_selection = False  # 防止重复处理选择变化
        self._last_selection_time = 0  # 用于防抖

    def pushStartVMD(self):
        try:
            vmd_path = "C:/Program Files/VMD/vmd.exe"
            rctl_path = "path/to/remote_ctl.tcl"  # 替换为 remote_ctl.tcl 的实际路径
            self.vmd_tcp = VMDTcp(rctl_path, vmd_path)
            if self.vmd_tcp.start() == 0:
                self.vmd_running = True
                self.startButton.setEnabled(False)
                self.stopButton.setEnabled(True)
                self.infoLabel.setText("VMD is running. Drag and drop a CSV file to load data.")
                print("VMD started successfully.")
            else:
                self.infoLabel.setText("Failed to connect to VMD.")
                print("Failed to connect to VMD.")
        except Exception as e:
            self.infoLabel.setText(f"Failed to start VMD: {str(e)}")
            print(f"Error starting VMD: {e}")

    def pushStopVMD(self):
        if self.vmd_tcp:
            try:
                self.vmd_tcp.stop()
                self.vmd_tcp = None
                self.vmd_running = False
                self.startButton.setEnabled(True)
                self.stopButton.setEnabled(False)
                self.infoLabel.setText("VMD stopped. Click 'Start VMD' to launch again.")
                print("VMD stopped successfully.")
            except Exception as e:
                self.infoLabel.setText(f"Failed to stop VMD: {str(e)}")
                print(f"Error stopping VMD: {e}")
        else:
            self.infoLabel.setText("VMD is not running.")

    def onSelectionChanged(self, selected, deselected):
        """
        处理表格选择变化事件，使用QTimer延迟处理，避免阻塞UI
        """
        import time
        
        # 防抖：如果正在处理或距离上次处理时间太短，则跳过
        current_time = time.time()
        if self._processing_selection or (current_time - self._last_selection_time < 0.1):
            return
        
        # 检查窗口状态并恢复（如果需要）
        try:
            parent_window = self.parent()
            while parent_window and not hasattr(parent_window, 'isVisible'):
                parent_window = parent_window.parent()
            
            if parent_window and (not parent_window.isVisible() or parent_window.isMinimized()):
                parent_window.show()
                parent_window.raise_()
                parent_window.activateWindow()
        except:
            pass
        
        # 立即处理UI更新，让UI有时间响应
        for i in range(3):
            QApplication.processEvents()
        
        # 使用QTimer延迟处理，让UI先更新
        QTimer.singleShot(50, self._processSelectionDelayed)
    
    def _processSelectionDelayed(self):
        """
        延迟处理选择变化
        """
        import time
        
        self._processing_selection = True
        self._last_selection_time = time.time()
        
        try:
            # 检查基本状态
            if not self.vmd_running or not self.vmd_tcp or self.df is None:
                return
            
            if not self.vmd_tcp.tcpClientSocket:
                self.infoLabel.setText("VMD connection lost. Please restart VMD.")
                return
            
            # 快速获取选中的列
            try:
                # 优先检查列标题选择（更快）
                selected_cols = []
                column_count = self.table.columnCount()
                
                for col in range(column_count):
                    if col > 0:  # 跳过第一列
                        try:
                            if self.table.isColumnSelected(col):
                                selected_cols.append(col)
                        except:
                            pass
                
                # 如果没有通过列标题选择，尝试获取选中的单元格列（限制数量避免阻塞）
                if not selected_cols:
                    try:
                        # 只获取前1000个索引，避免阻塞
                        selected_indexes = list(self.table.selectionModel().selectedIndexes())[:1000]
                        selected_cols = list(set(index.column() for index in selected_indexes))
                        
                        # 如果有很多选中项，检查是否整列被选中
                        if len(selected_indexes) >= 100:
                            col_counts = {}
                            for index in selected_indexes:
                                col = index.column()
                                col_counts[col] = col_counts.get(col, 0) + 1
                            # 如果某个列的选中项数接近总行数，认为是整列选择
                            row_count = self.table.rowCount()
                            for col, count in col_counts.items():
                                if col > 0 and count >= row_count * 0.9:  # 90%以上选中
                                    if col not in selected_cols:
                                        selected_cols.append(col)
                    except Exception as e:
                        pass
            except Exception as e:
                self.infoLabel.setText(f"Error getting selection: {str(e)}")
                return
            
            # 如果选择了列（点击列标题或选择整列）
            if selected_cols:
                unique_cols = list(set(selected_cols))
                
                # 只处理第一个选中的列，避免同时处理多列导致UI卡住
                for col_idx, col in enumerate(unique_cols):
                    if col_idx > 0:  # 只处理第一列
                        break
                    
                    if col > 0:  # 跳过第一列（resid列）
                        try:
                            col_name = self.df.columns[col]
                            
                            if col_name.startswith("Frame_"):
                                frame = int(col_name.split("_")[1])
                                
                                # 跳转到对应帧
                                result = self.vmd_tcp.send_command(VMDCommands.gotoFrame(frame))
                                if result.startswith("error"):
                                    self.infoLabel.setText(f"Failed to send command to VMD: {result}")
                                    break
                                
                                # 获取当前列的数据和resid列
                                column_data = self.df.iloc[:, col].values
                                resid_column = self.df['Resid'].values
                                
                                # 计算beta值（直接使用原始数值）
                                resid_values = VMDCommands.setColumnBetaValues(column_data, resid_column)
                                
                                if resid_values:
                                    # 计算原始数值的范围用于显示信息
                                    valid_values = [val for val in column_data if pd.notna(val) and isinstance(val, (int, float))]
                                    min_val = min(valid_values) if valid_values else 0
                                    max_val = max(valid_values) if valid_values else 0
                                    
                                    # 设置beta值（大命令，使用后台线程）
                                    beta_command = VMDCommands.setBetaValues(resid_values)
                                    result = self.vmd_tcp.send_command(beta_command)
                                    if result.startswith("error"):
                                        self.infoLabel.setText(f"Failed to set beta values: {result}")
                                        break
                                    
                                    # 设置beta着色显示（使用实际数值范围）
                                    coloring_command = VMDCommands.setupBetaColoring(min_val, max_val)
                                    result = self.vmd_tcp.send_command(coloring_command)
                                    if result.startswith("error"):
                                        self.infoLabel.setText(f"Failed to set coloring: {result}")
                                        break
                                    
                                    # 显示beta值信息
                                    info_command = VMDCommands.showBetaInfo(col_name, min_val, max_val, resid_values)
                                    self.vmd_tcp.send_command(info_command)  # 非关键命令，忽略错误
                                    
                                    self.infoLabel.setText(f"Frame {frame}: Beta coloring applied (range: {min_val:.3f} to {max_val:.3f})")
                        except Exception as e:
                            self.infoLabel.setText(f"Error processing column: {str(e)}")
                            break
            
        except Exception as e:
            self.infoLabel.setText(f"Selection error: {str(e)}")
        finally:
            self._processing_selection = False

    def dragEnterEvent(self, event: QDragEnterEvent):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
        else:
            event.ignore()

    def dropEvent(self, event: QDropEvent):
        if not self.vmd_running:
            self.infoLabel.setText("Please start VMD before dropping a file.")
            event.ignore()
            return

        urls = event.mimeData().urls()
        if not urls:
            event.ignore()
            return

        file_path = urls[0].toLocalFile()
        if file_path.lower().endswith('.csv'):
            self.loadCSV(file_path)
            event.acceptProposedAction()
        else:
            self.infoLabel.setText("Please drop a CSV file.")
            event.ignore()

    def loadCSV(self, file_path):
        try:
            if not os.path.exists(file_path):
                self.infoLabel.setText("CSV file not found!")
                return

            self.valid_comments, self.df, self.frame_info = read_excel_vmd(file_path)
            if self.df is None:
                self.infoLabel.setText("Error: Failed to read CSV file.")
                return

            # 调整列名以匹配文件中的大小写
            self.df.rename(columns={'Resid': 'Resid', 'Resname': 'Resname', 'Coordinates': 'Coordinates'}, inplace=True)
            # 忽略 resname 和 coordinates 列（保留原始列名）
            self.df = self.df.drop(columns=['Resname', 'Coordinates'], errors='ignore')
            
            # 如果有frame_info，将Time列标题替换为Frame列标题
            if self.frame_info:
                frame_cols = [col for col in self.df.columns if col != 'Resid']
                
                # 创建新的列名映射：Time列 -> Frame列
                column_mapping = {}
                for i, time_col in enumerate(frame_cols):
                    if i < len(self.frame_info) and self.frame_info[i]:
                        frame_value = self.frame_info[i]
                        column_mapping[time_col] = f"Frame_{frame_value}"
                
                
                # 重命名列
                self.df.rename(columns=column_mapping, inplace=True)
                frame_cols = [col for col in self.df.columns if col != 'Resid']
            else:
                frame_cols = [col for col in self.df.columns if col != 'Resid']
            
            self.displayData(frame_cols)
            
            frame_count = len(self.frame_info) if self.frame_info else 0
            self.infoLabel.setText(f"CSV loaded successfully. {frame_count} frames detected. Valid comment: {self.valid_comments}")
        except Exception as e:
            self.infoLabel.setText(f"Error loading CSV: {str(e)}")
            print(f"Error loading CSV: {e}")

    def displayData(self, frame_cols):
        print("displayData called")
        self.table.clear()
        self.table.setRowCount(len(self.df))
        self.table.setColumnCount(len(frame_cols) + 1)
        self.table.setHorizontalHeaderLabels(['Resid'] + frame_cols)

        for i, row in self.df.iterrows():
            self.table.setItem(i, 0, QTableWidgetItem(str(row['Resid'])))
            for j, frame in enumerate(frame_cols):
                self.table.setItem(i, j + 1, QTableWidgetItem(str(row[frame])))

        self.table.resizeColumnsToContents()
        self.table.resizeRowsToContents()
        self.table.setVisible(True)
        self.table.update()
import sys
import os
import time
import pandas as pd
from socket import socket, AF_INET, SOCK_STREAM
import subprocess
from PySide6.QtWidgets import QApplication, QWidget, QGridLayout, QLabel, QPushButton, QTableWidget, QTableWidgetItem, \
    QSizePolicy, QFrame, QHBoxLayout, QVBoxLayout
from PySide6.QtCore import Qt
from PySide6.QtGui import QDragEnterEvent, QDropEvent

# VMD 命令类
class VMDCommands:
    @staticmethod
    def gotoFrame(frame):
        return f"animate goto {frame}"

    @staticmethod
    def highlightResid(resids):
        """
        高亮显示选中的resid，使用VDW和Beta着色
        注意：这不会影响背景representation的coloring method
        """
        resid_str = " ".join(map(str, resids))
        # 使用rep 1来高亮，不影响rep 0的背景
        return f"mol delrep 1 0; mol selection resid {resid_str}; mol representation VDW; mol color Beta; mol addrep 0"
    
    @staticmethod
    def setBetaValues(resid_value_dict):
        """
        设置每个resid的beta值
        resid_value_dict: {resid: value, ...}
        """
        commands = []
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

    def send_command(self, cmd):
        if self.tcpClientSocket:
            self.tcpClientSocket.send((cmd + "\n").encode())
            # VMD命令通常是异步的，不需要等待响应
            return "sent"

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
        selected_rows = [index.row() for index in self.table.selectionModel().selectedRows()]
        selected_cols = [index.column() for index in self.table.selectionModel().selectedIndexes()]
        
        if self.vmd_running and self.vmd_tcp and self.df is not None:
            for row in selected_rows:
                for col in selected_cols:
                    if col > 0:  # 跳过第一列（resid列）
                        # 直接从列标题获取frame号
                        col_name = self.df.columns[col]
                        if col_name.startswith("Frame_"):
                            frame = int(col_name.split("_")[1])
                            resid = self.df.iloc[row]["Resid"]
                            self.vmd_tcp.send_command(VMDCommands.gotoFrame(frame))
                            self.vmd_tcp.send_command(VMDCommands.highlightResid([int(resid)]))
                        else:
                            pass
        else:
            pass

    def dragEnterEvent(self, event: QDragEnterEvent):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
        else:
            event.ignore()

    def dropEvent(self, event: QDropEvent):
        print("dropEvent triggered in VMDControlPanel")
        if not self.vmd_running:
            self.infoLabel.setText("Please start VMD before dropping a file.")
            event.ignore()
            return

        urls = event.mimeData().urls()
        if not urls:
            event.ignore()
            return

        file_path = urls[0].toLocalFile()
        print(f"File path: {file_path}")
        if file_path.lower().endswith('.csv'):
            self.loadCSV(file_path)
            event.acceptProposedAction()
        else:
            self.infoLabel.setText("Please drop a CSV file.")
            event.ignore()

    def loadCSV(self, file_path):
        print("loadCSV called")
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
import sys
from PySide6.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QPushButton, QStackedWidget
from PySide6.QtCore import Qt

# ----------------------------------------------------------------------
# 这是一个包含了“解决方案”的样式表
# 关键在于第一条规则：#pageContainer { background-color: #f0f0f0; }
# 它为容纳按钮的父容器提供了一个不透明的“画布”
# ----------------------------------------------------------------------
STYLESHEET = """
#pageContainer {
    /* 关键修复：给页面容器一个坚实的不透明背景 */
    background-color: #f8f8f2; 
    border-radius: 10px; /* 圆角只是为了让效果更清晰 */
}

QPushButton {
    /* 这是我们希望按钮最终呈现的样式 */
    background-color: #6272a4; /* 您想要的蓝色 */
    color: white;
    font-size: 16px;
    padding: 10px;
    border-radius: 5px;
    border: 2px solid #506090;
}
QPushButton:hover {
    background-color: #7081b6;
}
QPushButton:pressed {
    background-color: #506090;
}
"""

class MinimalApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Minimal Working Example")
        self.resize(500, 300)

        # 1. 模拟您的主窗口环境：无边框 + 透明背景
        self.setWindowFlags(Qt.FramelessWindowHint)
        self.setAttribute(Qt.WA_TranslucentBackground)

        # 2. 模拟您的核心布局：使用QStackedWidget作为中心
        container = QStackedWidget()
        # 我们给这个容器一个ObjectName，以便在QSS中定位它
        container.setObjectName("pageContainer")

        # 3. 创建一个页面和布局
        page = QWidget()
        layout = QVBoxLayout(page)

        # 4. 在页面上创建一个按钮（这是我们一直无法染色的按钮）
        button = QPushButton("如果这个按钮是蓝色的，我们就成功了")
        layout.addWidget(button, alignment=Qt.AlignCenter)

        # 5. 将页面添加到容器中，并将容器设为窗口的中心控件
        container.addWidget(page)
        self.setCentralWidget(container)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    
    # 将样式表应用到整个程序
    app.setStyleSheet(STYLESHEET)
    
    window = MinimalApp()
    window.show()
    sys.exit(app.exec())
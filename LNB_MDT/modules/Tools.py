from PySide6.QtCore import Qt
from PySide6.QtGui import QCursor
from PySide6.QtWidgets import *


class UISettings:
    FONT_SIZE = '15pt'
    BTN_BORDER_RADIUS = '15px'
    BTN_WIDTH = '30px'
    BTN_HEIGHT = '30px'
    LABEL_HEIGHT = 20
    SPIN_VALUE = 0
    SPIN_MIN = 0
    SPIN_MAX = 10000000
    SPIN_STEP = 1


class UIItemsMake:
    @classmethod
    def make_widget(cls):
        widget = QWidget()
        return widget

    @classmethod
    def make_btn(cls,
                 btnName: str
                 , callback=None
                 , layout=None
                 , **kwargs):

        settings = {
            "background_color": "#9faeda",  # 统一背景色
            "font_size": UISettings.FONT_SIZE,
            "border_radius": "5px",  # 统一圆角
            "width": UISettings.BTN_WIDTH,
            "height": UISettings.BTN_HEIGHT,
            "font_color": "black",  # 统一文字颜色
            "border_color": "#9faeda"  # 统一边框颜色
        }

        settings.update(**kwargs)
        btn = QPushButton(btnName)
        # 设置按钮样式，确保动态创建的按钮也有统一的样式
        style_sheet = (
            f"background-color: {settings['background_color']} !important;"
            f"color: {settings['font_color']} !important;"
            f"border: 2px solid {settings['border_color']} !important;"
            f"border-radius: {settings['border_radius']} !important;"
            f"font: 12pt \"华文细黑\";"
            f"padding: 5px;"
        )
        btn.setStyleSheet(style_sheet)
        btn.setCursor(QCursor(Qt.PointingHandCursor))
        # 确保按钮填充背景（macOS兼容）
        btn.setAutoFillBackground(True)
        if callback:
            btn.clicked.connect(callback)
        if layout:
            layout.addWidget(btn)
        return btn

    @classmethod
    def make_label(cls, text: str, height=UISettings.LABEL_HEIGHT, **kwargs):
        label = QLabel(text)
        label.setMaximumHeight(height)
        return label

    @classmethod
    def make_group_box(cls, title: str):
        group_box = QGroupBox(title)
        return group_box

    @classmethod
    def make_radio_check(cls, btn_type, text):
        radio_check = btn_type(text)
        return radio_check

    @classmethod
    def make_spin_box(cls, value=UISettings.SPIN_VALUE, min=UISettings.SPIN_MIN, 
                      max=UISettings.SPIN_MAX, step=UISettings.SPIN_STEP):
        spin_box = QSpinBox()
        spin_box.setValue(value)
        spin_box.setMinimum(min)
        spin_box.setMaximum(max)
        spin_box.setSingleStep(step)
        return spin_box


def create_warn_dialog(text, title="Warning"):
    app = QApplication.instance()
    if not app:
        app = QApplication([])

    dialog = QDialog()
    dialog.setWindowTitle(title)

    layout = QVBoxLayout()

    label = QLabel(text)
    layout.addWidget(label)

    ok_button = QPushButton("OK")
    ok_button.clicked.connect(dialog.accept)  # 连接按钮点击事件
    layout.addWidget(ok_button)

    dialog.setLayout(layout)

    dialog.exec()

from PySide6.QtCore import Qt
from PySide6.QtGui import QCursor
from PySide6.QtWidgets import *


class UISettings:
    FONT_COLOR = 'white'
    FONT_SIZE = '15pt'
    # widget
    WIDGET_COLOR = 'rgb(33, 37, 43)'  # 浅黑色
    # btn
    BTN_COLOR = '#6272a4'  # 浅蓝色
    BTN_BORDER_RADIUS = '15px'
    BTN_WIDTH = '30px'
    BTN_HEIGHT = '30px'
    #label
    LABEL_HEIGHT = 20
    # spin_box
    SPIN_VALUE = 0
    SPIN_MIN = 0
    SPIN_MAX = 10000000
    SPIN_STEP = 1


class UIItemsMake:
    # font

    @classmethod
    def make_widget(cls
                    , color=UISettings.WIDGET_COLOR
                    , font_size=UISettings.FONT_SIZE):
        widget = QWidget()
        return widget

    @classmethod
    def make_btn(cls,
                 btnName: str
                 , callback=None
                 , layout=None
                 , **kwargs):

        settings = {
            "background_color": UISettings.BTN_COLOR,
            "font_size": UISettings.FONT_SIZE,
            "border_radius": UISettings.BTN_BORDER_RADIUS,
            "width": UISettings.BTN_WIDTH,
            "height": UISettings.BTN_HEIGHT,
            "font_color": UISettings.FONT_COLOR
        }

        settings.update(**kwargs)
        btn = QPushButton(btnName)
        btn.setStyleSheet(
            f"background-color:{settings['background_color']};"
            f"font-size:{settings['font_size']};"
            f"border-radius:{settings['border_radius']};"
            f"width:{settings['width']};"
            f"height:{settings['height']};"
        )
        btn.setCursor(QCursor(Qt.PointingHandCursor))
        if callback:
            btn.clicked.connect(callback)
        if layout:
            layout.addWidget(btn)
        return btn

    @classmethod
    def make_label(cls
                   , text: str
                   , height=UISettings.LABEL_HEIGHT
                   , **kwargs):
        settings = {
            "font_size": UISettings.FONT_SIZE
        }
        settings.update(**kwargs)
        label = QLabel(text)
        label.setStyleSheet(
            f"font-size:{settings['font_size']};")
        label.setMaximumHeight(height)
        label.setMaximumHeight(height)
        return label

    @classmethod
    def make_group_box(cls
                       , title: str
                       , title_color='white'):
        group_box = QGroupBox(title)
        group_box.setStyleSheet("""QGroupBox {font-size:14pt;}
                                   QGroupBox::title {color: %s;}
                                   QGroupBox QCheckBox {color: white;}
                                   QGroupBox QRadioButton {color: white;}
                                """ % title_color)
        return group_box

    @classmethod
    def make_radio_check(cls
                         , btn_type
                         , text
                         , font_size=UISettings.FONT_SIZE
                         ):
        radio_check = btn_type(text)
        radio_check.setStyleSheet('font-size:%s;'
                                  % (font_size))
        return radio_check

    @classmethod
    def make_spin_box(cls
                      , value=UISettings.SPIN_VALUE
                      , min=UISettings.SPIN_MIN
                      , max=UISettings.SPIN_MAX
                      , step=UISettings.SPIN_STEP):
        spin_box = QSpinBox()
        spin_box.setStyleSheet("""
                            QSpinBox {
                                font: 16pt "华文细黑";
                                color: white;
                                border: 1px solid white;
                                background-color: rgb(33, 37, 43);
                            }
                            QSpinBox::up-button {
                                height: 10px; /* 设置一个合适的高度 */
                                width: 20px; /* 设置一个合适的宽度 */
                                subcontrol-position: top right; /* 保持按钮在右侧 */
                            }
                            QSpinBox::down-button {
                                height: 10px; /* 设置一个合适的高度 */
                                width: 20px; /* 设置一个合适的宽度 */
                                subcontrol-position: bottom right; /* 保持按钮在右侧 */
                            }
                        """)
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

# MAIN FILE
# ///////////////////////////////////////////////////////////////
from main import *

# GLOBALS
# ///////////////////////////////////////////////////////////////
GLOBAL_STATE = False
GLOBAL_TITLE_BAR = True

class UIFunctions(MainWindow):
    # MAXIMIZE/RESTORE
    # ///////////////////////////////////////////////////////////////
    def maximize_restore(self):
        global GLOBAL_STATE
        status = GLOBAL_STATE
        if status == False:
            self.showMaximized()
            GLOBAL_STATE = True
            self.ui.appMargins.setContentsMargins(0, 0, 0, 0)
            self.ui.maximizeRestoreAppBtn.setToolTip("Restore")
            self.ui.maximizeRestoreAppBtn.setIcon(QIcon(u":/icons/images/icons/icon_restore.png"))
            self.ui.frame_size_grip.hide()
            self.left_grip.hide()
            self.right_grip.hide()
            self.top_grip.hide()
            self.bottom_grip.hide()
        else:
            GLOBAL_STATE = False
            self.showNormal()
            self.resize(self.width()+1, self.height()+1)
            self.ui.appMargins.setContentsMargins(10, 10, 10, 10)
            self.ui.maximizeRestoreAppBtn.setToolTip("Maximize")
            self.ui.maximizeRestoreAppBtn.setIcon(QIcon(u":/icons/images/icons/icon_maximize.png"))
            self.ui.frame_size_grip.show()
            self.left_grip.show()
            self.right_grip.show()
            self.top_grip.show()
            self.bottom_grip.show()

    # RETURN STATUS
    # ///////////////////////////////////////////////////////////////
    def returStatus(self):
        return GLOBAL_STATE

    # SET STATUS
    # ///////////////////////////////////////////////////////////////
    def setStatus(self, status):
        global GLOBAL_STATE
        GLOBAL_STATE = status

    # TOGGLE MENU
    # ///////////////////////////////////////////////////////////////
    def toggleMenu(self, enable):
        if enable:
            # GET WIDTH
            width = self.ui.leftMenuBg.width()
            maxExtend = Settings.MENU_WIDTH
            standard = 60

            # SET MAX WIDTH
            if width == 60:
                widthExtended = maxExtend
            else:
                widthExtended = standard

            # ANIMATION
            self.animation = QPropertyAnimation(self.ui.leftMenuBg, b"minimumWidth")
            self.animation.setDuration(Settings.TIME_ANIMATION)
            self.animation.setStartValue(width)
            self.animation.setEndValue(widthExtended)
            self.animation.setEasingCurve(QEasingCurve.InOutQuart)
            self.animation.start()

    # TOGGLE LEFT BOX
    # ///////////////////////////////////////////////////////////////
    def toggleGeneRightBox(self, enable):
        if enable:
            # GET WIDTH
            width = self.ui.Gene_extraBox.width()
            maxExtend = Settings.GENER_RIGHT_BOX_WIDTH
            standard = 0

            # SET MAX WIDTH
            if width == 0:
                widthExtended = maxExtend
            else:
                widthExtended = standard

            # ANIMATION
            self.Generation_rihtBox = QPropertyAnimation(self.ui.Gene_extraBox, b"minimumWidth")
            self.Generation_rihtBox.setDuration(Settings.TIME_ANIMATION)
            self.Generation_rihtBox.setStartValue(width)
            self.Generation_rihtBox.setEndValue(widthExtended)
            self.Generation_rihtBox.setEasingCurve(QEasingCurve.InOutQuart)
            self.Generation_rihtBox.start()

    def toggleFigureColorBox(self, enable):
        if enable:
            # GET WIDTH
            width = self.ui.figure_color_extra_box.width()
            maxExtend = Settings.FIGURE_RIGHT_BOX_WIDTH
            standard = 0


            # SET MAX WIDTH
            if width == 0:
                widthExtended = maxExtend
            else:
                widthExtended = standard

            # ANIMATION
            self.Figure_colorBox = QPropertyAnimation(self.ui.figure_color_extra_box, b"minimumWidth")
            self.Figure_colorBox.setDuration(Settings.TIME_ANIMATION)
            self.Figure_colorBox.setStartValue(width)
            self.Figure_colorBox.setEndValue(widthExtended)
            self.Figure_colorBox.setEasingCurve(QEasingCurve.InOutQuart)
            self.Figure_colorBox.start()
        else:
            if self.ui.figure_color_extra_box.width() != 0:
                self.Figure_colorBox = QPropertyAnimation(self.ui.figure_color_extra_box, b"minimumWidth")
                self.Figure_colorBox.setDuration(Settings.TIME_ANIMATION)
                self.Figure_colorBox.setStartValue(self.ui.figure_color_extra_box.width())
                self.Figure_colorBox.setEndValue(0)
                self.Figure_colorBox.setEasingCurve(QEasingCurve.InOutQuart)
                self.Figure_colorBox.start()

    def toggleFigureShapeBox(self, enable):
        if enable:
            # GET WIDTH
            width = self.ui.figure_shape_extra_box.width()
            maxExtend = Settings.FIGURE_RIGHT_BOX_WIDTH
            standard = 0

            # SET MAX WIDTH
            if width == 0:
                widthExtended = maxExtend
            else:
                widthExtended = standard

            # ANIMATION
            self.Figure_shapeBox = QPropertyAnimation(self.ui.figure_shape_extra_box, b"minimumWidth")
            self.Figure_shapeBox.setDuration(Settings.TIME_ANIMATION)
            self.Figure_shapeBox.setStartValue(width)
            self.Figure_shapeBox.setEndValue(widthExtended)
            self.Figure_shapeBox.setEasingCurve(QEasingCurve.InOutQuart)
            self.Figure_shapeBox.start()
        else:
            if self.ui.figure_shape_extra_box.width() != 0:
                self.Figure_shapeBox = QPropertyAnimation(self.ui.figure_shape_extra_box, b"minimumWidth")
                self.Figure_shapeBox.setDuration(Settings.TIME_ANIMATION)
                self.Figure_shapeBox.setStartValue(self.ui.figure_shape_extra_box.width())
                self.Figure_shapeBox.setEndValue(0)
                self.Figure_shapeBox.setEasingCurve(QEasingCurve.InOutQuart)
                self.Figure_shapeBox.start()


    # TOGGLE RIGHT BOX
    # ///////////////////////////////////////////////////////////////
    def toggleRightBox(self, enable):
        if enable:
            # GET WIDTH
            width = self.ui.extraRightBox.width()
            maxExtend = Settings.RIGHT_BOX_WIDTH
            standard = 0

            # SET MAX WIDTH
            if width == 0:
                widthExtended = maxExtend
            else:
                widthExtended = standard

            # ANIMATION
            self.Analysis_rightBox = QPropertyAnimation(self.ui.extraRightBox, b"minimumWidth")
            self.Analysis_rightBox.setDuration(Settings.TIME_ANIMATION)
            self.Analysis_rightBox.setStartValue(width)
            self.Analysis_rightBox.setEndValue(widthExtended)
            self.Analysis_rightBox.setEasingCurve(QEasingCurve.InOutQuart)
            self.Analysis_rightBox.start()

    def selectMenu(getStyle):
        select = getStyle + Settings.MENU_SELECTED_STYLESHEET
        return select

    # DESELECT
    def deselectMenu(getStyle):
        deselect = getStyle.replace(Settings.MENU_SELECTED_STYLESHEET, "")
        return deselect

    # START SELECTION
    def selectStandardMenu(self, widget):
        for w in self.ui.topMenu.findChildren(QPushButton):
            if w.objectName() == widget:
                w.setStyleSheet(UIFunctions.selectMenu(w.styleSheet()))

    # RESET SELECTION
    def resetStyle(self, widget):
        for w in self.ui.topMenu.findChildren(QPushButton):
            if w.objectName() != widget:
                w.setStyleSheet(UIFunctions.deselectMenu(w.styleSheet()))

    # IMPORT THEMES FILES QSS/CSS
    # ///////////////////////////////////////////////////////////////
    def theme(self, file, useCustomTheme):
        if useCustomTheme:
            str = open(file, 'r').read()
            self.ui.styleSheet.setStyleSheet(str)

    # START - GUI DEFINITIONS
    # ///////////////////////////////////////////////////////////////
    def uiDefinitions(self):
        def dobleClickMaximizeRestore(event):
            # IF DOUBLE CLICK CHANGE STATUS
            if event.type() == QEvent.MouseButtonDblClick:
                QTimer.singleShot(250, lambda: UIFunctions.maximize_restore(self))
        self.ui.titleRightInfo.mouseDoubleClickEvent = dobleClickMaximizeRestore

        if Settings.ENABLE_CUSTOM_TITLE_BAR:
            #STANDARD TITLE BAR
            self.setWindowFlags(Qt.FramelessWindowHint)
            self.setAttribute(Qt.WA_TranslucentBackground)

            # MOVE WINDOW / MAXIMIZE / RESTORE
            def moveWindow(event):
                # IF MAXIMIZED CHANGE TO NORMAL
                if UIFunctions.returStatus(self):
                    UIFunctions.maximize_restore(self)
                # MOVE WINDOW
                if event.buttons() == Qt.LeftButton:
                    self.move(self.pos() + event.globalPos() - self.dragPos)
                    self.dragPos = event.globalPos()
                    event.accept()
            self.ui.titleRightInfo.mouseMoveEvent = moveWindow

            # CUSTOM GRIPS
            self.left_grip = CustomGrip(self, Qt.LeftEdge, True)
            self.right_grip = CustomGrip(self, Qt.RightEdge, True)
            self.top_grip = CustomGrip(self, Qt.TopEdge, True)
            self.bottom_grip = CustomGrip(self, Qt.BottomEdge, True)

        else:
            self.ui.appMargins.setContentsMargins(0, 0, 0, 0)
            self.ui.minimizeAppBtn.hide()
            self.ui.maximizeRestoreAppBtn.hide()
            self.ui.closeAppBtn.hide()
            self.ui.frame_size_grip.hide()

        # DROP SHADOW
        self.shadow = QGraphicsDropShadowEffect(self)
        self.shadow.setBlurRadius(17)
        self.shadow.setXOffset(0)
        self.shadow.setYOffset(0)
        self.shadow.setColor(QColor(0, 0, 0, 150))
        self.ui.bgApp.setGraphicsEffect(self.shadow)

        # RESIZE WINDOW
        # self.sizegrip = QSizeGrip(self.ui.frame_size_grip)
        # self.sizegrip.setStyleSheet("width: 20px; height: 20px; margin 0px; padding: 0px;")

        # MINIMIZE
        self.ui.minimizeAppBtn.clicked.connect(lambda: self.showMinimized())

        # MAXIMIZE/RESTORE
        self.ui.maximizeRestoreAppBtn.clicked.connect(lambda: UIFunctions.maximize_restore(self))

        # CLOSE APPLICATION
        self.ui.closeAppBtn.clicked.connect(lambda: self.close())

    def resize_grips(self):
        if Settings.ENABLE_CUSTOM_TITLE_BAR:
            self.left_grip.setGeometry(0, 10, 10, self.height())
            self.right_grip.setGeometry(self.width() - 10, 10, 10, self.height())
            self.top_grip.setGeometry(0, 0, self.width(), 10)
            self.bottom_grip.setGeometry(0, self.height() - 10, self.width(), 10)

    # ///////////////////////////////////////////////////////////////
    # END - GUI DEFINITIONS

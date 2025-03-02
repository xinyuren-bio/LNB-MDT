

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
    # self.setGeometry(200,200,200,400)
    # self.MainLayout = QVBoxLayout()
    # # 创建按钮
    # self.scrollWidget = QWidget()
    # self.GeneSrcollArea = QScrollArea()
    # self.GeneSrcollArea.setWidget(self.scrollWidget)
    # self.GeneSrcollArea.setWidgetResizable(True)
    #
    # self.mainLayout = QVBoxLayout()
    #
    # self.button_pc = QPushButton("PC ▼")
    # self.button_pe = QPushButton("PE ▼")
    # self.button_pg = QPushButton("PG ▼")
    # self.button_ps = QPushButton("PS ▼")
    # self.button_other = QPushButton("Other ▼")
    # # 创建用于显示内容的容器
    # self.container_pc = QWidget()
    # self.container_pe = QWidget()
    # self.container_pg = QWidget()
    # self.container_ps = QWidget()
    # self.container_other = QWidget()
    #
    # self.grid_pc = QGridLayout(self.container_pc)
    # self.grid_pe = QGridLayout(self.container_pe)
    # self.grid_pg = QGridLayout(self.container_pg)
    # self.grid_ps = QGridLayout(self.container_ps)
    # self.grid_other = QGridLayout(self.container_other)
    #
    # for i in pc:
    #     label = QLabel(i)
    #     spin = QSpinBox()
    #     setattr(self, i, spin)
    #     self.grid_pc.addWidget(label, pc.index(i), 0)
    #     self.grid_pc.addWidget(spin, pc.index(i), 1)
    # for i in pe:
    #     label = QLabel(i)
    #     spin = QSpinBox()
    #     setattr(self, i, spin)
    #     self.grid_pe.addWidget(label, pe.index(i), 0)
    #     self.grid_pe.addWidget(spin, pe.index(i), 1)
    # for i in pg:
    #     label = QLabel(i)
    #     spin = QSpinBox()
    #     setattr(self, i, spin)
    #     self.grid_pg.addWidget(label, pg.index(i), 0)
    #     self.grid_pg.addWidget(spin, pg.index(i), 1)
    # for i in ps:
    #     label = QLabel(i)
    #     spin = QSpinBox()
    #     setattr(self, i, spin)
    #     self.grid_ps.addWidget(label, ps.index(i), 0)
    #     self.grid_ps.addWidget(spin, ps.index(i), 1)
    # for i in other:
    #     label = QLabel(i)
    #     spin = QSpinBox()
    #     setattr(self, i, spin)
    #     self.grid_other.addWidget(label, other.index(i), 0)
    #     self.grid_other.addWidget(spin, other.index(i), 1)
    #
    #
    # # 主布局
    # self.mainLayout.addWidget(self.GeneSrcollArea)
    #
    # self.mainLayout.addWidget(self.button_pc)
    # self.mainLayout.addWidget(self.container_pc)
    #
    # self.mainLayout.addWidget(self.button_pe)
    # self.mainLayout.addWidget(self.container_pe)
    #
    # self.mainLayout.addWidget(self.button_pg)
    # self.mainLayout.addWidget(self.container_pg)
    #
    # self.mainLayout.addWidget(self.button_ps)
    # self.mainLayout.addWidget(self.container_ps)
    #
    # self.mainLayout.addWidget(self.button_other)
    # self.mainLayout.addWidget(self.container_other)
    #
    # # 中心部件
    # self.container_pc.hide()
    # self.container_pe.hide()
    # self.container_pg.hide()
    # self.container_ps.hide()
    # self.container_other.hide()
    #
    # # 连接信号与槽
    # self.button_pc.clicked.connect(lambda :self.toggle_container('container_pc'))
    # self.button_pe.clicked.connect(lambda :self.toggle_container('container_pe'))
    # self.button_pg.clicked.connect(lambda :self.toggle_container('container_pg'))
    # self.button_ps.clicked.connect(lambda :self.toggle_container('container_ps'))
    # self.button_other.clicked.connect(lambda :self.toggle_container('container_other'))
    #
    # self.MainLayout.addWidget(self.scrollWidget)
    # self.setLayout(self.MainLayout)
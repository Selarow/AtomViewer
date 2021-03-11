from dat import DAT
from xml_parser import XML
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import pyqtgraph as pg
import time
import sys




class Window(QMainWindow):
    def __init__(self):
        super().__init__()
        stylesheet = open("data/style.css", "r").read()
        self.setWindowTitle("Atom Viewer")
        self.setGeometry(50, 50, 1600, 980)
        self.setStyleSheet(stylesheet)
        self.UI()
        self.show()

        #VARIABLES
        self.xml = XML()
        self.dat = DAT()
        self.xml_transition = 1
        self.dat_transition = 1
        self.dat_visible = False
        self.colors = ["r", "y", "g", "c", "b", "m"]



    def UI(self):
        #WIDGETS/LAYOUTS
        main_widget = QWidget()
        main_layout = QGridLayout(main_widget)
        self.setCentralWidget(main_widget)

        control_top_widget = QWidget()
        control_top_layout = QHBoxLayout(control_top_widget)
        control_top_layout.setContentsMargins(5, 1, 5, 1)
        main_layout.addWidget(control_top_widget, 2, 1)

        control_bot_widget = QWidget()
        control_bot_layout = QHBoxLayout(control_bot_widget)
        control_bot_layout.setContentsMargins(5, 1, 5, 1)
        main_layout.addWidget(control_bot_widget, 3, 1)

        dat_widget = QWidget(self)
        dat_layout = QVBoxLayout(dat_widget)
        dat_widget.setGeometry(-700, 0, 700, 980)
        dat_widget.setObjectName("dat")

        dat_control_widget = QWidget()
        dat_control_layout = QHBoxLayout(dat_control_widget)
        dat_control_layout.setContentsMargins(5, 0, 5, 0)
        dat_layout.addWidget(dat_control_widget)


        #BUTTONS
        open_vasprun_button = QPushButton("Open Vasprun", self)
        open_vasprun_button.clicked.connect(self.open_vasprun)
        control_top_layout.addWidget(open_vasprun_button)

        dat_button = QPushButton("DAT", self)
        dat_button.clicked.connect(self.toggle_dat)
        control_top_layout.addWidget(dat_button)

        clear_vasp_button = QPushButton("Clear", self)
        clear_vasp_button.clicked.connect(self.clear_vasp)
        control_top_layout.addWidget(clear_vasp_button)

        save_button = QPushButton("Save", self)
        save_button.clicked.connect(self.save)
        control_bot_layout.addWidget(save_button)

        movie_button = QPushButton("Movie", self)
        movie_button.clicked.connect(self.movie)
        control_bot_layout.addWidget(movie_button)

        open_dat_button = QPushButton("Open .dat", self)
        open_dat_button.clicked.connect(self.open_dat)
        dat_control_layout.addWidget(open_dat_button)

        clear_dat_button = QPushButton("Clear", self)
        clear_dat_button.clicked.connect(self.clear_dat)
        dat_control_layout.addWidget(clear_dat_button)


        #ANIMATIONS
        self.close_anim = QPropertyAnimation(dat_widget, b"pos")
        self.close_anim.setEasingCurve(QEasingCurve.InOutCubic)
        self.close_anim.setEndValue(QPoint(-700,0))
        self.close_anim.setDuration(500)

        self.open_anim = QPropertyAnimation(dat_widget, b"pos")
        self.open_anim.setEasingCurve(QEasingCurve.InOutCubic)
        self.open_anim.setEndValue(QPoint(0,0))
        self.open_anim.setDuration(500)


        #GRAPHS VASPRUN
        font = QFontDatabase.addApplicationFont("data/lmromanb.ttf")
        labelStyle = {'color': '#969696', 'font-size': '16px', 'font-family': 'LM Roman 10'}

        self.graph_energy = pg.PlotWidget()
        self.graph_entrp = pg.PlotWidget()
        self.graph_total = pg.PlotWidget()

        self.graph_energy.showGrid(x = True, y = True, alpha = 0.8)
        self.graph_entrp.showGrid(x = True, y = True, alpha = 0.8)
        self.graph_total.showGrid(x = True, y = True, alpha = 0.8)

        self.graph_energy.setLabel("bottom", "t", units = "x", **labelStyle)
        self.graph_entrp.setLabel("bottom", "t", units = "x", **labelStyle)
        self.graph_total.setLabel("bottom", "t", units = "x", **labelStyle)

        self.graph_energy.setLabel("left", "energy", units = "y", **labelStyle)
        self.graph_entrp.setLabel("left", "entrp", units = "y", **labelStyle)
        self.graph_total.setLabel("left", "total", units = "y", **labelStyle)

        main_layout.addWidget(self.graph_energy, 0, 1)
        main_layout.addWidget(self.graph_entrp, 0, 0)
        main_layout.addWidget(self.graph_total, 1, 1)


        #GRAPHS DAT
        self.graph_T = pg.PlotWidget()
        self.graph_EK = pg.PlotWidget()
        self.graph_SP = pg.PlotWidget()
        self.graph_SK = pg.PlotWidget()

        self.graph_T.showGrid(x = True, y = True, alpha = 0.8)
        self.graph_EK.showGrid(x = True, y = True, alpha = 0.8)
        self.graph_SP.showGrid(x = True, y = True, alpha = 0.8)
        self.graph_SK.showGrid(x = True, y = True, alpha = 0.8)

        #self.graph_T.setLabel("bottom", "t", units = "x", **labelStyle)
        #self.graph_EK.setLabel("bottom", "t", units = "x", **labelStyle)
        #self.graph_SP.setLabel("bottom", "t", units = "x", **labelStyle)
        #self.graph_SK.setLabel("bottom", "t", units = "x", **labelStyle)

        self.graph_T.setLabel("left", "T", units = "y", **labelStyle)
        self.graph_EK.setLabel("left", "EK", units = "y", **labelStyle)
        self.graph_SP.setLabel("left", "SP", units = "y", **labelStyle)
        self.graph_SK.setLabel("left", "SK", units = "y", **labelStyle)

        dat_layout.addWidget(self.graph_T)
        dat_layout.addWidget(self.graph_EK)
        dat_layout.addWidget(self.graph_SP)
        dat_layout.addWidget(self.graph_SK)



    def open_vasprun(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Open File", "data", "XML Files (*.xml)")

        if file_name:
            start_time = time.time()
            self.xml.parse(file_name)
            print("EXEC TIME:", time.time() - start_time, "seconds")
            print("----------------------------------------")
            self.plot_vasprun()
            XML.num += 1
            DAT.POTIM = self.xml.POTIM
            self.xml_transition = self.xml.block



    def open_dat(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Open File", "data", "DAT Files (*.dat)")

        if file_name:
            self.dat.add(file_name)
            self.plot_dat()
            DAT.num += 1
            self.dat_transition = self.dat.block



    def toggle_dat(self):
        if self.dat_visible:
            self.close_anim.start()
            self.dat_visible = False

        else:
            self.open_anim.start()
            self.dat_visible = True



    def clear_vasp(self):
        DAT.POTIM = 1
        self.xml = XML()
        self.xml_transition = 1
        self.graph_energy.clear()
        self.graph_entrp.clear()
        self.graph_total.clear()



    def clear_dat(self):
        self.dat = DAT()
        self.dat_transition = 1
        self.graph_T.clear()
        self.graph_EK.clear()
        self.graph_SP.clear()
        self.graph_SK.clear()



    def save(self):
        file_name, _ = QFileDialog.getSaveFileName(self, "Save File", "data/data.txt", "Text Files (*.txt);;All Files (*)")

        if file_name:
            self.xml.save(file_name)
            print("DATA SAVED:", file_name)



    def movie(self):
        pass



    def plot_vasprun(self):
        t = []
        e_fr_energy = []
        e_wo_entrp = []
        total = []

        for i in range(self.xml_transition, self.xml.block):
            t.append(self.xml.t[i])
            e_fr_energy.append(self.xml.e_fr_energy[i])
            e_wo_entrp.append(self.xml.e_wo_entrp[i])
            total.append(self.xml.total[i])

        color = self.colors[XML.num%len(self.colors)]

        min_t = min(self.xml.t.values())
        max_t = max(self.xml.t.values())

        self.graph_energy.plot(t, e_fr_energy, pen = pg.mkPen(color, width = 2))
        self.graph_energy.setXRange(min_t, max_t, padding=0)
        self.graph_energy.setYRange(min(self.xml.e_fr_energy.values()), max(self.xml.e_fr_energy.values()), padding=0)
        
        self.graph_entrp.plot(t, e_wo_entrp, pen = pg.mkPen(color, width = 2))
        self.graph_entrp.setXRange(min_t, max_t, padding=0)
        self.graph_entrp.setYRange(min(self.xml.e_wo_entrp.values()), max(self.xml.e_wo_entrp.values()), padding=0)
        
        self.graph_total.plot(t, total, pen = pg.mkPen(color, width = 2))
        self.graph_total.setXRange(min_t, max_t, padding=0)
        self.graph_total.setYRange(min(self.xml.total.values()), max(self.xml.total.values()), padding=0)



    def plot_dat(self):
        color = self.colors[DAT.num%len(self.colors)]

        t = self.dat.t[self.dat_transition-1:self.dat.block]
        T = self.dat.T[self.dat_transition-1:self.dat.block]
        EK = self.dat.EK[self.dat_transition-1:self.dat.block]
        SP = self.dat.SP[self.dat_transition-1:self.dat.block]
        SK = self.dat.SK[self.dat_transition-1:self.dat.block]

        min_t = min(self.dat.t)
        max_t = max(self.dat.t)

        self.graph_T.plot(t, T, pen = pg.mkPen(color, width = 2))
        self.graph_T.setXRange(min_t, max_t, padding=0)
        self.graph_T.setYRange(min(self.dat.T), max(self.dat.T), padding=0)

        self.graph_EK.plot(t, EK, pen = pg.mkPen(color, width = 2))
        self.graph_EK.setXRange(min_t, max_t, padding=0)
        self.graph_EK.setYRange(min(self.dat.EK), max(self.dat.EK), padding=0)

        self.graph_SP.plot(t, SP, pen = pg.mkPen(color, width = 2))
        self.graph_SP.setXRange(min_t, max_t, padding=0)
        self.graph_SP.setYRange(min(self.dat.SP), max(self.dat.SP), padding=0)

        self.graph_SK.plot(t, SK, pen = pg.mkPen(color, width = 2))
        self.graph_SK.setXRange(min_t, max_t, padding=0)
        self.graph_SK.setYRange(min(self.dat.SK), max(self.dat.SK), padding=0) 




app = QApplication(sys.argv)
window = Window()
sys.exit(app.exec())

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
        self.num_xml = 0
        self.num_dat = 0
        self.block_trs = 1
        self.block_dat = 1
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

        self.graph00 = pg.PlotWidget()
        self.graph01 = pg.PlotWidget()
        self.graph10 = pg.PlotWidget()
        self.graph11 = pg.PlotWidget()

        self.graph00.showGrid(x = True, y = True, alpha = 0.8)
        self.graph01.showGrid(x = True, y = True, alpha = 0.8)
        self.graph10.showGrid(x = True, y = True, alpha = 0.8)
        self.graph11.showGrid(x = True, y = True, alpha = 0.8)

        self.graph00.setLabel("bottom", "i * POTIM", units = "x", **labelStyle)
        self.graph01.setLabel("bottom", "i * POTIM", units = "x", **labelStyle)
        self.graph10.setLabel("bottom", "i * POTIM", units = "x", **labelStyle)
        self.graph11.setLabel("bottom", "", units = "x", **labelStyle)

        self.graph00.setLabel("left", "energy", units = "y", **labelStyle)
        self.graph01.setLabel("left", "entrp", units = "y", **labelStyle)
        self.graph10.setLabel("left", "total", units = "y", **labelStyle)
        self.graph11.setLabel("left", "", units = "y", **labelStyle)

        main_layout.addWidget(self.graph00, 0, 0)
        main_layout.addWidget(self.graph01, 0, 1)
        main_layout.addWidget(self.graph10, 1, 0)
        main_layout.addWidget(self.graph11, 1, 1)


        #GRAPHS DAT
        self.graph0 = pg.PlotWidget()
        self.graph1 = pg.PlotWidget()
        self.graph2 = pg.PlotWidget()
        self.graph3 = pg.PlotWidget()

        self.graph0.showGrid(x = True, y = True, alpha = 0.8)
        self.graph1.showGrid(x = True, y = True, alpha = 0.8)
        self.graph2.showGrid(x = True, y = True, alpha = 0.8)
        self.graph3.showGrid(x = True, y = True, alpha = 0.8)

        #self.graph0.setLabel("bottom", "i * POTIM", units = "x", **labelStyle)
        #self.graph1.setLabel("bottom", "i * POTIM", units = "x", **labelStyle)
        #self.graph2.setLabel("bottom", "i * POTIM", units = "x", **labelStyle)
        #self.graph3.setLabel("bottom", "", units = "x", **labelStyle)

        self.graph0.setLabel("left", "T", units = "y", **labelStyle)
        self.graph1.setLabel("left", "EK", units = "y", **labelStyle)
        self.graph2.setLabel("left", "SP", units = "y", **labelStyle)
        self.graph3.setLabel("left", "SK", units = "y", **labelStyle)

        dat_layout.addWidget(self.graph0)
        dat_layout.addWidget(self.graph1)
        dat_layout.addWidget(self.graph2)
        dat_layout.addWidget(self.graph3)



    def open_vasprun(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Open File", "data", "XML Files (*.xml)")

        if file_name:
            start_time = time.time()
            self.xml.parse(file_name)
            print("EXEC TIME:", time.time() - start_time, "seconds")
            print("----------------------------------------")
            self.plot_vasprun()
            self.num_xml += 1
            self.block_trs = self.xml.block



    def open_dat(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Open File", "data", "DAT Files (*.dat)")

        if file_name:
            try:
                f = open(file_name, "r")
                t = []
                T = []
                EK = []
                SP = []
                SK = []

                for l in f:
                    line = l.split()
                    if line:
                        t.append(self.block_dat * self.xml.POTIM)
                        T.append(float(line[2]))
                        EK.append(float(line[10]))
                        SP.append(float(line[12]))
                        SK.append(float(line[14]))
                        self.block_dat += 1
                
                color = self.colors[self.num_dat%len(self.colors)]
                self.graph0.plot(t, T, pen = pg.mkPen(color, width = 2))
                self.graph1.plot(t, EK, pen = pg.mkPen(color, width = 2))
                self.graph2.plot(t, SP, pen = pg.mkPen(color, width = 2))
                self.graph3.plot(t, SK, pen = pg.mkPen(color, width = 2))
                self.num_dat += 1


            except:
                raise Exception("Invalid DAT file")
                #CLEAR?



    def toggle_dat(self):
        if self.dat_visible:
            self.close_anim.start()
            self.dat_visible = False

        else:
            self.open_anim.start()
            self.dat_visible = True



    def clear_vasp(self):
        self.xml = XML()
        self.num_xml = 0
        self.block_trs = 1
        self.graph00.clear()
        self.graph01.clear()
        self.graph10.clear()



    def clear_dat(self):
        self.num_dat = 0
        self.block_dat = 1
        self.graph0.clear()
        self.graph1.clear()
        self.graph2.clear()
        self.graph3.clear()



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

        for i in range(self.block_trs, self.xml.block):
            t.append(self.xml.t[i])
            e_fr_energy.append(self.xml.energy[i][0])
            e_wo_entrp.append(self.xml.energy[i][1])
            total.append(self.xml.energy[i][2])

        color = self.colors[self.num_xml%len(self.colors)]
        self.graph00.plot(t, e_fr_energy, pen = pg.mkPen(color, width = 2))
        self.graph01.plot(t, e_wo_entrp, pen = pg.mkPen(color, width = 2))
        self.graph10.plot(t, total, pen = pg.mkPen(color, width = 2))




app = QApplication(sys.argv)
window = Window()
sys.exit(app.exec())

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
        self.setWindowTitle("Atom Viewer")
        self.setGeometry(100, 50, 1600, 980)
        self.UI()
        self.show()

        #VARIABLES
        self.DAT = 0
        self.DAT_I = 1
        self.DAT_POTIM = 1
        self.colors = ["r", "y", "g", "c", "b", "m"]



    def UI(self):
        #LAYOUTS
        main_widget = QWidget()
        main_layout = QGridLayout()
        main_widget.setLayout(main_layout)
        self.setCentralWidget(main_widget)

        control_widget = QWidget()
        control_layout = QHBoxLayout()
        control_widget.setLayout(control_layout)
        main_layout.addWidget(control_widget, 2, 0)

        self.dialog = QDialog(self)
        dialog_layout = QVBoxLayout()
        self.dialog.setLayout(dialog_layout)
        self.dialog.setWindowTitle("DAT Editor")
        self.dialog.setGeometry(1000, 50, 700, 980)

        dialog_control_widget = QWidget()
        dialog_control_layout = QHBoxLayout()
        dialog_control_widget.setLayout(dialog_control_layout)
        dialog_layout.addWidget(dialog_control_widget)


        #BUTTONS
        open_vasprun_button = QPushButton("Open Vasprun", self)
        open_vasprun_button.clicked.connect(self.open_vasprun)
        control_layout.addWidget(open_vasprun_button)

        dat_button = QPushButton("DAT", self)
        dat_button.clicked.connect(self.dialog.show)
        control_layout.addWidget(dat_button)

        open_dat_button = QPushButton("Open Dat", self)
        open_dat_button.clicked.connect(self.open_dat)
        dialog_control_layout.addWidget(open_dat_button)

        clear_button = QPushButton("Clear", self)
        clear_button.clicked.connect(self.clear)
        dialog_control_layout.addWidget(clear_button)

        save_button = QPushButton("Save", self)
        save_button.clicked.connect(self.save)
        control_layout.addWidget(save_button)


        #GRAPHS VASPRUN
        font = QFontDatabase.addApplicationFont("data/cmunti.ttf")
        labelStyle = {'color': '#969696', 'font-size': '18px', 'font-family': 'CMU Serif'}

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

        self.graph0.setLabel("left", "energy", units = "y", **labelStyle)
        self.graph1.setLabel("left", "entrp", units = "y", **labelStyle)
        self.graph2.setLabel("left", "total", units = "y", **labelStyle)
        self.graph3.setLabel("left", "", units = "y", **labelStyle)

        dialog_layout.addWidget(self.graph0)
        dialog_layout.addWidget(self.graph1)
        dialog_layout.addWidget(self.graph2)
        dialog_layout.addWidget(self.graph3)



    def open_vasprun(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Open File", "data", "XML Files (*.xml)")

        if file_name:
            start_time = time.time()
            self.xml = XML()
            self.xml.parse(file_name)
            print("EXEC TIME:", time.time() - start_time, "seconds")
            self.DAT_POTIM = self.xml.POTIM
            self.plot_vasprun()



    def open_dat(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Open File", "data", "DAT Files (*.dat)")

        if file_name:
            try:
                f = open(file_name, "r")
                X = []
                T = []
                EK = []
                SP = []
                SK = []

                for l in f:
                    line = l.split()
                    if line:
                        X.append(self.DAT_I*self.DAT_POTIM)
                        T.append(float(line[2]))
                        EK.append(float(line[10]))
                        SP.append(float(line[12]))
                        SK.append(float(line[14]))
                        self.DAT_I += 1
                
                self.graph0.plot(X, T, pen = pg.mkPen(self.colors[self.DAT%len(self.colors)], width = 2))
                self.graph1.plot(X, EK, pen = pg.mkPen(self.colors[self.DAT%len(self.colors)], width = 2))
                self.graph2.plot(X, SP, pen = pg.mkPen(self.colors[self.DAT%len(self.colors)], width = 2))
                self.graph3.plot(X, SK, pen = pg.mkPen(self.colors[self.DAT%len(self.colors)], width = 2))
                self.DAT += 1


            except:
                raise Exception("Invalid DAT file")
                #CLEAR?



    def clear(self):
        self.DAT = 0
        self.DAT_I = 1
        self.graph0.clear()
        self.graph1.clear()
        self.graph2.clear()
        self.graph3.clear()



    def save(self):
        file_name, _ = QFileDialog.getSaveFileName(self, "Save File", "data/data.txt", "Text Files (*.txt);;All Files (*)")

        if file_name:
            try:
                self.xml.save(file_name)
                print("DATA SAVED:", file_name)
            
            except AttributeError:
                raise Exception("Vasprun not found")



    def plot_vasprun(self):
        e_fr_energy = []
        e_wo_entrp = []
        total = []

        for i in range(1, self.xml.i):
            e_fr_energy.append(self.xml.energy[i][0])
            e_wo_entrp.append(self.xml.energy[i][1])
            total.append(self.xml.energy[i][2])

        self.graph00.clear()
        self.graph01.clear()
        self.graph10.clear()
        self.graph00.plot(self.xml.t, e_fr_energy, pen = pg.mkPen(width = 2))
        self.graph01.plot(self.xml.t, e_wo_entrp, pen = pg.mkPen(width = 2))
        self.graph10.plot(self.xml.t, total, pen = pg.mkPen(width = 2))




app = QApplication(sys.argv)
stylesheet = open("data/style.css", "r").read()
app.setStyleSheet(stylesheet)
window = Window()
sys.exit(app.exec())

from xml_parser import Xml
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import pyqtgraph as pg
import time
import sys



class Window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Atom Visualisation")
        self.setGeometry(100, 100, 1280, 720)
        self.UI()
        self.show()

    
    def UI(self):
        main_widget = QWidget()
        main_layout = QGridLayout()
        main_widget.setLayout(main_layout)

        control_widget = QWidget()
        control_layout = QHBoxLayout()
        control_widget.setLayout(control_layout)
        main_layout.addWidget(control_widget, 2, 0)

        open_button = QPushButton("Open", self)
        open_button.clicked.connect(self.open)
        control_layout.addWidget(open_button)

        save_button = QPushButton("Save", self)
        save_button.clicked.connect(self.save)
        control_layout.addWidget(save_button)

        self.graph00 = pg.PlotWidget()
        self.graph01 = pg.PlotWidget()
        self.graph10 = pg.PlotWidget()
        self.graph11 = pg.PlotWidget()
        main_layout.addWidget(self.graph00, 0, 0)
        main_layout.addWidget(self.graph01, 0, 1)
        main_layout.addWidget(self.graph10, 1, 0)
        main_layout.addWidget(self.graph11, 1, 1)

        self.setCentralWidget(main_widget)
    

    def open(self):
        file_name = QFileDialog.getOpenFileName(self, "OpenFile")[0]

        try:
            start_time = time.time()
            xml = Xml()
            xml.parse(file_name)
            print("--- %s seconds ---" % (time.time() - start_time))

            print(xml.POTIM)

        except:
            raise Exception("Invalid file")


    def save(self):
        pass




app = QApplication(sys.argv)
window = Window()
sys.exit(app.exec())

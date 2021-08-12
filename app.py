from numpy import empty
from dat import DAT
from xml_parser import XML
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import pyqtgraph as pg
import time
import sys




class Window2(QMainWindow):
    def __init__(self):
        super().__init__()
        stylesheet = open("resources/style.css", "r").read()
        self.setWindowIcon(QIcon("resources/icon.png"))
        self.setWindowTitle("Details")
        self.setGeometry(100, 100, 400, 800)
        self.setStyleSheet(stylesheet)




class Window(QMainWindow):
    def __init__(self):
        super().__init__()
        stylesheet = open("resources/style.css", "r").read()
        self.setWindowIcon(QIcon("resources/icon.png"))
        self.setWindowTitle("Atom Viewer")
        self.setGeometry(50, 50, 1600, 980)
        self.setStyleSheet(stylesheet)
        self.UI()
        self.show()

        #VARIABLES
        self.xml = XML()
        self.dat = DAT()
        self.empty_vasp = True
        self.xml_transition = 1
        self.dat_transition = 1
        self.dat_visible = False
        self.colors = ["r", "y", "g", "c", "b", "m"]



    def UI(self):
        #MAIN WINDOW
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

        info_widget = QWidget()
        info_layout = QGridLayout(info_widget)
        info_layout.setSpacing(0)
        main_layout.addWidget(info_widget, 1, 0)

        param_widget = QWidget()
        param_layout = QHBoxLayout(param_widget)
        param_layout.setSpacing(0)
        main_layout.addWidget(param_widget, 2, 0, 2, 1)


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

        details_button = QPushButton("Details", self)
        details_button.clicked.connect(self.toggle_details)
        control_top_layout.addWidget(details_button)

        save_button = QPushButton("Save", self)
        save_button.clicked.connect(self.save)
        control_bot_layout.addWidget(save_button)

        movie_button = QPushButton("Movie", self)
        movie_button.clicked.connect(self.movie)
        control_bot_layout.addWidget(movie_button)

        json_button = QPushButton("Json", self)
        json_button.clicked.connect(self.json)
        control_bot_layout.addWidget(json_button)

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
        font = QFontDatabase.addApplicationFont("resources/lmromanb.ttf")
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


        #INFOS
        line = QLabel("")
        info_layout.addWidget(line, 0, 0)
        vasp_title = QLabel("Vasprun")
        dat_title = QLabel("Dat")
        l_num = QLabel("Files")
        l_last = QLabel("Last")
        l_block = QLabel("Blocks")
        info_layout.addWidget(vasp_title, 0, 1, 1, 10)
        info_layout.addWidget(dat_title, 0, 11, 1, 10)
        info_layout.addWidget(l_num, 1, 0, 1, 1)
        info_layout.addWidget(l_last, 2, 0, 1, 1)
        info_layout.addWidget(l_block, 3, 0, 1, 1)

        self.vasp_num = QLabel("")
        self.dat_num = QLabel("")
        self.vasp_block = QLabel("")
        self.dat_block = QLabel("")
        self.vasp_last = QLabel("")
        self.dat_last = QLabel("")
        info_layout.addWidget(self.vasp_num, 1, 1, 1, 10)
        info_layout.addWidget(self.dat_num, 1, 11, 1, 10)
        info_layout.addWidget(self.vasp_last, 2, 1, 1, 10)
        info_layout.addWidget(self.dat_last, 2, 11, 1, 10)
        info_layout.addWidget(self.vasp_block, 3, 1, 1, 10)
        info_layout.addWidget(self.dat_block, 3, 11, 1, 10)


        #PARAMS
        num_atoms_title = QLabel("Atoms:  ")
        num_types_title = QLabel("Types:  ")
        potim_title = QLabel("POTIM:  ")
        num_atoms_title.setObjectName("param")
        num_types_title.setObjectName("param")
        potim_title.setObjectName("param")

        self.num_atoms = QLabel("")
        self.num_types = QLabel("")
        self.potim = QLabel("")
        self.num_atoms.setObjectName("num")
        self.num_types.setObjectName("num")
        self.potim.setObjectName("num")

        param_layout.addWidget(num_atoms_title)
        param_layout.addWidget(self.num_atoms)
        param_layout.addWidget(num_types_title)
        param_layout.addWidget(self.num_types)
        param_layout.addWidget(potim_title)
        param_layout.addWidget(self.potim)



        #DETAILS WINDOW
        self.w = Window2()
        details_widget = QWidget()
        details_layout = QVBoxLayout(details_widget)

        tabs = QTabWidget()
        self.input_widget = QLineEdit(self)
        self.table_basis = QTableWidget()
        self.table_position = QTableWidget()
        self.table_force = QTableWidget()
        self.table_stress = QTableWidget()
        
        tabs.addTab(self.table_basis, "Basis")
        tabs.addTab(self.table_position, "Position")
        tabs.addTab(self.table_force, "Force")
        tabs.addTab(self.table_stress, "Stress")

        ok_button = QPushButton("Ok", self)
        ok_button.clicked.connect(self.update_details)

        details_layout.addWidget(self.input_widget)
        details_layout.addWidget(ok_button)
        details_layout.addWidget(tabs)
        self.w.setCentralWidget(details_widget)

        self.table_basis.setColumnCount(3)
        self.table_basis.setRowCount(3)
        header = self.table_basis.horizontalHeader()
        self.table_basis.setHorizontalHeaderLabels(["", "", ""])
        self.table_basis.setVerticalHeaderLabels(["", "", ""])
        header.setSectionResizeMode(0, QHeaderView.Stretch)
        header.setSectionResizeMode(1, QHeaderView.Stretch)
        header.setSectionResizeMode(2, QHeaderView.Stretch)

        self.table_stress.setColumnCount(3)
        self.table_stress.setRowCount(3)
        header = self.table_stress.horizontalHeader()
        self.table_stress.setHorizontalHeaderLabels(["", "", ""])
        self.table_stress.setVerticalHeaderLabels(["", "", ""])
        header.setSectionResizeMode(0, QHeaderView.Stretch)
        header.setSectionResizeMode(1, QHeaderView.Stretch)
        header.setSectionResizeMode(2, QHeaderView.Stretch)

        self.table_position.setColumnCount(3)
        self.table_position.setRowCount(0)
        header = self.table_position.horizontalHeader() 
        self.table_position.setHorizontalHeaderLabels(["x", "y", "z"])
        header.setSectionResizeMode(0, QHeaderView.Stretch)
        header.setSectionResizeMode(1, QHeaderView.Stretch)
        header.setSectionResizeMode(2, QHeaderView.Stretch)

        self.table_force.setColumnCount(3)
        self.table_force.setRowCount(0)
        header = self.table_force.horizontalHeader()
        self.table_force.setHorizontalHeaderLabels(["x", "y", "z"])
        header.setSectionResizeMode(0, QHeaderView.Stretch)
        header.setSectionResizeMode(1, QHeaderView.Stretch)
        header.setSectionResizeMode(2, QHeaderView.Stretch)



    def open_vasprun(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Open File", "", "XML Files (*.xml)")

        if file_name:
            start_time = time.time()
            self.xml.parse(file_name)
            print("EXEC TIME:", time.time() - start_time, "seconds")
            print("----------------------------------------")
            self.plot_vasprun()
            self.empty_vasp = False
            self.update_details()
            XML.num += 1
            DAT.POTIM = self.xml.POTIM
            tmp = self.xml.block - 1 if self.xml_transition == 1 else self.xml.block - 1 - self.xml_transition
            self.vasp_num.setText(str(XML.num))
            self.vasp_last.setText(str(tmp))
            self.vasp_block.setText(str(self.xml.block - 1))
            self.num_atoms.setText(str(self.xml.num_atoms))
            self.num_types.setText(str(self.xml.num_types))
            self.potim.setText(str(self.xml.POTIM))
            self.xml_transition = self.xml.block - 1

            if self.xml.NSW != tmp:
                self.vasp_last.setStyleSheet("color: red")
                self.vasp_block.setStyleSheet("color: red")

            else:
                self.vasp_last.setStyleSheet("")



    def open_dat(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Open File", "", "DAT Files (*.dat)")

        if file_name:
            self.dat.add(file_name)
            self.plot_dat()
            DAT.num += 1
            tmp = 0 if self.dat_transition == 1 else self.dat_transition
            self.dat_num.setText(str(DAT.num))
            self.dat_block.setText(str(self.dat.block - 1))
            self.dat_last.setText(str(self.dat.block - 1 - tmp))
            self.dat_transition = self.dat.block - 1



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
        self.empty_vasp = True
        self.vasp_num.setText("")
        self.vasp_last.setText("")
        self.vasp_block.setText("")
        self.num_atoms.setText("")
        self.num_types.setText("")
        self.potim.setText("")
        self.vasp_last.setStyleSheet("")
        self.vasp_block.setStyleSheet("")
        self.graph_energy.clear()
        self.graph_entrp.clear()
        self.graph_total.clear()



    def toggle_details(self):
        if self.w.isVisible():
            self.w.hide()
        
        else:
            self.w.show()



    def clear_dat(self):
        self.dat = DAT()
        self.dat_transition = 1
        self.dat_num.setText("")
        self.dat_last.setText("")
        self.dat_block.setText("")
        self.graph_T.clear()
        self.graph_EK.clear()
        self.graph_SP.clear()
        self.graph_SK.clear()



    def save(self):
        file_name, _ = QFileDialog.getSaveFileName(self, "Save File", "data.dat", "Dat Files (*.dat);;All Files (*)")

        if file_name:
            self.xml.save(file_name)
            print("DATA SAVED:", file_name)



    def movie(self):
        file_name, _ = QFileDialog.getSaveFileName(self, "Save File", "movie.xyz", "Dat Files (*.xyz);;All Files (*)")

        if file_name:
            self.xml.movie(file_name)
            print("DATA SAVED:", file_name)



    def json(self):
        file_name, _ = QFileDialog.getSaveFileName(self, "Save File", "data.json", "Dat Files (*.json);;All Files (*)")

        if file_name:
            self.xml.save_json(file_name)
            print("DATA SAVED:", file_name)



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



    def update_details(self):
        if self.empty_vasp == False:
            
            if self.input_widget.text().isnumeric():
                self.input_widget.setStyleSheet("color: #000000;")
                num_block = int(self.input_widget.text())

                if 1 <= num_block < self.xml.block:
                    self.input_widget.setStyleSheet("color: #000000;")
                    self.table_position.setRowCount(self.xml.num_atoms)
                    self.table_force.setRowCount(self.xml.num_atoms)

                    #BASIS/STRESS
                    for i in range(3):
                        for j in range(3):
                            self.table_basis.setItem(i, j, QTableWidgetItem(str(self.xml.basis[num_block][i][j])))
                            self.table_stress.setItem(i, j, QTableWidgetItem(str(self.xml.stress[num_block][i][j])))

                    #POSITION/FORCE
                    for i in range(self.xml.num_atoms):
                        self.table_position.setItem(i, 0, QTableWidgetItem(str(self.xml.position[num_block][i][0])))
                        self.table_position.setItem(i, 1, QTableWidgetItem(str(self.xml.position[num_block][i][1])))
                        self.table_position.setItem(i, 2, QTableWidgetItem(str(self.xml.position[num_block][i][2])))

                        self.table_force.setItem(i, 0, QTableWidgetItem(str(self.xml.force[num_block][i][0])))
                        self.table_force.setItem(i, 1, QTableWidgetItem(str(self.xml.force[num_block][i][1])))
                        self.table_force.setItem(i, 2, QTableWidgetItem(str(self.xml.force[num_block][i][2])))

                    print("UPDATED")
                
                else:
                    self.input_widget.setStyleSheet("color: #d10003;")
            
            else:
                self.input_widget.setStyleSheet("color: #d10003;")




if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = Window()
    sys.exit(app.exec())

import xml.etree.ElementTree as ET

class XML:
    def parse(self, file):
        xml = ET.parse(file)

        self.date = xml.find("./generator/i[@name='date']").text
        self.time = xml.find("./generator/i[@name='time']").text
        self.num_atoms = int(xml.find("./atominfo/atoms").text.strip())
        self.num_types = int(xml.find("./atominfo/types").text.strip())
        self.divisions = xml.find("./kpoints/generation/v[@name='divisions']").text ###
        self.PREC = xml.find("./incar/i[@name='PREC']").text.strip()
        self.ALGO = xml.find("./incar/i[@name='ALGO']").text.strip()
        self.NSIM = int(xml.find("./incar/i[@name='NSIM']").text.strip())
        self.ENCUT = float(xml.find("./incar/i[@name='ENCUT']").text.strip())
        self.POTIM = float(xml.find("./incar/i[@name='POTIM']").text.strip())
        self.TEBEG = float(xml.find("./incar/i[@name='TEBEG']").text.strip())
        self.SMASS = float(xml.find("./incar/i[@name='SMASS']").text.strip())
        self.LREAL = xml.find("./incar/i[@name='LREAL']").text.strip()


        #INITIALPOS
        self.initial_basis = []
        for v in xml.find("./structure[@name='initialpos']/crystal/varray[@name='basis']"):
            self.initial_basis.append([float(k) for k in v.text.split()])

        self.initial_positions = []
        for v in xml.find("./structure[@name='initialpos']/varray[@name='positions']"):
            self.initial_positions.append([float(k) for k in v.text.split()])

        self.initial_volume = float(xml.find("./structure[@name='initialpos']/crystal/i[@name='volume']").text.strip())


        #FINALPOS
        self.final_basis = []
        for v in xml.find("./structure[@name='finalpos']/crystal/varray[@name='basis']"):
            self.final_basis.append([float(k) for k in v.text.split()])
        
        self.final_positions = []
        for v in xml.find("./structure[@name='finalpos']/varray[@name='positions']"):
            self.final_positions.append([float(k) for k in v.text.split()])
        
        self.final_volume = float(xml.find("./structure[@name='finalpos']/crystal/i[@name='volume']").text.strip())


        #CALCULATIONS
        self.basis = dict()
        self.position = dict()
        self.force = dict()
        self.stress = dict()
        self.volume = dict()
        self.energy = dict()
        self.t = []
        i = 1

        for calculation in xml.findall("./calculation"):
            #BASIS
            i_basis = []
            for v in calculation.find("./structure/crystal/varray[@name='basis']"):
                i_basis.append([float(k) for k in v.text.split()])
            self.basis[i] = i_basis

            #POSITION
            i_position = []
            for v in calculation.find("./structure/varray[@name='positions']"):
                i_position.append([float(k) for k in v.text.split()])
            self.position[i] = i_position

            #FORCE
            i_force = []
            for v in calculation.find("./varray[@name='forces']"):
                i_force.append([float(k) for k in v.text.split()])
            self.force[i] = i_force

            #STRESS
            i_stress = []
            for v in calculation.find("./varray[@name='stress']"):
                i_stress.append([float(k) for k in v.text.split()])
            self.stress[i] = i_stress

            #VOLUME
            self.volume[i] = float(calculation.find("./structure/crystal/i[@name='volume']").text.strip())

            #ENERGY
            self.energy[i] = [float(calculation.find("./energy/i[@name='e_fr_energy']").text.strip()),
                            float(calculation.find("./energy/i[@name='e_wo_entrp']").text.strip()),
                            float(calculation.find("./energy/i[@name='total']").text.strip())]

            self.t.append(i*self.POTIM)
            i += 1
        
        self.num_blocs = len(self.position)


    def save(self, file):
        f = open(file,"a")

        for i in range(1, self.num_blocs+1):
            f.write(str(self.energy[i][0]))
            f.write(":")
            f.write(str(self.num_atoms))
            f.write(":")
            f.write(str(self.num_types))
            f.write(":")
            f.write(str(self.position[i]))
            f.write("\n")
        
        f.close()

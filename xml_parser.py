import xml.etree.ElementTree as ET




class XML:
    def __init__(self):
        self.basis = dict()
        self.position = dict()
        self.force = dict()
        self.stress = dict()
        self.volume = dict()
        self.energy = dict()
        self.t = []
        self.i = 1



    def parse(self, file):
        try:
            for _, element in ET.iterparse(file):
                #GENERATOR
                if element.tag == "generator":
                    self.date = element.find("i[@name='date']").text
                    self.time = element.find("i[@name='time']").text
                    element.clear()

                #ATOMINFO
                elif element.tag == "atominfo":
                    self.num_atoms = int(element.find("atoms").text)
                    self.num_types = int(element.find("types").text)
                    element.clear()

                #KPOINTS
                elif element.tag == "kpoints":
                    self.divisions = element.find("generation/v[@name='divisions']").text.split()
                    element.clear()

                #INCAR
                elif element.tag == "incar":
                    self.PREC = element.find("i[@name='PREC']").text
                    self.ALGO = element.find("i[@name='ALGO']").text
                    self.NSIM = int(element.find("i[@name='NSIM']").text)
                    self.NSW = int(element.find("i[@name='NSW']").text)
                    self.ENCUT = float(element.find("i[@name='ENCUT']").text)
                    self.POTIM = float(element.find("i[@name='POTIM']").text)
                    self.TEBEG = float(element.find("i[@name='TEBEG']").text)
                    self.SMASS = float(element.find("i[@name='SMASS']").text)
                    self.LREAL = element.find("i[@name='LREAL']").text
                    element.clear()


                #CALCULATION
                elif element.tag == "calculation":
                    #BASIS
                    i_basis = []
                    for v in element.find("structure/crystal/varray[@name='basis']"):
                        i_basis.append([float(k) for k in v.text.split()])
                    self.basis[self.i] = i_basis

                    #POSITION
                    i_position = []
                    for v in element.find("structure/varray[@name='positions']"):
                        i_position.append([float(k) for k in v.text.split()])
                    self.position[self.i] = i_position

                    #FORCE
                    i_force = []
                    for v in element.find("varray[@name='forces']"):
                        i_force.append([float(k) for k in v.text.split()])
                    self.force[self.i] = i_force

                    #STRESS
                    i_stress = []
                    for v in element.find("varray[@name='stress']"):
                        i_stress.append([float(k) for k in v.text.split()])
                    self.stress[self.i] = i_stress

                    #VOLUME
                    self.volume[self.i] = float(element.find("structure/crystal/i[@name='volume']").text)

                    #ENERGY
                    self.energy[self.i] = [float(element.find("energy/i[@name='e_fr_energy']").text),
                                           float(element.find("energy/i[@name='e_wo_entrp']").text),
                                           float(element.find("energy/i[@name='total']").text)]

                    #VARIABLES
                    self.t.append(self.i*self.POTIM)
                    self.i += 1
                    element.clear()


                #INITIALPOS/FINALPOS
                elif element.tag == "structure" and "name" in element.attrib and (element.attrib["name"] == "initialpos" or element.attrib["name"] == "finalpos"):
                    tmp_basis = []
                    tmp_position = []
                    
                    if element.attrib["name"] == "initialpos":
                        tmp_i = 0

                    else:
                        tmp_i = self.i

                    #BASIS
                    for v in element.find("crystal/varray[@name='basis']"):
                        tmp_basis.append([float(k) for k in v.text.split()])

                    #POSITION
                    for v in element.find("varray[@name='positions']"):
                        tmp_position.append([float(k) for k in v.text.split()])

                    #VOLUME
                    self.volume[tmp_i] = float(element.find("crystal/i[@name='volume']").text)
                    
                    self.basis[tmp_i] = tmp_basis
                    self.position[tmp_i] = tmp_position
                    element.clear()
                    
            print("SUCCESS: Vasprun loaded")


        except ET.ParseError:
            print("WARNING: Incomplete/invalid XML file")



    def save(self, file):
        f = open(file, "a")

        for i in range(1, self.i):
            f.write(f"{self.energy[i][0]}:{self.num_atoms}:{self.num_types}:{self.position[i]}\n")
        
        f.close()

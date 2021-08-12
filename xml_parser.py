import xml.etree.ElementTree as ET
import json




class XML:
    num = 0

    def __init__(self):
        self.basis = dict()
        self.position = dict()
        self.force = dict()
        self.stress = dict()
        self.volume = dict()
        self.e_fr_energy = dict()
        self.e_wo_entrp = dict()
        self.total = dict()
        self.t = dict()
        self.types = []
        self.block = 1
        XML.num = 0



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
                    for e in element.find("array[@name='atoms']/set"):
                        self.types.append(e.find("c").text)
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
                    tmp_basis = []
                    for v in element.find("structure/crystal/varray[@name='basis']"):
                        tmp_basis.append([float(k) for k in v.text.split()])
                    self.basis[self.block] = tmp_basis

                    #POSITION
                    tmp_position = []
                    for v in element.find("structure/varray[@name='positions']"):
                        tmp_position.append([float(k) for k in v.text.split()])
                    self.position[self.block] = tmp_position

                    #FORCE
                    tmp_force = []
                    for v in element.find("varray[@name='forces']"):
                        tmp_force.append([float(k) for k in v.text.split()])
                    self.force[self.block] = tmp_force

                    #STRESS
                    tmp_stress = []
                    for v in element.find("varray[@name='stress']"):
                        tmp_stress.append([float(k) for k in v.text.split()])
                    self.stress[self.block] = tmp_stress

                    #VOLUME
                    self.volume[self.block] = float(element.find("structure/crystal/i[@name='volume']").text)

                    #ENERGY
                    self.e_fr_energy[self.block] = float(element.find("energy/i[@name='e_fr_energy']").text)
                    self.e_wo_entrp[self.block] = float(element.find("energy/i[@name='e_wo_entrp']").text)
                    self.total[self.block] = float(element.find("energy/i[@name='total']").text)

                    #VARIABLES
                    self.t[self.block] = self.block * self.POTIM
                    self.block += 1
                    element.clear()


                #INITIALPOS/FINALPOS
                """
                elif element.tag == "structure" and "name" in element.attrib and (element.attrib["name"] == "initialpos" or element.attrib["name"] == "finalpos"):
                    tmp_basis = []
                    tmp_position = []
                    
                    if element.attrib["name"] == "initialpos":
                        tmp_block = 0

                    else:
                        tmp_block = self.block

                    #BASIS
                    for v in element.find("crystal/varray[@name='basis']"):
                        tmp_basis.append([float(k) for k in v.text.split()])

                    #POSITION
                    for v in element.find("varray[@name='positions']"):
                        tmp_position.append([float(k) for k in v.text.split()])

                    #VOLUME
                    self.volume[tmp_block] = float(element.find("crystal/i[@name='volume']").text)
                    
                    self.basis[tmp_block] = tmp_basis
                    self.position[tmp_block] = tmp_position
                    element.clear()
                """
                    
            print("SUCCESS: Vasprun loaded")


        except ET.ParseError:
            print(f"WARNING: Incomplete/invalid XML file ({file})")



    def save(self, file):
        f = open(file, "a")

        for i in range(1, self.block):
            f.write(f"{self.e_fr_energy[i]}:{self.num_atoms}:{self.num_types}:{self.position[i]}\n")
        
        f.close()



    def movie(self, file):
        f = open(file, "a")

        for i in range(1, self.block):
            f.write(f"{i}\n")
            f.write("\n")

            for j in range(len(self.position[i])):
                basis = self.basis[i]
                x, y, z = self.position[i][j]
                bx, by, bz = basis[0][0], basis[1][1], basis[2][2]
                f.write(f"{self.types[j]}   {x*bx}   {y*by}   {z*bz}\n")

            f.write("\n")
            f.write("\n")

        f.close()



    def save_json(self, file):
        data = []

        for i in range(1, self.block):
            istep = {}
            istep["energy"] = self.e_wo_entrp[i]
            istep["stress"] = self.stress[i]
            istep["volume"] = self.volume[i]
            
            atoms = []

            for j in range(self.num_atoms):
                atom_dict = {}
                atom_dict["id"] = j
                atom_dict["symbol"] = self.types[j]
                atom_dict["ordered_distances"] = {}
                atom_dict["forces"] = {}
                atoms.append(atom_dict)
            
            istep["atoms"] = atoms
            data.append(istep)
        
        f = open(file, "a")
        json.dump(data, f)
        f.close()

class DAT:
    num = 0
    POTIM = 1

    def __init__(self):
        self.t = []
        self.T = []
        self.EK = []
        self.SP = []
        self.SK = []
        self.block = 1
        self.POTIM = 1
        DAT.num = 0



    def add(self, file):
        try:
            f = open(file, "r")

            for l in f:
                line = l.split()
                if line:
                    self.t.append(self.block * DAT.POTIM)
                    self.T.append(float(line[2]))
                    self.EK.append(float(line[10]))
                    self.SP.append(float(line[12]))
                    self.SK.append(float(line[14]))
                    self.block += 1


        except:
            raise Exception("Invalid DAT file")

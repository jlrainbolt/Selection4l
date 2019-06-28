class LHEevent():
    
    def __init__(self):
        self.Particles = []
        
    def fillEvent(self, lheLines):
        # check that this is a good event
        #if lheLines[0].find("<event>") == -1 or lheLines[-1].find("</event>") == -1:
            #print "THIS IS NOT A LHE EVENT"
            #return 0
        for i in range(2,len(lheLines)-1):
            self.Particles.append(self.readParticle(lheLines[i]))
        return 1

    def readParticle(self, lheLine):
        dataIN = lheLine[:-1].split(" ")
        dataINgood = []
        for entry in dataIN:
            if entry != "": dataINgood.append(entry)
        return {'ID': int(dataINgood[0]),
                'mIdx': int(dataINgood[2])-1,
                'mIdx2': int(dataINgood[3])-1,
                'Px' : float(dataINgood[6]),
                'Py' : float(dataINgood[7]),
                'Pz' : float(dataINgood[8]),
                'E' : float(dataINgood[9]),
                'M' : float(dataINgood[10])}


class LHEfile():
    
    def __init__(self, fileINname):
        self.eventList = []
        self.fileINname = fileINname
        self.MaxEv = -99
        self.Model = "NONE"
        self.MU = 0
        self.WU = 0
        self.gUbe = 0
        self.gUbmu = 0
        self.xsec = 0

    def setMax(self, maxVal):
        self.Max = maxVal

    def readInfo(self):
        fileIN = open(self.fileINname)
        for line in fileIN:
            if line.find("Begin MODEL") != -1:
                self.Model = next(fileIN).split()[0]
            if line.find("gubmu") != -1:
                self.gUbmu = float(line.split()[1])
            if line.find("gube") != -1:
                self.gUbe = float(line.split()[1])
            if line.find("mub") != -1:
                self.MU = float(line.split()[1])
            if line.find("DECAY") != -1:
                if line.split()[1] == "99":
                    self.WU = float(line.split()[2])
            if line.find("Integrated weight") != -1:
                self.xsec = float(line.split()[5])

        
    def readEvents(self):
        fileIN = open(self.fileINname)
        newEVENT = False
        oneEvent = []
        for line in fileIN:
            if newEVENT: oneEvent.append(line)
            if line.find("<mgrwt>") != -1:
                # the event block ends
                newEVENT = False
            if line.find("</event>") != -1:
                newEVENT = False
                self.eventList.append(oneEvent)
                oneEvent = []
                if len(self.eventList) >= self.Max and self.Max>0: break
            if line.find("<event>") != -1:
                # the event block starts
                newEVENT = True
                oneEvent.append(line)
        fileIN.close()
        return self.eventList

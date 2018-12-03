class uniProtein:
    def __init__(self, id,name,description,organism,fullName,seq):
        self.id=id
        self.name=name
        self.description=description.replace("'","")
        self.organism=organism
        self.fullName=fullName.replace("'","")
        self.seq=seq
    
    def printUniProtein(self, printSeq=True):
        print("Id:\t%s"%self.id)
        print("Name:\t%s"%self.name)
        print("Full name:\t%s"%self.fullName)
        print("Organism:\t%s"%self.organism)
        print("Description:\t%s"%self.description)
        if printSeq:
            print("Sequence: %s"%self.seq)
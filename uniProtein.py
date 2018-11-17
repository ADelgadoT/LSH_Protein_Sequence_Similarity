class uniProtein:
    def __init__(self, id,name,description,organism,fullName,seq):
        self.id=id
        self.name=name
        self.description=description.replace("'","")
        self.organism=organism
        self.fullName=fullName.replace("'","")
        self.seq=seq
    
    def printUniProtein(self):
        print("Id %s",self.id)
        print("Name %s",self.name)
        print("Description %s",self.description)
        print("Organism %s",self.organism)
        print("FullName %s",self.fullName)
        print("Sequence %s",self.seq)
import sqlite3

class UniprotDB(object):

    def __init__(self, databaseName):
        self.db = sqlite3.connect(databaseName)
        self.c = self.db.cursor()

    def close(self):
        self.c.execute("DROP TABLE proteins")
        self.db.close()
        
    def createTables(self):
        self.c.execute("CREATE TABLE IF NOT EXISTS proteins(proteinID TEXT PRIMARY KEY, name TEXT, fullName TEXT, description TEXT, sequence TEXT, organism TEXT)")
        self.db.commit()
    
    def addProtein(self, proteinID, name, fullName, description, sequence, organism):
        self.c.execute("INSERT INTO proteins(proteinID, name, fullName, description, sequence, organism) VALUES ('"+ proteinID + "','" + name +"','"+ fullName +"','" + description +"','" + sequence +"','" + organism +"')")
        self.db.commit()
    
    def extractProteinSeqFromSpecieCsv(self, organism):
        output = open('Protein_Specie.csv', 'w+')
        self.c.execute("SELECT proteinID, sequence FROM proteins WHERE organism ='"+ organism +"'")
        while True: 
            row = self.c.fetchone()
            if row == None: 
                break
            print(row[0] + ';' + row[1], file = output)
   
    def extractProteinSeqFromSpecie(self, organism):
        self.c.execute("SELECT proteinID, sequence FROM proteins WHERE organism ='"+ organism +"'")
        proteins = list(self.c.fetchall())
        return proteins

    def extractProteinSeqFromSpecie25(self, organism):
        self.c.execute("SELECT proteinID, sequence FROM proteins WHERE organism ='"+ organism +"' LIMIT 2500")
        proteins = list(self.c.fetchall())
        return proteins
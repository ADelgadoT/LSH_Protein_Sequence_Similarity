import sqlite3

class ResultsDB(object):

	def __init__(self, databaseName):
		self.db = sqlite3.connect(databaseName)
		self.c = self.db.cursor()

	def close(self):
		self.db.close()

	def deleteBLASTresults(self):
		self.c.execute("DROP TABLE blastresults")
		self.db.commit()
	
	def deleteLSHresults(self):
		self.c.execute("DROP TABLE lshresults")
		self.db.commit()
		
	def createTables(self):
		self.c.execute("CREATE TABLE IF NOT EXISTS blastresults(queryID TEXT, matchID TEXT, identity REAL, alignmentlength INTEGER, PRIMARY KEY (queryID, matchID))")
		self.c.execute("CREATE TABLE IF NOT EXISTS lshresults(queryID TEXT, matchID TEXT, jaccard REAL, PRIMARY KEY (queryID, matchID))")
		self.db.commit()
		
	def addBLASTresult(self, queryID, matchID, identity, alignmentlength):
		try:
			self.c.execute("INSERT INTO blastresults(queryID, matchID, identity, alignmentlength) VALUES ('"+ queryID + "','" + matchID +"','"+ str(identity) +"','" + str(alignmentlength) +"')")
			self.db.commit()
		except sqlite3.IntegrityError as err:
			print('Multiple matches found for: %s - %s. Only the first match was kept.' % (queryID, matchID))
		
	def addLSHresult(self, queryID, matchID, jaccard):
		self.c.execute("INSERT INTO lshresults(queryID, matchID, jaccard) VALUES ('"+ queryID + "','" + matchID +"','"+ str(jaccard) +"')")
		self.db.commit()
		
	def extractBLASTresults(self):
		self.c.execute("SELECT queryID, matchID, identity, alignmentlength FROM blastresults")
		results = list(self.c.fetchall())
		return results
		
	def extractLSHresults(self):
		self.c.execute("SELECT * FROM lshresults")
		results = list(self.c.fetchall())
		return results
		
	def extractIntersect(self):
		self.c.execute("""SELECT lshresults.queryID, lshresults.matchID, blastresults.identity, 
						blastresults.alignmentlength, lshresults.jaccard FROM lshresults 
						INNER JOIN blastresults ON (lshresults.queryID = blastresults.queryID
						AND lshresults.matchID = blastresults.matchID)
						""")
		results = list(self.c.fetchall())
		return results
"""	
	def extractBLASTresultsFromFile(self, organism, filename):
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
		self.c.execute("SELECT proteinID, sequence FROM proteins WHERE organism ='"+ organism +"'")
		proteins = list(self.c.fetchall())
		return proteins


"""
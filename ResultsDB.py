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
		
	def deleteTable(self, tablename):
		self.c.execute("DROP TABLE " + tablename)
		self.db.commit()
		
	def createBLASTtable(self):
		self.c.execute("CREATE TABLE IF NOT EXISTS blastresults(queryID TEXT, matchID TEXT, identity REAL, alignmentlength INTEGER, PRIMARY KEY (queryID, matchID))")
		self.db.commit()
		
	def createLSHtable(self, tablename):
		self.c.execute("CREATE TABLE IF NOT EXISTS "+ tablename +"(queryID TEXT, matchID TEXT, jaccard REAL, PRIMARY KEY (queryID, matchID))")
		self.db.commit()
		
	def addBLASTresult(self, queryID, matchID, identity, alignmentlength):
		self.c.execute("INSERT INTO blastresults(queryID, matchID, identity, alignmentlength) VALUES ('"+ queryID + "','" + matchID +"','"+ str(identity) +"','" + str(alignmentlength) +"')")
		self.db.commit()
		
	def addLSHresult(self, queryID, matchID, jaccard, tablename="lshresults"):
		self.c.execute("INSERT INTO "+ tablename +"(queryID, matchID, jaccard) VALUES ('"+ queryID + "','" + matchID +"','"+ str(jaccard) +"')")
		self.db.commit()
		
	def extractBLASTresults(self, identity=0.0, alignmentlength=0):
		self.c.execute("SELECT queryID, matchID, identity, alignmentlength FROM blastresults WHERE (identity >= "+ identity +"AND alignmentlength >="+ alignmentlength +")")
		results = list(self.c.fetchall())
		return results
		
	def extractBLASTresults(self, identity=0.0, alignmentlength=0):
		self.c.execute("SELECT queryID, matchID, identity, alignmentlength FROM blastresults WHERE (identity >= "+ identity +"AND alignmentlength >="+ alignmentlength +")")
		results = list(self.c.fetchall())
		return results
		
	def extractLSHresults(self, tablename="lshresults", jaccard=0.0):
		self.c.execute("SELECT * FROM " + tablename + "WHERE jaccard >=" + jaccard)
		results = list(self.c.fetchall())
		return results
		
	def extractLSHcount(self, tablename="lshresults", jaccard=0.0):
		self.c.execute("SELECT COUNT(*) FROM " + tablename + "WHERE jaccard >=" + jaccard)
		results = list(self.c.fetchall())
		return results
		
	def extractCount(self, tablename="lshresults"):
		self.c.execute("SELECT COUNT(*) FROM " + tablename)
		results = list(self.c.fetchall())
		return results
		
	def extractIntersect(self, lshtablename = "lshresults", identity = 0.0, alignmentlength = 0, jaccard = 0):
		self.c.execute("""SELECT lshreslts.queryID, lshreslts.matchID, blastresults.identity, 
						blastresults.alignmentlength, lshreslts.jaccard FROM blastresults INNER JOIN""" \
						+ lshtablename+ \
						""" lshreslts ON (lshreslts.queryID = blastresults.queryID
						AND lshreslts.matchID = blastresults.matchID)
						WHERE blastresults.identity >= """ \
						+ identity + \
						"blastresults.alignmentlength >=" \
						+ alignmentlength + \
						"lshreslts.jaccard >=" \
						+ jaccard
						)
		results = list(self.c.fetchall())
		return results
		
		
	def extractIntersect(self, lshtablename = "lshresults", identity = 0.0, alignmentlength = 0, jaccard = 0):
		self.c.execute("""SELECT COUNT(*) FROM blastresults INNER JOIN""" \
						+ lshtablename+ \
						""" lshreslts ON (lshreslts.queryID = blastresults.queryID
						AND lshreslts.matchID = blastresults.matchID)
						WHERE blastresults.identity >= """ \
						+ identity + \
						"blastresults.alignmentlength >=" \
						+ alignmentlength+ \
						"lshreslts.jaccard >=" \
						+ jaccard
						)
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
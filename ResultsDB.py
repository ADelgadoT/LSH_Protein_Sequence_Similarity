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
		try:
			self.c.execute("INSERT INTO blastresults(queryID, matchID, identity, alignmentlength) VALUES ('"+ queryID + "','" + matchID +"','"+ str(identity) +"','" + str(alignmentlength) +"')")
			self.db.commit()
		except sqlite3.IntegrityError as err:
			print('Multiple matches found for: %s - %s. Only the first match was kept.' % (queryID, matchID))
		
	def addLSHresult(self, queryID, matchID, jaccard, tablename="lshresults"):
		self.c.execute("INSERT INTO "+ tablename +"(queryID, matchID, jaccard) VALUES ('"+ queryID + "','" + matchID +"','"+ str(jaccard) +"')")
		self.db.commit()
		
	def extractBLASTresults(self, identity=0.0, alignmentlength=0):
		self.c.execute("SELECT queryID, matchID, identity, alignmentlength FROM blastresults WHERE (identity >= "+ str(identity) +" AND alignmentlength >= "+ str(alignmentlength) +")")
		results = list(self.c.fetchall())
		return results
		
	def extractBLASTalignmentlengths(self, identity=0.0, alignmentlength=0):
		self.c.execute("SELECT alignmentlength FROM blastresults WHERE (identity >= "+ str(identity) +" AND alignmentlength >= "+ str(alignmentlength) +")")
		results = [l[0] for l in self.c.fetchall()]
		return results
		
	def extractBLASTidentities(self, identity=0.0, alignmentlength=0):
		self.c.execute("SELECT identity FROM blastresults WHERE (identity >= "+ str(identity) +" AND alignmentlength >= "+ str(alignmentlength) +")")
		results = [id[0] for id in self.c.fetchall()]
		return results
		
	def extractBLASTcount(self, identity=0.0, alignmentlength=0):
		self.c.execute("SELECT COUNT(*) FROM blastresults WHERE (identity >= "+ str(identity) +" AND alignmentlength >= "+ str(alignmentlength) +")")
		results = self.c.fetchone()[0]
		return results
		
	def extractBLASTcountPerProtein(self, queryID, identity=0.0, alignmentlength=0):
		self.c.execute("SELECT COUNT(*) FROM blastresults WHERE (identity >= "+ str(identity) +" AND alignmentlength >= "+ str(alignmentlength) + " AND queryID = '" +queryID+"')")
		results = self.c.fetchone()[0]
		return results
		
		
	def extractLSHresults(self, tablename="lshresults", jaccard=0.0):
		self.c.execute("SELECT * FROM " + tablename + " WHERE jaccard >= " + str(jaccard))
		results = list(self.c.fetchall())
		return results
		
	def extractLSHcount(self, tablename="lshresults", jaccard=0.0):
		self.c.execute("SELECT COUNT(*) FROM " + tablename + " WHERE jaccard >= " + str(jaccard))
		results = self.c.fetchone()[0]
		return results
		
	def extractLSHcountPerProtein(self, queryID, tablename="lshresults", jaccard=0.0):
		self.c.execute("SELECT COUNT(*) FROM " + tablename + " WHERE jaccard >= " + str(jaccard) + " AND queryID = '" + queryID +"'")
		results = self.c.fetchone()[0]
		return results
		
	def extractCount(self, tablename="lshresults"):
		self.c.execute("SELECT COUNT(*) FROM " + tablename)
		results = self.c.fetchone()[0]
		return results
		

	
	def extractIntersect(self, lshtablename = "lshresults", identity = 0.0, alignmentlength = 0, jaccard = 0):
		self.c.execute("SELECT lshreslts.queryID, lshreslts.matchID, blastresults.identity, \
						blastresults.alignmentlength, lshreslts.jaccard \
						FROM blastresults INNER JOIN "+lshtablename+\
						" lshreslts ON (lshreslts.queryID = blastresults.queryID AND lshreslts.matchID = blastresults.matchID)" +\
						" WHERE (blastresults.identity >= "+ str(identity) + \
						" AND blastresults.alignmentlength >= " + str(alignmentlength)+ \
						" AND lshreslts.jaccard >= " + str(jaccard) +")")
					
		results = self.c.fetchone()[0]
		return results
	
	def extractIntersectCount(self, lshtablename = "lshresults", identity = 0.0, alignmentlength = 0, jaccard = 0):
		self.c.execute("SELECT COUNT(*) FROM blastresults INNER JOIN "+lshtablename+\
						" lshreslts ON (lshreslts.queryID = blastresults.queryID AND lshreslts.matchID = blastresults.matchID)" +\
						" WHERE (blastresults.identity >= "+ str(identity) + \
						" AND blastresults.alignmentlength >= " + str(alignmentlength)+ \
						" AND lshreslts.jaccard >= " + str(jaccard) +")")
					
		results = self.c.fetchone()[0]
		return results

	def extractIntersectCountPerProtein(self, queryID, lshtablename = "lshresults", identity = 0.0, alignmentlength = 0, jaccard = 0):
		self.c.execute("SELECT COUNT(*) FROM blastresults INNER JOIN "+lshtablename+\
						" lshreslts ON (lshreslts.queryID = blastresults.queryID AND lshreslts.matchID = blastresults.matchID)" +\
						" WHERE (blastresults.identity >= "+ str(identity) + \
						" AND blastresults.alignmentlength >= " + str(alignmentlength)+ \
						" AND lshreslts.jaccard >= " + str(jaccard) +\
						" AND blastresults.queryID = '" + queryID +\
						"')")
					
		results = self.c.fetchone()[0]
		return results		
	
		
	def extractIntersectIdentity(self, lshtablename = "lshresults", identity = 0.0, alignmentlength = 0, jaccard = 0):
		self.c.execute("SELECT blastresults.identity FROM blastresults INNER JOIN "+lshtablename+\
						" lshreslts ON (lshreslts.queryID = blastresults.queryID AND lshreslts.matchID = blastresults.matchID)" +\
						" WHERE (blastresults.identity >= "+ str(identity) + \
						" AND blastresults.alignmentlength >= " + str(alignmentlength)+ \
						" AND lshreslts.jaccard >= " + str(jaccard) +")")
					
		results = [id[0] for id in self.c.fetchall()]
		return results


	def extractIntersectAlignmentLength(self, lshtablename = "lshresults", identity = 0.0, alignmentlength = 0, jaccard = 0):
		self.c.execute("SELECT blastresults.alignmentlength FROM blastresults INNER JOIN "+lshtablename+\
						" lshreslts ON (lshreslts.queryID = blastresults.queryID AND lshreslts.matchID = blastresults.matchID)" +\
						" WHERE (blastresults.identity >= "+ str(identity) + \
						" AND blastresults.alignmentlength >= " + str(alignmentlength)+ \
						" AND lshreslts.jaccard >= " + str(jaccard) +")")
					
		results = [l[0] for l in self.c.fetchall()]
		return results

	def extractIntersectJaccard(self, lshtablename = "lshresults", identity = 0.0, alignmentlength = 0, jaccard = 0):
		self.c.execute("SELECT lshreslts.jaccard FROM blastresults INNER JOIN "+lshtablename+\
						" lshreslts ON (lshreslts.queryID = blastresults.queryID AND lshreslts.matchID = blastresults.matchID)" +\
						" WHERE (blastresults.identity >= "+ str(identity) + \
						" AND blastresults.alignmentlength >= " + str(alignmentlength)+ \
						" AND lshreslts.jaccard >= " + str(jaccard) +")")
					
		results = [id[0] for id in self.c.fetchall()]
		return results
		
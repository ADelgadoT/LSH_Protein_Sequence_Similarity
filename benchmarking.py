import sys
from ProteinsManager import ProteinsManager
from UniprotDB import UniprotDB
from uniProtein import uniProtein
from lsh import LSH
from ResultsDB import ResultsDB
from collections import OrderedDict
import matplotlib.pyplot as plt
import pandas as pd

uniDB = UniprotDB("Uniprot_DB.sqlite")
minhash3 = LSH(0.2,128)
minhash4 = LSH(0.2,128)
minhash5 = LSH(0.2,128)

#DB creation:
"""
protManager = ProteinsManager()
uniDB.createTables()
protManager.loadProteins("Ecolx.xml",uniDB)
protManager.loadProteins("PseA7.xml",uniDB)
"""

# Create the minhashes
proteins = uniDB.extractProteins()
minhashes3, lsh3 = minhash3.calculateLSH(proteins, 3)
minhashes4, lsh4 = minhash4.calculateLSH(proteins, 4)
minhashes5, lsh5 = minhash5.calculateLSH(proteins, 5)

resultsDB = ResultsDB("Results_DB.sqlite")

"""
#Create the BLAST table in the results database
handle = open('all_results_nofilter.txt', 'r')
resultsDB = ResultsDB("Results_DB.sqlite")
resultsDB.createBLASTtable()
resultsDB.deleteBLASTresults()
resultsDB.createBLASTtable()
for line in handle:
	line = line[:-1].split('\t')
	# Extract accessions from 'sp|A0A0R6L508|MCR1_ECOLX'-like string
	line[0] = line[0].split('|')[1]
	line[1] = line[1].split('|')[1]
	print(line)
	# Filter self-matches, add to the database
	if line[0] != line[1]:
		resultsDB.addBLASTresult(line[0], line[1], line[2], line[3])
"""

# Create the LSH tables in the results database
resultsDB.createLSHtable("lshresults3")
resultsDB.deleteTable("lshresults3")
resultsDB.createLSHtable("lshresults3")
for query in minhash3.minhashes.keys():
	matches = minhash3.queryProtein(query)
	for match in matches:
		# Filter self-matches
		if query != match:
			jaccard = minhash3.estimateJaccard(query, match)
			resultsDB.addLSHresult(query, match, jaccard, "lshresults3")
			
# Create the LSH tables in the results database
resultsDB.createLSHtable("lshresults4")
resultsDB.deleteTable("lshresults4")
resultsDB.createLSHtable("lshresults4")
for query in minhash4.minhashes.keys():
	matches = minhash4.queryProtein(query)
	for match in matches:
		# Filter self-matches
		if query != match:
			jaccard = minhash4.estimateJaccard(query, match)
			resultsDB.addLSHresult(query, match, jaccard, "lshresults4")

# Create the LSH tables in the results database
resultsDB.createLSHtable("lshresults5")
resultsDB.deleteTable("lshresults5")
resultsDB.createLSHtable("lshresults5")
for query in minhash5.minhashes.keys():
	matches = minhash5.queryProtein(query)
	for match in matches:
		# Filter self-matches
		if query != match:
			jaccard = minhash5.estimateJaccard(query, match)
			resultsDB.addLSHresult(query, match, jaccard, "lshresults5")			

print("BLAST results with 80% identity and length > 100, LSH results with Jaccard > 0.5")
print("BLAST: ", resultsDB.extractCount('blastresults', 80.0, 100))
print("LSH3: ", resultsDB.extractLSHcount('lshresults3', 0.5))
print("Intersect: ", resultsDB.extractIntersectCount('lshresults3', 80.0, 100, 0.5))
print("LSH4: ", resultsDB.extractLSHcount('lshresults4'), 0.5)
print("Intersect: ", resultsDB.extractIntersectCount('lshresults4', 80.0, 100, 0.5))
print("LSH5: ", resultsDB.extractLSHcount('lshresults5'), 0.5)
print("Intersect: ", resultsDB.extractIntersectCount('lshresults5', 80.0, 100, 0.5))


"""

#tp = intersect
#fp = lshResults - intersect
#fn = blastResults
#precision = tp/(tp+fp)
#recall = tp/(tp+fn)
#print("Precision:", precision)
#print("Recall:", recall)

fig1, ax1 = plt.subplots()
ax1.set_title('Identity distribution')
ax1.boxplot([intersect.iloc[:,2], blastResults.iloc[:,2]], labels=["LSH", "BLAST"])
ax1.set_ylabel('Alignment Identity (%)')
plt.savefig('test.png')
"""
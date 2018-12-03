import sys
from ProteinsManager import ProteinsManager
from UniprotDB import UniprotDB
from uniProtein import uniProtein
from lsh import LSH
from ResultsDB import ResultsDB
from collections import OrderedDict
import matplotlib.pyplot as plt
import pandas as pd
import sqlite3


uniDB = UniprotDB("Uniprot_DB.sqlite")
#Construct the protein database
"""
uniDB.deleteProteins()
protManager = ProteinsManager()
uniDB.createTables()
protManager.loadProteins("Ecolx.xml",uniDB)
protManager.loadProteins("PseA7.xml",uniDB)
"""


minhash3 = LSH(0.3,32)
minhash4 = LSH(0.3,64)
minhash4b = LSH(0.3,96)
minhash5 = LSH(0.3,128)

# Create the minhashes
proteins = uniDB.extractProteins()
"""
minhashes3, lsh3 = minhash3.calculateLSH(proteins, 3)
minhashes4, lsh4 = minhash4.calculateLSH(proteins, 3)
minhashes4b, lsh4b = minhash5.calculateLSH(proteins, 3)
minhashes5, lsh5 = minhash6.calculateLSH(proteins, 3)

minhash3.saveLSH(32)
minhash4.saveLSH(64)
minhash4b.saveLSH(96)
minhash5.saveLSH(128)
"""

minhash3.loadLSH(32)
minhash4.loadLSH(64)
minhash4b.loadLSH(96)
minhash5.loadLSH(128)


resultsDB = ResultsDB("Results_DB_permutations.sqlite")

"""
#Create the BLAST table in the results database
handle = open('all_results_nofilter.txt', 'r')
resultsDB = ResultsDB("Results_DB_permutations.sqlite")
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

# Create the LSH tables in the results database
resultsDB.createLSHtable("lshresults3")
resultsDB.deleteTable("lshresults3")
resultsDB.createLSHtable("lshresults3")
for query in minhash3.minhashes.keys():
	print(query)
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
resultsDB.createLSHtable("lshresults4b")
resultsDB.deleteTable("lshresults4b")
resultsDB.createLSHtable("lshresults4b")
for query in minhash4b.minhashes.keys():
	matches = minhash4b.queryProtein(query)
	for match in matches:
		# Filter self-matches
		if query != match:
			jaccard = minhash4b.estimateJaccard(query, match)
			resultsDB.addLSHresult(query, match, jaccard, "lshresults4b")			
			
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
"""
			




"""
for tablename in ['lshresults3','lshresults4','lshresults5']:
	intersect = resultsDB.extractIntersectCount(tablename, 0, 30, 0.3)
	lshresults = resultsDB.extractLSHcount(tablename, 0.3)
	blastresults = resultsDB.extractBLASTcount(0, 30)
	tp = intersect
	fp = lshresults - intersect
	fn = blastresults
	tn = (len(proteins)-1) ** 2 - tp - fp - fn
	precision = tp/(tp+fp)
	recall = tp/(tp+fn)
	accuracy = (tp + tn)/(tp + tn + fp + fn)
	print("%s \n blastresults: %i\n lshresults: %i \n intersection: %i \nprecision: %0.3f recall: %0.3f accuracy: %0.3f\n" \
		% (tablename, blastresults, lshresults, intersect, precision, recall, accuracy))

		
for x in resultsDB.extractBLASTresults(80.0, 100):		
	print(x)
	
print("XXXXXXXXXXXXXXXXXXXXXXX LSH RESULTS XXXXXXXXXXXXX")
for x in resultsDB.extractLSHresults('lshresults3', 0.5):
	print(x)
"""

all_precisions = []
all_recalls = []
all_identities = []
all_alignments = []
identity_th, alignment_th, jaccard_th = 80.0, 100, 0.5 
for tablename in ['lshresults4b']:
	precisions = []
	recalls = []
	
	overall_intersect = resultsDB.extractIntersectCount(tablename, identity_th, alignment_th, jaccard_th)
	nothreshold_intersect = resultsDB.extractIntersectCount(tablename, 0, 0, jaccard_th)
	overall_lshresults = resultsDB.extractLSHcount(tablename, jaccard_th)
	overall_blastresults = resultsDB.extractBLASTcount(identity_th, alignment_th)
	#print("%s \n blastresults: %i\n lshresults: %i \n intersection: %i \n nr of LSH results matched with a BLAST identity and alignmentlength: %s \n nr of LSH results not matched with BLAST: %s \n" \
	#	% (tablename, overall_blastresults, overall_lshresults, overall_intersect, nothreshold_intersect, overall_lshresults - nothreshold_intersect))
	
	all_identities.append(resultsDB.extractIntersectIdentity(tablename, 0, 0, jaccard_th)) 
	all_alignments.append(resultsDB.extractIntersectAlignmentLength(tablename, 0, 0, jaccard_th))
	all_jaccards = resultsDB.extractIntersectJaccard(tablename, 0, 0, jaccard_th)
	
	for query in proteins:
		intersect = resultsDB.extractIntersectCountPerProtein(query[0], tablename, identity_th, alignment_th, jaccard_th)
		lshresults = resultsDB.extractLSHcountPerProtein(query[0],tablename, jaccard_th)
		blastresults = resultsDB.extractBLASTcountPerProtein(query[0], identity_th, alignment_th)
		tp = intersect
		fp = lshresults - intersect
		fn = blastresults - intersect
		#tn = (len(proteins)-1) ** 2 - tp - fp - fn
		precision = tp/(tp+fp) if (tp+fp) != 0 else -1
		recall = tp/(tp+fn) if (tp+fn) != 0 else -1
		#accuracy = (tp + tn)/(tp + tn + fp + fn)
					
		if lshresults == 0 and blastresults == 0:
			precision = -1
			recall = -1
		if precision != -1:
			precisions.append(precision)
		if recall != -1:
			recalls.append(recall)
		#print(query[0], intersect, lshresults, blastresults, precision, recall)
	
	all_precisions.append(precisions)
	all_recalls.append(recalls)
	print("%s \n blastresults: %i\n lshresults: %i \n intersection: %i \n nr of LSH results matched with a BLAST identity and alignmentlength: %s \n nr of LSH results not matched with BLAST: %s \n" \
		% (tablename, overall_blastresults, overall_lshresults, overall_intersect, nothreshold_intersect, overall_lshresults - nothreshold_intersect))
	
		
for i in range(len(all_precisions)):
	plt.hist(all_precisions[i], bins=100)
	plt.savefig('precisions%i.png' % i)
	plt.close()
	plt.hist(all_recalls[i], bins=100)
	plt.savefig('recalls%i.png' % i)
	plt.close()
	print("Precision LSH%i: "% i, sum(all_precisions[i])/len(all_precisions[i]))
	print("Recall LSH%i: " % i,sum(all_recalls[i])/len(all_recalls[i]))
	print("Lengths precisions, recalls: %i, %i" %(len(all_precisions[i]), len(all_recalls[i])))

# Add the blast results to the identities and alignments
all_identities.append(resultsDB.extractBLASTidentities(identity_th, alignment_th))
all_alignments.append(resultsDB.extractBLASTalignmentlengths(identity_th, alignment_th))
	
print(len(proteins), "proteins in total")

import numpy as np
from numpy.polynomial.polynomial import polyfit
b, m = polyfit(all_jaccards,all_identities[0], 1)




fig1, ax1 = plt.subplots()
#ax1.set_title('Jaccard ~ identity')
ax1.scatter(all_jaccards,all_identities[0], s=[x/100 for x in all_alignments[0]], c='b')
#ax1.plot(all_jaccards, b + m * np.array(all_jaccards), '-')
ax1.set_ylabel('Alignment Identity (%)')
ax1.set_xlabel('Jaccard Similarity Score')


legend_sizes = np.array([47, 393, 1377])
# get the indices for each of the legend sizes
indices = [np.where(all_alignments[0]==v)[0][0] for v in legend_sizes]
print(indices)
# plot each point again, and its value as a label
legendlabels = [50, 400, 1400]
for n, i in enumerate(indices):
    ax1.scatter(all_jaccards[i],all_identities[0][i],s=[x/100 for x in [all_alignments[0][i]]],
                  label='{:.0f}'.format(legendlabels[n]), c='b')
# add the legend
ax1.legend(scatterpoints=1, loc='lower right', title='Alignment Length')
ax1.scatter(all_jaccards,all_identities[0], s=[x/100 for x in all_alignments[0]])

plt.savefig('jaccard_identity_alignment.png')
plt.close()

fig1, ax1 = plt.subplots()
#ax1.set_title('Jaccard ~ alignmentlength')
ax1.scatter(all_jaccards,all_alignments[0], s=3)
ax1.set_ylabel('Alignment Length (bp)')
ax1.set_xlabel('Jaccard Similarity Score')
plt.savefig('jaccard_length.png')
plt.close()


"""	
fig1, ax1 = plt.subplots()
ax1.set_title('Identity distribution')
ax1.boxplot(all_identities, labels=["LSH 32 perm.", "LSH 64 perm.","LSH 96 perm.","LSH 128 perm.","BLAST"])
ax1.set_ylabel('Alignment Identity (%)')
plt.savefig('identity_perm.png')
plt.close()

fig1, ax1 = plt.subplots()
ax1.set_title('Alignment length distribution')
ax1.boxplot(all_alignments, labels=["LSH 32 perm.", "LSH 64 perm.","LSH 96 perm.","LSH 128 perm.","BLAST"])
ax1.set_ylabel('Alignment Length (basepairs)')
plt.savefig('alignmentlength_perm.png')
plt.close()
"""
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
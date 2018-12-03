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


uniDB = UniprotDB("Uniprot_DB_ec_pa_human.sqlite")
#Construct the protein database
"""
uniDB.deleteProteins()
protManager = ProteinsManager()
uniDB.createTables()
protManager.loadProteins("Ecolx.xml",uniDB)
protManager.loadProteins("PseA7.xml",uniDB)
"""
protManager = ProteinsManager()
protManager.loadProteins("Human.xml",uniDB)

minhash3 = LSH(0.3,96)
minhash4 = LSH(0.5,96)
minhash5 = LSH(0.5,128)


# Create the minhashes
proteins = uniDB.extractProteins()

minhashes3, lsh3 = minhash3.calculateLSH(proteins, 3)
minhashes4, lsh4 = minhash4.calculateLSH(proteins, 3)
minhashes5, lsh5 = minhash5.calculateLSH(proteins, 3)

minhash3.saveLSH(963)
minhash4.saveLSH(965)
minhash5.saveLSH(1285)
"""
minhash3.loadLSH(963)
minhash4.loadLSH(1283)
minhash5.loadLSH(1285)
"""

resultsDB = ResultsDB("Results_DB_ec_pa_human.sqlite")

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

"""
"""
all_precisions = []
all_recalls = []
all_identities = []
all_alignments = []
identity_th, alignment_th, jaccard_th = 80.0, 100, 0.5 
for tablename in ['lshresults3']:
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
	
identity_th, alignment_th, jaccard_th = 80.0, 100, 0.4 
for tablename in ['lshresults4']:
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
	

identity_th, alignment_th, jaccard_th = 80.0, 100, 0.3 
for tablename in ['lshresults5']:
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
	
fig1, ax1 = plt.subplots()
ax1.set_title('Identity distribution')
ax1.boxplot(all_identities, labels=["LSH 3-shingling", "LSH 4-shingling","LSH 5-shingling","BLAST"])
ax1.set_ylabel('Alignment Identity (%)')
plt.savefig('identity.png')
plt.close()

fig1, ax1 = plt.subplots()
ax1.set_title('Alignment length distribution')
ax1.boxplot(all_alignments, labels=["LSH 3-shingling", "LSH 4-shingling","LSH 5-shingling","BLAST"])
ax1.set_ylabel('Alignment Length (basepairs)')
plt.savefig('alignmentlength.png')
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
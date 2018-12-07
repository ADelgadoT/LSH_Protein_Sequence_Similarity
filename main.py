import sys
from ProteinsManager import ProteinsManager
from UniprotDB import UniprotDB
from uniProtein import uniProtein
from lsh import LSH
from ResultsDB import ResultsDB
from collections import OrderedDict
import matplotlib.pyplot as plt
import pandas as pd
import time

class Analyzer(object):

	def run(self):
		print(\
		"""Local Sensitivity Hashing-based protein similarity search.
	Options: E[X]it, [L]oad Database, [D]elete Database,
	[C]alculate LSH, [RC] Recalculate LSH, [LL] Load LSH, [S]ave LSH
	[Q]uery LSH, Query [A]ll LSH, Read [B]LAST, Compare [R]esults,
		""")
		mode=input('Choose option:')
		
		uniDB = UniprotDB("Uniprot_DB.sqlite")
		minhash = LSH(0.5,96)
		
		while(mode!='Exit' and mode!='X'):

			if (mode=='Delete Database' or mode=='D'):
				uniDB.deleteProteins()

			if (mode=='Load Database' or mode=='L'):
				protManager = ProteinsManager()
				uniDB.createTables()
				filename = input('XML filename (e.g. Ecolx.xml or PseA7.xml or Human.xml): ')
				protManager.loadProteins(filename,uniDB)

			if (mode=='Calculate LSH' or mode=='C'):
				uniDB = UniprotDB("Uniprot_DB.sqlite")
				proteins = uniDB.extractProteins()
				minhashes, lsh = minhash.calculateLSH(proteins, 3)
				print("Calculated")

			if (mode=='Recalculate LSH' or mode=='RC'):
				jaccardThreshold = float(input("Specify a Jaccard similarity threshold (default: 0.5): "))
				permutations = int(input("Specify the number of permutations(default: 96) : "))
				shinglesize = int(input("Specify the shingle size (default: 3): "))
				minhash = LSH(jaccardThreshold, permutations)
				proteins = uniDB.extractProteins()
				minhashes, lsh = minhash.calculateLSH(proteins, shinglesize)
				print("Recalculated")

			if (mode=='Query LSH' or mode=='Q'):
				protein = input('Protein accession: ')
				start_time = time.time()
				result = minhash.queryProtein(protein)
				if result is not None:
					jaccResultsDict = minhash.checkJaccardResultsOfProtein(protein, result)
					# Return the results in sorted order, big to small Jaccard score
					sorted_jaccResultsDict = OrderedDict(sorted(jaccResultsDict.items(), key=lambda x: -x[1]))
					for jaccRes in sorted_jaccResultsDict.items():
						print("\nMatch with Jaccard:",jaccRes[1])
						information = uniDB.extractProteinInformation(jaccRes[0])
						proteininfo = uniProtein(*information)
						proteininfo.printUniProtein(printSeq=False)
				print("Runtime of query search: %s seconds " % (time.time() - start_time))

			if (mode=='Calculate All' or mode=='CA'):
				start_time = time.time()
				uniDB = UniprotDB("Uniprot_DB.sqlite")
				#uni_DB.close()
				proteins = uniDB.extractProteins()
				#minhash.calculateLSH([protein[1] for protein in proteins])
				minhashes, lsh = minhash.calculateLSH(proteins, 3)
				for protein in proteins:
					print("Protein ", protein[0])
					result = minhash.queryProtein(protein[0])
					if result is not None:
						jaccResultsDict = minhash.checkJaccardResultsOfProtein(protein[0], result)
						sorted_jaccResultsDict = OrderedDict(sorted(jaccResultsDict.items(), key=lambda x: -x[1]))
						for jaccRes in sorted_jaccResultsDict.items():
							print(jaccRes[0]," - Jaccard: ",jaccRes[1])
				print("Runtime of query all: %s seconds " % (time.time() - start_time))
					
			if (mode=='Query All LSH' or mode=='A'):
				resultsDB = ResultsDB("Results_DB.sqlite")
				resultsDB.createLSHtable("lshresults")
				resultsDB.deleteTable("lshresults")
				resultsDB.createLSHtable("lshresults")
				for query in minhash.minhashes.keys():
					matches = minhash.queryProtein(query)
					for match in matches:
						# Filter self-matches
						if query != match:
							jaccard = minhash.estimateJaccard(query, match)
							resultsDB.addLSHresult(query, match, jaccard, "lshresults")
				print(resultsDB.extractLSHresults("lshresults"))
					
			if (mode=='Read BLAST Results'  or mode=='B'):
				filename = input('Filename: ')
				handle = open(filename, 'r')
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
				print(resultsDB.extractBLASTresults())

			if (mode=='Compare Results' or mode=='R'):

				# Database with all LSH and BLASTp results
                resultsDB = ResultsDB("Results_DB.sqlite")
				identity_th, alignment_th, jaccard_th = 80.0, 100, 0.5 
				precisions = []
				recalls = []
				
                # Load in all protein ids to loop over
				uniDB = UniprotDB("Uniprot_DB.sqlite")
				proteins = uniDB.extractProteins()
				# Store all precisions and recalls per query, to calculate the average
				for query in proteins:
					intersect = resultsDB.extractIntersectCountPerProtein(query[0], 'lshresults', identity_th, alignment_th, jaccard_th)
					lshresults = resultsDB.extractLSHcountPerProtein(query[0],'lshresults', jaccard_th)
					blastresults = resultsDB.extractBLASTcountPerProtein(query[0], identity_th, alignment_th)
					tp = intersect
					fp = lshresults - intersect
					fn = blastresults - intersect
					precision = tp/(tp+fp) if (tp+fp) != 0 else -1
					recall = tp/(tp+fn) if (tp+fn) != 0 else -1
					# Exclude results without any similar proteins / division by zero
					if precision != -1:
						precisions.append(precision)
					if recall != -1:
						recalls.append(recall)

				print("Comparison of BLAST and LSH results:\n Number of proteins queried: %i \n Average precision: %0.3f Average recall: %0.3f\n" \
					% (len(proteins), sum(precisions)/len(precisions), sum(recalls)/len(recalls)))

					
				
			if (mode=='Save LSH' or mode=='S'):
				number = int(input('Suffix number: '))
				minhash.saveLSH(number)

			if (mode=='Load LSH' or mode=='LL'):
				number = int(input('Suffix number: '))
				minhash.loadLSH(number)

			mode = input('Choose option: ')
	   

if __name__ == '__main__':
	Analyzer().run()

import sys
from ProteinsManager import ProteinsManager
from UniprotDB import UniprotDB
from uniProtein import uniProtein
from lsh import LSH
from ResultsDB import ResultsDB
from collections import OrderedDict
import matplotlib.pyplot as plt
import pandas as pd
#from Bio import SeqIO

class Analyzer(object):

	def run(self):

		mode=input('Choose option:')
		
		uniDB = UniprotDB("Uniprot_DB.sqlite")
		minhash = LSH(0.2,128)
		
		while(mode!='Exit'):
			#print(mode)
			if (mode=='Delete Database'):
				uniDB.deleteProteins()

			if (mode=='Load Database'):
				#DB creation:
				protManager = ProteinsManager()
				#uniDB = UniprotDB("Uniprot_DB.sqlite")
				uniDB.createTables()
				#loadProteins("uniprot_thaliana.xml",uni_DB)
				protManager.loadProteins("Ecolx.xml",uniDB)
				protManager.loadProteins("PseA7.xml",uniDB)

			#uniDB = UniprotDB("Uniprot_DB.sqlite")
			#uni_DB.close()
			#uniDB.extractProteinSeqFromSpecieCsv("Arabidopsis thaliana (Mouse-ear cress)")

			if (mode=='Calculate LSH'):
				uniDB = UniprotDB("Uniprot_DB.sqlite")
				#uni_DB.close()
				proteins = uniDB.extractProteins()
				#minhash.calculateLSH([protein[1] for protein in proteins])
				minhashes, lsh = minhash.calculateLSH(proteins, 3)
				print(minhashes.keys())

			if (mode=='Recalculate LSH'):
				jaccardThreshold = float(input("Specify a Jaccard similarity threshold: "))
				minhash = LSH(jaccardThreshold,128)
				uniDB = UniprotDB("Uniprot_DB.sqlite")
				#uni_DB.close()
				proteins = uniDB.extractProteins()
				#minhash.calculateLSH([protein[1] for protein in proteins])
				minhashes, lsh = minhash.calculateLSH(proteins, 3)
				print(minhashes.keys())

			if (mode=='Query'):
				protein = input('Protein accession: ')
				result = minhash.queryProtein(protein)
				jaccResultsDict = minhash.checkJaccardResultsOfProtein(protein, result)
				sorted_jaccResultsDict = OrderedDict(sorted(jaccResultsDict.items(), key=lambda x: -x[1]))
				for jaccRes in sorted_jaccResultsDict.items():
					print(jaccRes[0]," - Jaccard: ",jaccRes[1])
					
			if (mode=='LSH Query All'):
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
					
			if (mode=='Read BLAST Results'):
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

			if (mode=='Compare results'):
				resultsDB = ResultsDB("Results_DB.sqlite")
				print("BLAST: ", len(resultsDB.extractBLASTresults()))
				print("LSH: ", len(resultsDB.extractLSHresults()))
				print("Intersect: ", len(resultsDB.extractIntersect()))
				
				lshResults = pd.DataFrame(resultsDB.extractLSHresults())
				blastResults = pd.DataFrame(resultsDB.extractBLASTresults())
				intersect = pd.DataFrame(resultsDB.extractIntersect())
				
				tp = intersect.shape[0]
				fp = lshResults.shape[0] - intersect.shape[0]
				fn = blastResults.shape[0]
				precision = tp/(tp+fp)
				recall = tp/(tp+fn)
				print("Precision:", precision)
				print("Recall:", recall)
				
				fig1, ax1 = plt.subplots()
				ax1.set_title('Identity distribution')
				ax1.boxplot([intersect.iloc[:,2], blastResults.iloc[:,2]], labels=["LSH", "BLAST"])
				ax1.set_ylabel('Alignment Identity (%)')
				plt.savefig('test.png')
					
			if (mode=='Save LSH'):
				minhash.saveLSH()
			if (mode=='Load LSH'):
				minhash.loadLSH()

			mode = input('Choose option: ')
	   

if __name__ == '__main__':
	Analyzer().run()

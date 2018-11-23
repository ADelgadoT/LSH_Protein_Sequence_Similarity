import sys
from ProteinsManager import ProteinsManager
from UniprotDB import UniprotDB
from uniProtein import uniProtein
from lsh import LSH
from collections import OrderedDict
#from Bio import SeqIO

class Analyzer(object):

	def run(self):

		mode=input('Choose option:')
		
		uniDB = UniprotDB("Uniprot_DB.sqlite")
		minhash = LSH(0.32,128)
		
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
				minhashes, lsh = minhash.calculateLSH(proteins)
				print(minhashes.keys())

			if (mode=='Recalculate LSH'):
				minhash = LSH(0.32,128)
				uniDB = UniprotDB("Uniprot_DB.sqlite")
				#uni_DB.close()
				proteins = uniDB.extractProteins()
				#minhash.calculateLSH([protein[1] for protein in proteins])
				minhashes, lsh = minhash.calculateLSH(proteins)
				print(minhashes.keys())

			if (mode=='Query'):
				protein=input('Protein:')
				result = minhash.queryProtein(protein)
				jaccResultsDict = minhash.checkJaccardResultsOfProtein(protein, result)
				sorted_jaccResultsDict = OrderedDict(sorted(jaccResultsDict.items(), key=lambda x: -x[1]))
				for jaccRes in sorted_jaccResultsDict.items():
					print(jaccRes[0]," - Jaccard: ",jaccRes[1])
			if (mode=='Save LSH'):
				minhash.saveLSH()
			if (mode=='Load LSH'):
				minhash.loadLSH()

			mode=input('Choose option:')
       

if __name__ == '__main__':
    Analyzer().run()

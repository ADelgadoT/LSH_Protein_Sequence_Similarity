import sys
from ProteinsManager import ProteinsManager
from UniprotDB import UniprotDB
from uniProtein import uniProtein
from lsh import LSH
#from Bio import SeqIO

class Analyzer(object):

	def run(self):

		mode=input('Choose menu options:')
		
		uniDB = UniprotDB("Uniprot_DB.sqlite")
		minhash = LSH(0.42,128)
		
		while(mode!='Exit'):
			print(mode)
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
				minhash = LSH(0.42,128)
				uniDB = UniprotDB("Uniprot_DB.sqlite")
				#uni_DB.close()
				proteins = uniDB.extractProteins()
				#minhash.calculateLSH([protein[1] for protein in proteins])
				minhashes, lsh = minhash.calculateLSH(proteins)
				print(minhashes.keys())


			if (mode=='Query'):
				protein=input('Protein:')
				minhash.queryProtein(protein)
			if (mode=='Save LSH'):
				minhash.saveLSH()
			if (mode=='Load LSH'):
				minhash.loadLSH()

			mode=input('Choose menu options:')
       

if __name__ == '__main__':
    Analyzer().run()

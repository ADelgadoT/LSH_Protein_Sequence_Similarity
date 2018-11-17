import sys
from ProteinsManager import ProteinsManager
from UniprotDB import UniprotDB
from uniProtein import uniProtein
from lsh import LSH
#from Bio import SeqIO

class Analyzer(object):

	def run(self):

		#uniDB = UniprotDB("Uniprot_DB.sqlite")
		#uniDB.close()

		#DB creation:
		protManager = ProteinsManager()
		#uniDB = UniprotDB("Uniprot_DB.sqlite")

		#uniDB.createTables()
		#loadProteins("uniprot_thaliana.xml",uni_DB)
		#protManager.loadProteins("uniprot-filtered-organism_all.xml",uniDB)

		#uniDB = UniprotDB("Uniprot_DB.sqlite")
		#uni_DB.close()
		#uniDB.extractProteinSeqFromSpecieCsv("Arabidopsis thaliana (Mouse-ear cress)")

		uniDB = UniprotDB("Uniprot_DB.sqlite")
		#uni_DB.close()
		proteins = uniDB.extractProteinSeqFromSpecie25("Arabidopsis thaliana (Mouse-ear cress)")

		minhash = LSH(0.5,128)
		minhash.calculateLSH([protein[1] for protein in proteins])


       

if __name__ == '__main__':
    Analyzer().run()
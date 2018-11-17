from Bio import SeqIO
from uniProtein import uniProtein

class ProteinsManager(object):
    
	def loadProteins(self,fileName, uni_DB):
	  descriptions = []
	  record = SeqIO.parse(open(fileName), 'uniprot-xml')
	  #print(len(list(record)))
	  for r in record: 
	      print(r.id)
	      #print(r.name)
	      #print(r.description)
	      #print(r.annotations['organism'])
	      #print(r.annotations['recommendedName_fullName'])
	      #print(r.seq)
	      prot = uniProtein(r.id,r.name,r.description,r.annotations['organism'],r.annotations['recommendedName_fullName'][0], str(r.seq))
	      #print(prot)
	      #prot.printUniProtein()
	      uni_DB.addProtein(prot.id, prot.name, prot.fullName, prot.description, prot.seq, prot.organism)
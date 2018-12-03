from datasketch import MinHash, MinHashLSH
from nltk import ngrams
import pickle

class LSH:

	def __init__(self, threshold=0.5, num_perm=96):
		self.threshold = threshold
		self.num_perm = num_perm
		self.lsh = MinHashLSH(threshold=self.threshold, num_perm=self.num_perm)
		self.minhashes = dict()

	def calculateLSH(self,data,n):
		# Create MinHash objects
		print("Initializing the LSH for %i samples" % len(data))
		for c, i in enumerate(data):
			print("data" + str(i[0]))
			minhash = MinHash(num_perm=self.num_perm)
			for d in ngrams(i[1], n):
				minhash.update("".join(d).encode('utf-8'))
			self.lsh.insert(i[0], minhash)
			self.minhashes[i[0]] = minhash
		print(self.minhashes.keys())
		return self.minhashes, self.lsh

	def queryProtein(self,protein):
		if (self.minhashes.get(protein) is None):
			print("Candidates for ",protein," not found")
			return None
		else:
			result = self.lsh.query(self.minhashes[protein])
			#print("Candidates for ",protein,": ",result)
			return result

	def queryAll(self):
		resultsDict = dict()
		for protein in self.minhashes.keys():
			result = self.lsh.query(self.minhashes[protein])
			resultsDict[protein] = result
		return resultsDict

	def estimateJaccard(self, query, match):
		return self.minhashes[query].jaccard(self.minhashes[match])
		
	def checkJaccardResultsOfProtein(self, query, matches):
		#print(self.minhashes)
		jaccResultsDict = dict()
		for match in matches:
			jaccRes = self.minhashes[query].jaccard(self.minhashes[match])
			jaccResultsDict[match] = jaccRes
		#print("Jaccard values ",query,": ",str(jaccResultsDict.items))
		return jaccResultsDict

	def saveLSH(self, number = 0):
		try:
			with open('lsh%i.pickle' % number, 'wb') as handle:
				pickle.dump(self.lsh, handle, protocol=pickle.HIGHEST_PROTOCOL)
			with open('minhashes%i.pickle' % number, 'wb') as handle:
				pickle.dump(self.minhashes, handle, protocol=pickle.HIGHEST_PROTOCOL)
			print("Saved")
		except:
			print("Error occurred when trying to save LSH")
			return

	def loadLSH(self, number = 0):
		try:
			with open('lsh%i.pickle' % number, 'rb') as handle:
				self.lsh = pickle.load(handle)
			with open('minhashes%i.pickle' % number, 'rb') as handle:
				self.minhashes = pickle.load(handle)
			print("Loaded")
		except:
			print("Error occurred when trying to load LSH")
			return



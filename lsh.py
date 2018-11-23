from datasketch import MinHash, MinHashLSH
from nltk import ngrams
import pickle

class LSH:

	def __init__(self, threshold=0.5, num_perm=128):
		self.threshold = threshold
		self.num_perm = num_perm
		self.lsh = MinHashLSH(threshold=self.threshold, num_perm=self.num_perm)
		self.minhashes = dict()

	def calculateLSH(self,data):
		# Create MinHash objects
		print(len(data))
		for c, i in enumerate(data):
			print("data" + str(i[0]))
			minhash = MinHash(num_perm=128)
			for d in ngrams(i[1], 3):
				minhash.update("".join(d).encode('utf-8'))
			self.lsh.insert(i[0], minhash)
			self.minhashes[i[0]] = minhash
		print(self.minhashes.keys())
		#print(self.lsh)
		#data = pickle.dumps(self.lsh)
		#with open('filename.pickle', 'wb') as handle:
		#	pickle.dump(self.lsh, handle, protocol=pickle.HIGHEST_PROTOCOL)
		#for i in self.minhashes.keys():
		#  result = self.lsh.query(self.minhashes[i])
		#  print("Candidates with Jaccard similarity > 0.5 for input", i, ":", result)
		#print(type(self.minhashes))
		return self.minhashes, self.lsh

	def queryProtein(self,protein):
		#print(self.minhashes)
		result = self.lsh.query(self.minhashes[protein])
		print("Candidates for ",protein,": ",result)
		return result

	def saveLSH(self):
		with open('lsh.pickle', 'wb') as handle:
			pickle.dump(self.lsh, handle, protocol=pickle.HIGHEST_PROTOCOL)
		with open('minhashes.pickle', 'wb') as handle:
			pickle.dump(self.minhashes, handle, protocol=pickle.HIGHEST_PROTOCOL)
		#print(self.minhashes)
		print("Saved")

	def loadLSH(self):
		with open('lsh.pickle', 'rb') as handle:
			self.lsh = pickle.load(handle)
		with open('minhashes.pickle', 'rb') as handle:
			self.minhashes = pickle.load(handle)
		#print(self.lsh)
		#print(self.minhashes)
		print("Loaded")



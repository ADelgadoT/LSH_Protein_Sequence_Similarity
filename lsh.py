from datasketch import MinHash, MinHashLSH
from nltk import ngrams

class LSH:

	def __init__(self, threshold=0.5, num_perm=128):
		self.threshold = threshold
		self.num_perm = num_perm
		self.lsh = MinHashLSH(threshold=self.threshold, num_perm=self.num_perm)

	def calculateLSH(self,data):
		# Create MinHash objects
		minhashes = {}
		for c, i in enumerate(data):
		  minhash = MinHash(num_perm=128)
		  for d in ngrams(i, 3):
		    minhash.update("".join(d).encode('utf-8'))
		  self.lsh.insert(c, minhash)
		  minhashes[c] = minhash
		print(minhashes.keys())
		for i in range(len(minhashes.keys())):
		  result = self.lsh.query(minhashes[i])
		  print("Candidates with Jaccard similarity > 0.5 for input", i, ":", result)


- AUTHORS
		* Anna Delgado Tejedor (s162178)
		* Balbina Virgili Rocosa (s181866)
		* Valentijn Floris Broeken (s172092)


- IMPLEMENTATION
	All of our code has been developed in Python3 except the obtainment of the BLASTp results to verify the reliability of the implemented method, which have been developed with Bash.

- PROJECT STRUCTURE
	.
	|
	├── Report.pdf
	├── README.txt
	│ 
	└── Code
	    ├── benckmarking.py
	    ├── benckmarking_permutations.py
	    ├── lsh.py
	    ├── main.py
	    ├── ProteinsManager.py
	    ├── ResultsDB.py
	    ├── UniprotDB.py
	    ├── uniProtein.py
 	    │
	    ├── Uniprot_DB.sqlite
	    ├── Results_DB.sqlite
 	    │
	    ├── lsh.pickle
	    ├── minhashes.pickle
	    │
	    ├── Ecolx.xml
	    ├── PseA7.xml
	    ├── uniprot_sprot.xml
	    │
	    ├── all_results_nofilter.txt
	    └── query_match_identity_alignmentLength.txt


To execute the main functionalities of our developed code, the following command must be executed.

	$ python3 main.py

With it, the developed interactive command line program will be executed. The main functionalities that can be used are the following ones:

	1-. Load Database or L
	2-. Delete Database or D
	3-. Calculate LSH or C
	4-. Recalculate LSH or RC
	5-. Query LSH or Q
	6-. Query All LSH or A
	7-. Read BLAST or B
	8-. Compare Results or R
	9-. Save LSH or S
	10-. Load LSH or LL
	11-. Exit or E

We consider that most of them are very self-described but if you want to read more details about the functionalities of our program you can read our written report.

The most useful happy-paths for our functionalities are the following ones:

A: Without precomputed results

	> Delete Database
	> Load Database
	> Calculate LSH
	> Save LSH
	> Query
   	    > XXXX
	> Compare results
	> Query
   	    > XXXX
	> Exit

B: With precomputed database

	> Calculate LSH
	> Save LSH
	> Query
   	    > XXXX
	> Compare Results
	> Exit

C: With precomputed results

	> Load LSH
	> Query
           > XXXX
	> Compare results
	> Exit

To be able to execute the benchmarking developed, the following commands must be executed.

	$ python3 benchmarking.py
	$ python3 benchmarking_permutations.py


** WARNING **

We HAVE NOT included the entire SwissProt xml file with all the proteins as it has a weight of more than 6GB. However, one can find it here:
	https://www.uniprot.org/uniprot/

To make things easier, we have already included on the uploaded files the pre-loaded databases with the information of all proteins that we have used on our experiments. These files are Uniprot_DB.sqlite and Results_DB.sqlite. 

Additionally, we have also saved the already calculated LSH results for all the proteins saved on theses databases.

So you are able to directly execute the C path of functionalities to quickly see the results obtained.

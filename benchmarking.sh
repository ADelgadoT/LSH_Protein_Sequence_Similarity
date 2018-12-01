#Benchmarking procedure: 

#Database creation: multi-fasta file is expected as input
makeblastdb -in PseA7_Ecoli.fsa -dbtype prot -out PseA7_Ecoli

#Blastp execution: individual FASTA files are expected as input
time for file in /home/anna_delgado/Documentos/Computational/single_fasta/*.txt 
do
blastp -query $file -db PseA7_Ecoli -out $file.blastp -outfmt 6 
done

#Merge all results: 
cat *.blastp > all_results

#Filter the results based on identity (%) and alignment length: 
perl -ane 'BEGIN{$prev=""};  if ($F[2] > 80 & $F[3] > 100) { print $_}' all_results > filtered_blast.txt


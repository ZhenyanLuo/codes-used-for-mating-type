MACSE v2.07.0
iqtree version 2.1.4-beta
```
for i in Pgt Pca Pst Pt ; do for f in int_nuc*.fasta ; do grep ${i} ${f}|sed -e 's/>//' >>${i}_${f%_aln.fasta}.list ; done ; done

for i in *_bW-HD.list; 
do species=`echo ${i}|cut -d'_' -f1`; 
seqtk subseq int_nuc_bW-HD_aln.fasta ${i} >${species}/${species}_bW-HD.fasta; 
done

for i in *_bE-HD.list; 
do species=`echo ${i}|cut -d'_' -f1`; 
seqtk subseq int_nuc_bE-HD_aln.fasta ${i} >${species}/${species}_bE-HD.fasta; 
done

export PATH=/home/jylin/software/trimal/source:$PATH
for i in Pca Pt Pst Pgt ; 
do      
for f in ${i}/*_bW-HD.fasta ;     
do 
java -jar /home/jylin/miniconda3/envs/mapping/share/macse-2.07-0/macse_v2.07.jar -prog alignSequences -out_NT ${f%.fasta}_macse.aln -seq ${f} ; 
trimal -in ${f%.fasta}_macse.aln -out ${f%.fasta}_trim.fasta  -automated1 -fasta;
java -jar /home/jylin/miniconda3/envs/mapping/share/macse-2.07-0/macse_v2.07.jar -prog alignSequences -out_NT ${f%.fasta}_realn.aln -seq ${f%.fasta}_trim.fasta;
rm ${f%.fasta}_trim.fasta ${f%.fasta}_macse.aln;
done ;  
done

for i in Pca Pt Pst Pgt ; 
do      
for f in ${i}/*_bE-HD.fasta ;     
do 
java -jar /home/jylin/miniconda3/envs/mapping/share/macse-2.07-0/macse_v2.07.jar -prog alignSequences -out_NT ${f%.fasta}_macse.aln -seq ${f} ; 
trimal -in ${f%.fasta}_macse.aln -out ${f%.fasta}_trim.fasta  -automated1 -fasta;
java -jar /home/jylin/miniconda3/envs/mapping/share/macse-2.07-0/macse_v2.07.jar -prog alignSequences -out_NT ${f%.fasta}_realn.aln -seq ${f%.fasta}_trim.fasta;
rm ${f%.fasta}_trim.fasta ${f%.fasta}_macse.aln;
done ;  
done


find . -type f -name "*_realn.aln"|xargs sed -i 's/!/-/g'

for i in Pca Pt Pst Pgt ; do for f in ${i}/*_realn.aln ; do iqtree2 -s ${f} -B 1000 -keep-ident --treels -T AUTO ; done ; done

#Find the best tree with SH test
for i in Pca Pt Pst Pgt ;
do
iqtree2 -s ${i}/${i}_bW-HD_realn.aln -z ${i}/${i}_bW-HD_realn.aln.treels -zb 10000 --keep-ident --prefix ${i}/${i}_bW-HD_test -n 0 -T AUTO --redo
cat ${i}/${i}_bW-HD_test.trees |cut -d'=' -f2|sort -k1,1nr|head -100|cut -d']' -f2|sed -E 's/:[0-9]+(\.[0-9]+)?//g; s/\)[0-9]+(\.[0-9]+)?/)/g'|sort -k1,1 -u > ${i}/${i}_bW-HD.best100
iqtree2 -s ${i}/${i}_bE-HD_realn.aln -z ${i}/${i}_bE-HD_realn.aln.treels -zb 10000 --keep-ident --prefix ${i}/${i}_bE-HD_test -n 0 -T AUTO --redo
cat ${i}/${i}_bE-HD_test.trees |cut -d'=' -f2|sort -k1,1nr|head -100|cut -d']' -f2|sed -E 's/:[0-9]+(\.[0-9]+)?//g; s/\)[0-9]+(\.[0-9]+)?/)/g'|sort -k1,1 -u > ${i}/${i}_bE-HD.best100
cat  ${i}/${i}_bW-HD.best100 ${i}/${i}_bE-HD.best100 >${i}/${i}_HD.to_test.tree
rm ${i}/${i}_HD.to_test.tree
done

#Perform AU-test
for i in Pca Pt Pst Pgt ;
do
aln1=${i}/${i}_bE-HD_realn.aln
aln2=${i}/${i}_bW-HD_realn.aln
trees=${i}/${i}_HD.to_test.trees
iqtree2 -s ${aln1} -z ${tree} -au --keep-ident -pre ${i}/bE-au-test -zb 10000 -n 0
iqtree2 -s ${aln2} -z ${tree} -au --keep-ident -pre ${i}/bW-au-test -zb 10000 -n 0
done

```
both blastn and blastp were used in this step to identify MAT genes, the query sequences were downloaded from Cuomo et al. 2017
```
#Make blastdb for blast
for i in *.faa ; do makeblastdb -dbtype prot -in ${i} ; done
for i in *.cds ; do makeblastdb -dbtype nucl -in ${i} ; done
#Use identified MAT genes as hints to identify proteins from Cuomo paper
for i in *faa ; do blastp -query MAT.faa -db ${i} -outfmt 6 -evalue 0.05 ; done
for i in *faa ; do blastn -query MAT.cds -db ${i} -outfmt 6 -evalue 0.05 ; done
#For detecting small pipetides
hmmbuild --amino  mfa.faa.hmm mfa.faa.sto
#Translate nucleotide sequence into amino acid 
for i in *.fna; do esl-translate ${i} > ${i%.fna}_6frame.faa; done
#Make name shorter
for i in *_6frame.faa; do cat ${i}|cut -d' ' -f1 >${i}.otf; done
#Convert output
for i in *_6frame.faa; do hmmsearch mfa.faa.hmm ${i} > ${i%faa}hmm_out; done
for i in *.hmm_out; do cat ${i}| awk '/Query:/{flag=1;next}/Domain annotation/{flag=0}flag'|grep -v '\--'|awk -OFS='\t' 'NR > 1 { print }'|sed -e 's/frame.*//g' -e 's/Description//g' -e 's/Sequence/ORF/g' -e 's/length=//g' -e 's/E-value/E-value_full/' -e 's/score/score_full/' -e 's/bias/bias_full/' -e 's/ \+ /\t/g' -e 's/ /\t/g' >${i}.tsv; done
#Make a file which have orresponding information of orf 
for i in *_6frame.faa; do grep '>'|sed -e 's/>//g' -e 's/source=//g' -e 's/coords=//g' -e 's/\../\'$'\t/g' -e 's/length.*//g' >${i%.faa}.map ; done
#Convert hmm output into bed file format
for i in *.tsv; do cat ${i}|awk '{print $10 "\t" $11 "\t" $9}'|sed -e 's/\../ /g'|awk '{if($2>$3) print $1"\t"$3"\t"$2"\t"$4"\t""0""\t""-"; else print $1"\t"$2"\t"$3"\t"$4"\t""0""\t""+"}' >${i}_hmm.bed; done
```
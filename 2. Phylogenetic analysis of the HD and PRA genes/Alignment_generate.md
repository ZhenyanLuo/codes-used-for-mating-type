MACSE v2.07.0
trimal v.1.2

Run MACSE->Trimmal->MACSE->Replace '!' generate by MACSE as gaps
```
for i in *.fasta ; do java -jar macse.jar -prog alignSequences -out_NT macse_${i} -seq ${i} ; done
for i in macse_*.fasta ; do trimal -in ${i} -gt 0.9 -fasta -out trim_${i} ; done
for i in trim_macse_*.fasta ; do java -jar /home/jylin/miniconda3/envs/mapping/share/macse-2.07-0/macse_v2.07.jar -prog alignSequences -out_NT realn_${i} -seq ${i} ; done
sed -i 's/!/-/g' realn_*.fasta
```
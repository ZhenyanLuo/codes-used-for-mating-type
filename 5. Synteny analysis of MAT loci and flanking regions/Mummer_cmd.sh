conda activate /mnt/data/dayhoff/home/u6575017/.conda/envs/base_tools
for i in *.cds ; do makeblastdb -in ${i} -dbtype nucl ; done
for i in *.cds ; do blastn -db ${i} -query ${i} -outfmt 6 -evalue 0.05 >${i%.cds}.o6; done
cat *.o6 >MAT.o6

#Identify conserved nucleotide regions with Mummer
nucmer --prefix=PR PR_chr.fna PR_chr.fna --maxmatch
nucmer --prefix=HD HD_chr.fna HD_chr.fna --maxmatch
for i in *.delta; do show-coords -c -l -T  ${i} >${i%.delta}.coords ; done
```
for ID in `cat SRA.list*`;
do
parallel-fastq-dump --sra-id ${ID} --threads 32 \
                --tmpdir ./tmp --outdir fastq --gzip --skip-technical  --readids \
                --read-filter pass --dumpbase --split-files --clip
done

for i in *_pass_1.fastq.gz;
do
fastqc -t 12 ${i} ${i%pass_1.fastq.gz}pass_2.fastq.gz
done
```

```
for i in *_pass_1.fastq.gz;
do
R1=${i}
R2=${i%_pass_1.fastq.gz}_pass_2.fastq.gz
R1_t=t_${R1}
R2_t=t_${R2}
R1_un=un_${R1}
R2_un=un_${R2}
trimmomatic PE -threads 32 -phred33 ${R1} ${R2} trimmed/${R1_t} unpair/${R1_un} trimmed/${R2_t} unpair/${R2_un} ILLUMINACLIP:adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
fastqc -t 32 trimmed/t_${R1} trimmed/t_${R2}
done
```
#Kallisto was used to quantify the expression level of each gene
#kallisto 0.50.0

```
#SRA data of PRJEB12497 was from Pst87/66, Pst134E was used as reference (but HD genes were from Pst87/66)
kallisto index -i Pst134E_fixed.index Pst134E_fixed.cds.fasta
for id in `cat list`; do kallisto quant -i Pst134E_fixed.index -o ${id} -b 100 -t 6 --plaintext t_${id}_pass_1.fastq.gz t_${id}_pass_2.fastq.gz; done
#SRA data of PRJNA396589 was from Pst104E, Pst134E's STE3.2-2 was added to cds file of Pst104E since original STE3.2-2 is incomplete in Pst104E annotation
kallisto index -i Pst104E.index Pst104E.cds.fasta
for id in `cat list`; do kallisto quant -i Pst104E.cds.idx -o ${id} -b 100 -t 6 --plaintext t_${id}_pass_1.fastq.gz t_${id}_pass_2.fastq.gz; done
#Pgt use Pgt 210 as reference
kallisto index -i Pgt210.cds.idx Pgt210.cds.fasta
for id in `cat list`; do kallisto quant -i Pgt210.cds.idx -b 50 -l 200 -s 20 -t 12 --plaintext --single --single-overhang -o ${id} t_${id}_pass_1.fastq.gz; done
#Pca use Pca12NC29 as referece
kallisto index -i Pca12NC29.cds.idx Pca12NC29.cds.fasta
for id in `cat list`; do kallisto quant -i Pca12NC29.cds.idx -o ${id} -b 100 -t 6 --plaintext t_${id}_pass_1.fastq.gz t_${id}_pass_2.fastq.gz; done
```
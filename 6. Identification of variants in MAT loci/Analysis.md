SRA data was downloaded in this way:
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
read1=${i}
read2=${i%_1.fastq.gz}_2.fastq.gz
t_read1=${i%_pass_1.fastq.gz}_t_1.fastq.gz
t_read2=${i%_pass_1.fastq.gz}_t_2.fastq.gz
trimmomatic PE -threads 16 \
        ${read1} ${read2} \
        ${t_read1} tmp/${read1}.tmp \
        ${t_read2} tmp/${read2}.tmp \
        ILLUMINACLIP:adapter.fa:2:30:10:2:keepBothReads \
        LEADING:3 TRAILING:3 MINLEN:36
fastqc -t 16 ${t_read1} ${t_read2}
done
```

```
R1=${i}_t_1.fastq.gz
R1_n=${i}_R1_rename.fastq
R2=${i}_t_2.fastq.gz
R2_n=${i}_R2_rename.fastq
ref=${ref}
module load bbmap
module load seqtk
module load bwa-mem2/2.2.1
module load samtools
module load java/jdk-13.33
module load parallel/20191022
bwa-mem2 index ${ref}
ls *.fastq.gz | parallel 'gunzip {}'
sed -E \"s/^((@|\+)[S,E]RR[^.]+\\.[^.]+)\\.(1|2)/\\1/\" ${R1%.gz} >${R1_n}
sed -E \"s/^((@|\+)[S,E]RR[^.]+\\.[^.]+)\\.(1|2)/\\1/\" ${R2%.gz} >${R2_n}
bwa-mem2 mem -R \"@RG\\tID:${i}\\tSM:${i}\\tPL:ILLUMINA\\tLB:${i}_lib1\" -t 16 ${ref} ${R1_n} ${R2_n} >tmp.sam
cat tmp.sam | samtools sort -O BAM -@ 16 -o ${i}.bam -
samtools index -@ 16 ${i}.bam
samtools stats -@ 16 ${i}.bam > ${i}.stats
```
```
input=${i}
output=${i%.bam}.picard.bam
SRA=`echo ${i}|cut -d'/' -f3|cut -d '.' -f1`
java -jar /g/data/fa63/share/software/picard/build/libs/picard.jar MarkDuplicates -I ${input} -O ${output} -M ${output%.bam}.metrics --REMOVE_DUPLICATES true
samtools index ${output}

for i in `cat picard.bam.list`;
do
species=`echo ${i}|cut -d'/' -f2`
list=`cat ../MAT.list|grep ${species}|sort -k1,1|awk '{print $1}'|uniq`
samtools view -b -h -@ 32 {i}.bam ${chr} >${i}_merged.sub.bam
samtools index -@ 32 ${i}_merged.sub.bam
done

for i in `cat picard.bam.list`;
do
species=`echo ${i}|cut -d'/' -f2`
list=`cat MAT.list|grep ${species}|sort -k1,1|awk '{print $1}'|uniq`
chr=`echo $list`
bam_name=`echo ${i}|cut -d'/' -f3`
samtools view -b -h -q 30 -@ 24 ${bam_name} ${chr} >${bam_name%.bam}_filter.sub.bam
samtools index -@ 24 ${bam_name%.bam}_filter.sub.bam
```

```
#Find variants
freebayes-parallel <(fasta_generate_regions.py \${s}.fna.fai 100000) 32 -f \${s}.fna \${i} > \${i%.bam}.vcf

for s in Pgt Pst Pca Pt;
do
for i in `cat ${s}.bam.list`;
do
#Remove low quality SNPs
bcftools view -e 'QUAL<40' ${i%.bam}.vcf.gz -O z -o ${i%.bam}.Q40.vcf.gz
bcftools index ${i%.bam}.Q40.vcf.gz
done
done

for i in *.vcf ; do bgzip ${i}; bcftools index ${i}.gz ; done

for s in Pgt Pst Pca Pt;
do
for i in `cat ${s}.bam.list`;
do
cat ${s}.fna |bcftools consensus ${i%.bam}.Q40.vcf.gz >${i%.bam}.fasta
done
done


cat PRA.gff |awk -v OFS='\t' '{if($4<$5)print$1,$4-1000,$5+1000,$9;else print $1,$5,$4,$9}' >PRA.bed

for s in Pgt Pst Pca Pt;
do
for i in `cat ${s}.bam.list`;
do
bedtools getfasta -fi ${i%.bam}.fasta -bed PRA.bed -name >${i%.bam}.PRA.fasta
done
done

sed -i 's/\t/=/' name.file

for s in Pgt Pst Pca Pt;
do
for i in `cat ${s}/name.file` ; 
do 
old=`echo ${i}|cut -d'=' -f2`;
new=`echo ${i}|cut -d'=' -f1`;
sed -i "s/>ID=/>${s}_${new}_/" ${s}/${old}.picard_filter.sub.PRA.fasta;
mv ${s}/${old}.picard_filter.sub.PRA.fasta ${s}/${s}_${old}.PRA.fasta
done
```


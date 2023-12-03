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

Generate consensus sequence
```
for i in `cat list.tmp` ; do freebayes-parallel <(fasta_generate_regions.py Leaf_rust_LR.scaffolds.fa.fai 100000) 16 -f Leaf_rust_LR.scaffolds.fa ${i}/Pt/Pt.bwa_mem2.${i}.picard.dup_removed.bam > ${i}/Pt/Pt.bwa_mem2.${i}.picard.dup_removed.vcf ; done
#Only look at specific region
for i in `cat list.tmp` ; do vcfintersect -b Pt76_MAT_10kb.bed ${i}/Pt/Pt.bwa_mem2.${i}.picard.dup_removed.vcf > ${i}/Pt/Pt.bwa_mem2.${i}_MAT_100kb.freebayes.vcf ; done
#Only look at specific region
for i in `cat list.tmp` ; do samtools view -b -L Pt76_MAT_10kb.bed ${i}/Pt/Pt.bwa_mem2.${i}.picard.dup_removed.bam > ${i}/Pt/Pt.bwa_mem2.${i}_MAT_10kb.bam ; done
# call variants
for i in `cat list.tmp` ; do bcftools mpileup -Ou -f Leaf_rust_LR.scaffolds.fa ${i}/Pt/Pt.bwa_mem2.${i}_MAT.bam | bcftools call -mv -Oz -o ${i}/Pt/Pt.bwa_mem2.${i}_MAT_10kb_calls.vcf.gz ; done
for i in `cat list.tmp` ; do bcftools index ${i}/Pt/Pt.bwa_mem2.${i}_MAT_10kb_calls.vcf.gz ; done
# apply variants to create consensus sequence
for i in `cat list.tmp` ; do cat Leaf_rust_LR.scaffolds.fa | bcftools consensus -H 1 ${i}/Pt/Pt.bwa_mem2.${i}_MAT_10kb_calls.vcf.gz > ${i}/Pt/Pt.bwa_mem2.${i}_MAT_10kb_consensus.fna ; samtools faidx ${i}/Pt/Pt.bwa_mem2.${i}_MAT_10kb_consensus.fna ; done
#Get specific region
for i in `cat list.tmp` ; do bedtools getfasta -fi ${i}/Pt/Pt.bwa_mem2.${i}_MAT_10kb_consensus.fna -bed Pt_76_MAT_10kb.bed > ${i}/Pt/Pt.bwa_mem2.${i}_MAT_10kb_consensus_sub.fna ; done
#Make a seperate folder
for i in `cat list.tmp` ; do mkdir -p ${i}/Pt/separated_chr ; ln ${i}/Pt/Pt.bwa_mem2.${i}_MAT_10kb_consensus_sub.fna ${i}/Pt/separated_chr/MAT_10kb_consensus_sub.fna ; done
for i in `cat list.tmp` ; do sed -i "s/chr/${i}_chr/g" ${i}/Pt/separated_chr/MAT_10kb_consensus_sub.fna ; done
for i in `cat list.tmp` ; do cat ${i}/Pt/separated_chr/MAT_10kb_consensus_sub.fna >> Mapping/merged_MAT_10kb_consensu_sub.fna; done
```

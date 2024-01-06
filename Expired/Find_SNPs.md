
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
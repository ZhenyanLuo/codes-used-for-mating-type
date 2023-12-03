Juicerbox



Since Pgt and Pst have published centromeric region, here only do for Pt and Pca, both haplotypes of both species extractly do in the same way

Restriction enyzeme is double checked by mapping Hi_C reads to genome to check cut site
```
set -xue
cd /g/data/xf3/zl1602/Pca_HiC
module load bwa-mem2/2.2.1
bwa-mem2 mem -t 24 Pca203_hapA.fa henningsen-203_S3HiC_R1.fastq.gz  henningsen-203_S3HiC_R2.fastq.gz -o Pca203_hapA_hic.sam
module load samtools
samtools view -@ 24 -bS Pca203_hapA_hic.sam |samtools sort - -@ 24 -n -o Pca203_hapA_hic.bam
samtools index -@ 24 Pca203_hapA_hic.bam

```

```
#Prepare input

git clone https://github.com/aidenlab/juicer.git --branch 1.6
ln -s ../juicer/CPU scripts
cd scripts/common/
wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
ln -s juicer_tools.1.9.9_jcuda.0.8.jar  juicer_tools.jar

#Prepare splited Hi-C reads
split -a 3 -l 90000000 -d --additional-suffix=_R2.fastq ../fastq/henningsen-203_S3HiC_R2.fastq
split -a 3 -l 90000000 -d --additional-suffix=_R1.fastq ../fastq/henningsen-203_S3HiC_R1.fastq
#Get chromosome length
bioawk -c fastx '{print $name"\t"length($seq)}' references/Pca203_hapA.fasta >chrom.sizes
#Generate restriction site
python3 /g/data/xf3/miniconda/share_software/juicer-1.6/misc/generate_site_positions.py Sau3AI Pca203_hapA.fasta Pca203_hapA.fasta
```

```
set -xue
module load java/jdk-8.40
module load python2/2.7.17
module load parallel/20191022
module load bwa/0.7.17
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate juicer
cd /scratch/xf3/zl1602/HiC-2/Pca/Pca_hapA
export wdir=/scratch/xf3/zl1602/HiC-2/Pca/Pca_hapA
scripts/juicer.sh -S early -D ${wdir} -g ${wdir}/references/Pca203_hapA.fa -z ${wdir}/references/Pca203_hapA.fa -y ${wdir}/restriction_sites/Pca203_hapA.fa_Sau3AI.txt -p $wdir/chrom.sizes -t 32
```

```
set -xue
module load java/jdk-8.40
module load python2/2.7.17
module load parallel/20191022
source /g/data/xf3/miniconda/etc/profile.d/conda.sh
conda activate 3d-dna
export JAVA_TOOL_OPTIONS="-Xmx150g"
export PATH=$PATH:/g/data/xf3/zl1602/Pt_HiC/Jucier/3d-dna
cd /scratch/xf3/zl1602/HiC-2/Pca/Pca_hapA
run-asm-pipeline.sh -r -1  -m haploid --build-gapped-map references/Pca203_hapA.fa aligned/merged_nodups.txt
```
.0.hic .0.asm and .0.assembly files were used for visualizing Hi-C heatmap in Juicerbox
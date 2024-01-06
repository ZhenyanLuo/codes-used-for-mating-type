#Find SNPs
```
#Generate draft assembly
cd ${wd}
module load fastp
module load bwa-mem2/2.2.1
module load samtools
module load fastqc
module load java/jdk-13.33
module load parallel/20191022
module load python3
export JAVA_HOME=/g/data/fa63/share/miniconda3/envs/SRA
export JAVA_LD_LIBRARY_PATH=/g/data/fa63/share/miniconda3/envs/SRA/lib/server
export PATH=/g/data/fa63/share/miniconda3/bin:/g/data/fa63/share/miniconda3/condabin:/home/106/zl1602/.local/bin:/home/106/zl1602/bin:/opt/pbs/default/bin:/opt/nci/bin:/opt/bin:/opt/Modules/v4.3.0/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/g/data/fa63/share/software/SPAdes-3.15.5-Linux/bin
export OMP_NUM_THREADS=16
R1=${i}_1.fastq.gz
R2=${i}_2.fastq.gz
PBS_d='$PBS_JOBFS'
cp ${R1} ${R2} ${PBS_d}
cd ${PBS_d}
#For samples which average length less than 101, use 71 bp k-mer instead
spades.py -t 16 -k 101 -1 ${R1} -2 ${R2} -o ${i}_spades
cd ${PBS_d}/${i}_spades
cp contigs.fasta ${wd}/${i}_contigs_r1.fasta
cp scaffolds.fasta ${wd}/${i}_scaffolds_r1.fasta
```

```
set -xue
cd ${wd}
export JAVA_HOME=/g/data/fa63/share/miniconda3/envs/SRA
export JAVA_LD_LIBRARY_PATH=/g/data/fa63/share/miniconda3/envs/SRA/lib/server
source /g/data/fa63/share/miniconda3/etc/profile.d/conda.sh
conda activate SRA
export PATH=$PATH:/g/data/fa63/share/software/SPAdes-3.15.5-Linux/bin
export OMP_NUM_THREADS=4
R1=${i}_1.fastq.gz
R2=${i}_2.fastq.gz
contig=${i}_contigs_r1.fasta
PBS_d='$PBS_JOBFS'
cd ${wd}
cp ${contig} ${R1} ${R2} ${PBS_d}
cp HD.fasta ${PBS_d}
cd ${PBS_d}
module load blast
module load bbmap
module load seqtk
module load bwa-mem2/2.2.1
module load samtools
module load java/jdk-13.33
module load parallel/20191022
module load python3
module load bedtools
makeblastdb -in ${contig} -dbtype nucl
blastn -db ${contig} -query HD.fasta -outfmt 6 -evalue 0.01 > ${i}_spades_r1_contigs.blastn
samtools faidx ${contig}
cat ${i}_spades_r1_contigs.blastn|awk -v OFS=\"\\\\t\" '{if (\$10 < \$9) print \$2, \$10, \$9; else print \$2, \$9, \$10}' |sort -k1,1 -k2,2n >tmp.bed1
bedtools merge -i tmp.bed1 >tmp.bed2
bedtools slop -i tmp.bed2 -g ${i}_contigs_r1.fasta.fai -b 1500 >tmp.bed3
bedtools merge -d 6000 -i tmp.bed3 >tmp.bed4
bedtools getfasta -fi ${i}_contigs_r1.fasta -bed tmp.bed4 >tmp.fasta
cat tmp.fasta  HD.fasta >${i}_spades_r1_contigs_HD.fasta
ls *.fastq.gz | parallel 'gunzip {}'
sed -E \"s/^((@|\+)[S,E]RR[^.]+\\.[^.]+)\\.(1|2)/\\1/\" ${R1%.gz} | cut -d' ' -f1 | sed -e '1~2 s/$//g' > ${i}_R1_rename.fastq
sed -E \"s/^((@|\+)[S,E]RR[^.]+\\.[^.]+)\\.(1|2)/\\1/\" ${R2%.gz} | cut -d' ' -f1 | sed -e '1~2 s/$//g' > ${i}_R2_rename.fastq
conda activate /scratch/fa63/zl1602/KAT
kat filter kmer -t 4 ${i}_spades_r1_contigs_HD.fasta
kat filter seq -t 4 -T 0.2 --seq ${i}_R1_rename.fastq --seq2 ${i}_R2_rename.fastq kat.filter.kmer-in.jf27
gzip kat.filter.kmer.in.R1.fastq
gzip kat.filter.kmer.in.R2.fastq
#For samples which average length less than 101, use 71 bp k-mer instead
spades.py -k 101 -1 kat.filter.kmer.in.R1.fastq.gz -2 kat.filter.kmer.in.R2.fastq.gz -o ${i}_spades_kat
mkdir -p ${wd}/KAT/${i}
rm ${i}_R1_rename.fastq
rm ${i}_R2_rename.fastq
cp kat.filter.kmer.in.R1.fastq.gz kat.filter.kmer.in.R2.fastq.gz ${wd}/KAT/${i}
cd ${i}_spades_kat
cp contigs.fasta ${wd}/${i}_contigs_r2.fasta
cp scaffolds.fasta ${wd}/${i}_scaffolds_r2.fasta
```
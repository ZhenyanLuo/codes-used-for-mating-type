#make draft assembly
for i in `cat list`
do
wd=/scratch/fa63/zl1602/SRA/Pst/trimmed/
basic_cmd="#!/bin/bash
#PBS -q normal
#PBS -l mem=100GB
#PBS -l walltime=48:00:00
#PBS -l ncpus=16
#PBS -l jobfs=200GB
#PBS -l wd
#PBS -l storage=scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
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
"
R1=${i}_1.fastq.gz
R2=${i}_2.fastq.gz
PBS_d='$PBS_JOBFS'
cmd="
cp ${R1} ${R2} ${PBS_d}
cd ${PBS_d}
spades.py -t 16 -k 101 -1 ${R1} -2 ${R2} -o ${i}_spades
cd ${PBS_d}/${i}_spades
cp contigs.fasta ${wd}/${i}_contigs_r1.fasta
cp scaffolds.fasta ${wd}/${i}_scaffolds_r1.fasta
"
echo -e "${basic_cmd}\n${cmd}" >${i}_spades.pbs.sh
done









#create kmer matrix from library with KAT
for i in `cat list`
do
wd=/scratch/fa63/zl1602/SRA/Pgt/Pgt_trimmed2/Pgt_trimmed
basic_cmd="#!/bin/bash
#PBS -q normal
#PBS -l mem=100GB
#PBS -l walltime=12:00:00
#PBS -l ncpus=4
#PBS -l jobfs=400GB
#PBS -l wd
#PBS -l storage=scratch/fa63+gdata/fa63
#PBS -P xf3
set -xue
cd /scratch/fa63/zl1602/SRA/Pca/trimmed
export JAVA_HOME=/g/data/fa63/share/miniconda3/envs/SRA
export JAVA_LD_LIBRARY_PATH=/g/data/fa63/share/miniconda3/envs/SRA/lib/server
source /g/data/fa63/share/miniconda3/etc/profile.d/conda.sh
conda activate SRA
export PATH=$PATH:/g/data/fa63/share/software/SPAdes-3.15.5-Linux/bin
export OMP_NUM_THREADS=4
"
R1=${i}_1.fastq.gz
R2=${i}_2.fastq.gz
contig=${i}_contigs_r1.fasta
PBS_d='$PBS_JOBFS'
cmd="
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
spades.py -k 101 -1 kat.filter.kmer.in.R1.fastq.gz -2 kat.filter.kmer.in.R2.fastq.gz -o ${i}_spades_kat
mkdir -p ${wd}/KAT/${i}
rm ${i}_R1_rename.fastq
rm ${i}_R2_rename.fastq
cp kat.filter.kmer.in.R1.fastq.gz kat.filter.kmer.in.R2.fastq.gz ${wd}/KAT/${i}
cd ${i}_spades_kat
cp contigs.fasta ${wd}/${i}_contigs_r2.fasta
cp scaffolds.fasta ${wd}/${i}_scaffolds_r2.fasta
"
echo -e "${basic_cmd}\n${cmd}" >${i}_kat.pbs.sh
done

sed -i 's/\t/=/' label
 for i in *_scaffolds_r2.fasta ; 
 do samtools faidx ${i}; cat ${i}.fai|sort -k1,1 -k2,2n|awk '{if ($2 >=2000)print $1}' >${i%_scaffolds_r2.fasta}.long.list; 
 seqtk subseq ${i} ${i%_scaffolds_r2.fasta}.long.list >${i%_scaffolds_r2.fasta}.long.fasta  ; 
 done

 for i in `cat label` ;
 do old=`echo ${i}|cut -d'=' -f2`; new=`echo ${i}|cut -d'=' -f1`; sed -i "s/NODE/Pgt_${new}/" ${old}.long.fasta ;
  mv  ${old}.long.fasta Pgt_${new}.long.fasta; 
  done
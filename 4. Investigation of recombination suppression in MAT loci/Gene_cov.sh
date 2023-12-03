#run with bash calculate_gene_cov.sh {genome file} {annotation file} {window size you want}
#Make a genomefile which contains the length of each contig
samtools faidx $1 
awk 'BEGIN{FS=OFS="\t"}''{print $1,$2}' $1.fai > tmp.genomefile
#Make windows by using bedtools makewindows
bedtools makewindows -g tmp.genomefile -w $3 > tmp.windows
#Only get gene from the annotation file
awk 'BEGIN{FS=OFS="\t"}''{if($3=="gene") print $1,$4,$5}' $2|grep -v '#'|sort -k1,1 -k2,2n >tmp_sorted.bed
#Calculate the coverage of genes in each window
bedtools coverage -a tmp.windows -b tmp_sorted.bed >$1.genedensity
rm tmp.genomefile tmp.windows tmp_sorted.bed
#samples of output
#chr1  0   100  3  30  100 0.3000000
#$1 is contig name, $2 is start on the contig, $3 is end, $4 is number of genes in this window, $5 is total length of all genes in this window, $6 is wondow size, $7 is coverage rate
#Use $1, $2/$3 and $7 to make density/coverage plot
#This script can also use to plot repeat density but need to use bedtools merge and modify the awk step. 

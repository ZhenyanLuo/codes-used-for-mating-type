#Create soft link
for i in *.gff3; do ln -s $PWD/${i} ${PWD}/analysis/${i} ; done
#Sort and merge repet with maximum distance between distance
for i in *_repet_o.gff3; do cat ${i}|sort -k1,1 -k4,4n|awk -v OFS='\t' '{print $1,$4,$5,".",".",$7}'|sort -k1,1 -k2,2n >${i%_repet_o.gff3}.tmp; bedtools merge -d 10 -s -i ${i%_repet_o.gff3}.tmp >${i%_repet_o.gff3}.merged ; done
#Use bedtools window to split genome
for i in *.genomefile; do bedtools makewindows -g ${i} -w 10000 >${i%.genomefile}.windows ; done
for i in *.merged ; do bedtools merge -i ${i} >${i%.merged}_nostrand.merged ; done
for i in *_nostrand.merged; do bedtools coverage -nonamecheck -a ~/gggenome/genome/genomefile/${i%_nostrand.merged}.windows -b ${i} > ${i%_nostrand.merged}.repeatdensity ; done

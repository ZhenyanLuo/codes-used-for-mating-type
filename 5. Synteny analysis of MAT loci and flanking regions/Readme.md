This subfolder contains scripts used for generating main Fig4, and other supporting information related to synteny analysis.

Repeat information are from the previous folder (4. Investigation of recombination suppressed region), all nucleotide sequences conserved genome in this part of study were from nucmer, with maxmatch tag. All conserved genes identified in this analysis used the same setting from blastn: -evalue 0.5, conserved gene pairs were screened by only reserve gene pairs have at least 70% hits. 

show-coords function in mummer was used to convert delta file to human-readable coords file. Macrosynteny plot was generated after filtering matching shorter than 1000 bp or identity lower than 90%. Please cite the plotting script if you use it for somewhere else, thanks. 

Synteny plot were generated with gggenome, with 40 genes flanking the target gene (HD/PRA/STE3.2-1).

Note: 
Coverage of each TE orders were calculated from repet original output, with only keep column three with 'match' tag, since 'match part' contains TE fragments which overlay 'match' tag. Due to repet pipeline generate lot's of unclassified output, an R code was used to classified these unclassified/have more than one potential classification based on the hits from hits identified by REPET from RepBase 22.05, if hits to reference TE sequences with more than 70%, this record will be classfied accordingly, mroe detail please check gggenome_4fungi.ipynb (chunk 12). 
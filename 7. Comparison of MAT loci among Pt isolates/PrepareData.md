Since only gff3 files of Pt 19NSW04 Pt 20QLD87 available,  gffread from Evidencemodeler was used to extract cds and proteins


```
gffread -g 20QLD87.scaffolds.fasta -x rename/20QLD87.cds -C -W -y rename/20QLD87.protein 20QLD87.gff3
gffread -g 19NSW04.scaffolds.fasta -x rename/19NSW04.cds -C -W -y rename/19NSW04.protein 19NSW04.gff3
```

Make whole-genome alignment with nucmer
```
nucmer --prefix=Pt_HD Pt_HD_chr.fna Pt_HD_chr.fna --maxmatch
nucmer --prefix=Pt_PR Pt_PR_chr.fna Pt_PR_chr.fna --maxmatch
for i in *.delta; do show-coords -c -l -T  ${i} >${i%.delta}.coords ; done
```
Identify MAT genes
```
for i in Pt*.cds ; do makeblastdb -in ${i} -dbtype nucl ; done
for i in Pt*.cds ; do blastn -db ${i} -query MAT_ref.cds -outfmt 6 -evalue 0.05 >${i%.cds}_blastn; done
```

RepeatMasker was used to annotate repeats using Pt_refTEs.fa (ref database from TEdenovo)
To make the analysis consistent, Pt76 was reannotated with the same database one more time
```
RepeatMasker -pa 24 -e rmblast -gff -xsmall -lib Pt76_refTEs.fa 19NSW04_chr_upper.fasta
RepeatMasker -pa 24 -e rmblast -gff -xsmall -lib Pt76_refTEs.fa 20QLD87_chr_upper.fasta
RepeatMasker -pa 24 -e rmblast -gff -xsmall -lib Pt76_refTEs.fa Pt15_chr_upper.fasta
RepeatMasker -pa 24 -e rmblast -gff -xsmall -lib Pt76_refTEs.fa Pt76_chr_upper.fasta
```


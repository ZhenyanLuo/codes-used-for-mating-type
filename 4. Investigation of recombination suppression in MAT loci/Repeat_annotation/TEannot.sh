REPET v3.0

#link the repeat database built in TEdenovo step
ln -s $1_Blaster_GrpRecPil_Map_TEclassif_Filtered/$1_sim_denovoLibTEs_filtered.fa $1_refTEs.fa
#Step1
srun TEannot.py -P $1 -C TEannot.cfg -S 1 > annot_step1.log
#Step2
srun TEannot.py -P $1 -C TEannot.cfg -S 2 -a BLR > annot_step2.log
srun TEannot.py -P $1 -C TEannot.cfg -S 2 -a RM >> annot_step2.log
#srun TEannot.py -P $1 -C TEannot.cfg -S 2 -a CEN >> annot_step2.log
srun TEannot.py -P $1 -C TEannot.cfg -S 2 -a BLR -r >> annot_step2.log
srun TEannot.py -P $1 -C TEannot.cfg -S 2 -a RM -r >> annot_step2.log
#srun TEannot.py -P $1 -C TEannot.cfg -S 2 -a CEN -r >> annot_step2.log
#Step3
srun TEannot.py -P $1 -C TEannot.cfg -S 3 -c BLR+RM > annot_step3.log
#Step4
srun TEannot.py -P $1 -C TEannot.cfg -S 4 -s TRF > annot_step4.log
srun TEannot.py -P $1 -C TEannot.cfg -S 4 -s Mreps >> annot_step4.log
srun TEannot.py -P $1 -C TEannot.cfg -S 4 -s RMSSR >> annot_step4.log
#Step5
srun TEannot.py -P $1 -C TEannot.cfg -S 5 > annot_step5.log
#Step7
srun TEannot.py -P $1 -C TEannot.cfg -S 7 >annot_step7.log
#Step8
srun TEannot.py -P $1 -C TEannot.cfg -S 8 -o GFF3 >annot_step8.log
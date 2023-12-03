#Repet v3.0

ln -s /share/banks/repbase20.05_aaSeq_cleaned_TE.fsa . 
ln -s /share/banks/repbase20.05_ntSeq_cleaned_TE.fsa .
ln -s /share/banks/ProfilesBankForREPET_Pfam27.0_GypsyDB.hmm .
#Run PreProcess.py
PreProcess.py -i {genome} -S 1 
launch_TEdenovo.py -P {genome name} -C TEdenovo.cfg -v 3 > launchTEdenovo.txt

for c in $(cat /fastdata/bop20pp/supergene_AB_maleref/csvs/Zand1to10.csv) ; do  for i in 2-4 3-6 7-12; do  qsub c1.phaser.sh $c $i; done; done

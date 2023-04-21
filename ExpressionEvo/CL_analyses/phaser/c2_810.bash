for c in $(cat /fastdata/bop20pp/supergene_AB_maleref/csvs/Zand1to10.csv) ; do  for i in 10-15 8-13; do  qsub c2.phaser_ase.sh $c $i; done; done

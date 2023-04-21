for c in $(cat /fastdata/bop20pp/supergene_AB_maleref/chicken/csvs/Zand1to10.csv) ; do  for i in 3-6 4-8 5-9 6-11; do  qsub c2.phaser_ase.sh $c $i; done; done

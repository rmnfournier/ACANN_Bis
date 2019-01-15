#!/bin/bash
my_folder='harm_osc'
nb_simulation=`ls $my_folder/sim* | wc -l`
echo "$nb_simulation"
let "nb_simulation = nb_simulation-1" 
for i in `seq 0 $nb_simulation`
do
    line=`head -n 1 $my_folder/simulation$i.csv`
    echo $line>>simulation_$my_folder.csv 
done

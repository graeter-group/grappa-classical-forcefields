#!/bin/bash

# only accounting for
# sidechain-sidechain, Ace-sidechain and sidechain-Nme interactions
# not Nter-Cter, Cter-sidechain, Nter-sidechain,...

AALETTERS=('A' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'K' 'L' 'M' 'N' 'P' 'Q' 'R' 'S' 'T' 'V' 'W' 'Y' 'O' 'J' '1' '2' '3' '4' '5' '6')
AALETTERS_ter=('A' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'K' 'L' 'M' 'N' 'P' 'Q' 'R' 'S' 'T' 'V' 'W' 'Y' '1' '2') 



root_dir="datasets_interactions"
ss=$root_dir/'ss-interaction'
ns=$root_dir/'ns-interaction'
cs=$root_dir/'cs-interaction'
nc=$root_dir/'nc-interaction'


mkdir $ss $ns $cs $nc -p

for ((i=0; i<${#AALETTERS[@]}; i++))
do
    for ((j=0; j<${#AALETTERS[@]}; j++))
    do
        #ss
        pepgen "L${AALETTERS[i]}${AALETTERS[j]}L" $ss/"${AALETTERS[i]}${AALETTERS[j]}" --pdb_only -r
    done
    #ns cap
    pepgen "${AALETTERS[i]}L" $ns/"B${AALETTERS[i]}" --pdb_only -r 
    #cs cap
    pepgen "L${AALETTERS[i]}" $cs/"${AALETTERS[i]}Z" --pdb_only -r
done

for ((i=0; i<${#AALETTERS_ter[@]}; i++))
do
    for ((j=0; j<${#AALETTERS_ter[@]}; j++))
    do
        #nc
        pepgen "${AALETTERS_ter[i]}${AALETTERS_ter[j]}" $nc/"${AALETTERS_ter[i]}${AALETTERS_ter[j]}" --pdb_only -r -t
    done
    #nc cap
    pepgen "B${AALETTERS_ter[i]}" $ns/"B${AALETTERS_ter[i]}" --pdb_only -r -t
    #nc cap
    pepgen "${AALETTERS_ter[i]}Z" $ns/"${AALETTERS_ter[i]}Z" --pdb_only -r -t
done

for ((i=0; i<${#AALETTERS[@]}; i++))
do
    for ((j=0; j<${#AALETTERS_ter[@]}; j++))
    do
        #cs
        pepgen "BL${AALETTERS[i]}${AALETTERS_ter[j]}" $cs/"${AALETTERS[i]}${AALETTERS_ter[j]}" --pdb_only -r -t
    done
done

for ((i=0; i<${#AALETTERS_ter[@]}; i++))
do
    for ((j=0; j<${#AALETTERS[@]}; j++))
    do
        #ns
        pepgen "${AALETTERS_ter[i]}${AALETTERS[j]}LZ" $ns/"${AALETTERS_ter[i]}${AALETTERS[j]}" --pdb_only -r -t        
    done
done

#nc cap cap
pepgen "BZ" $nc/"BZ" --pdb_only -r -t


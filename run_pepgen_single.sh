#!/bin/bash

AALETTERS=('A' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'K' 'L' 'M' 'N' 'P' 'Q' 'R' 'S' 'T' 'V' 'W' 'Y' 'O' 'J' '1' '2' '3' '4' '5' '6')

if [ $# -eq 0 ]; then
    # no argument supplied
    root_dir="datasets"
    sidechain=$root_dir/'sidechain'
    nter=$root_dir/'nter'
    cter=$root_dir/'cter'
    caps=$root_dir/'caps'
elif [ $# -eq 5 ]; then
    root_dir=$1
    sidechain=$root_dir/$2
    nter=$root_dir/$3
    cter=$root_dir/$4
    caps=$root_dir/$5
else
    echo Unexpected number of arguments $#
    exit 1
fi


mkdir $sidechain $nter $cter $caps -p

for ((i=0; i<${#AALETTERS[@]}; i++))
do
    pepgen "L${AALETTERS[i]}L" $sidechain/"${AALETTERS[i]}" --pdb_only -r
    pepgen "${AALETTERS[i]}"LZ $nter/"${AALETTERS[i]}" -t --pdb_only -r
    pepgen BL"${AALETTERS[i]}" $cter/"${AALETTERS[i]}" -t --pdb_only -r

done

pepgen LL $caps/B --pdb_only -r
pepgen LL $caps/Z --pdb_only -r
#!/bin/bash

#hyperparameters
grappa_tag='grappa-1.3.0'
charge_model='amber99'

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

#run
outputdir=topologies_${grappa_tag}_${charge_model}/
echo $outputdir
mkdir -p $outputdir

for datasetdir in $sidechain $nter $cter $caps
do
    datasetdir_arr=("$datasetdir")
    dataset_name=$(basename "$datasetdir") 
    mkdir -p $outputdir/$dataset_name
    for topdir in $datasetdir_arr/*    
    do
        echo $topdir
        seq=$(basename "$topdir")
        echo $seq

        grappa_gmx -f $topdir/topol.top -o $outputdir/$dataset_name/$seq.top -t $grappa_tag 
    done
done
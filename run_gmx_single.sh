#!/bin/bash

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
# Glob all files in the directory
subdirs=("$sidechain"/* "$nter"/* "$cter"/* "$caps"/*)

# Print the list of files
for subdir in "${subdirs[@]}"
do
    echo $subdir
    cd $subdir
    echo "$subdir"
    cp -r /hits/fast/mbm/hartmaec/workdir/FF99SBILDNPX_OpenMM/amber99sb-star-ildnp_mod.ff .
    cp amber99sb-star-ildnp_mod.ff/residuetypes.dat .
    mol=$(basename "$subdir")
    lys=''
    arg=''
    asp=''
    glu=''
    gln=''
    his=''
    for (( i=0; i<${#mol}; i++ ));
        do

        aa="${mol:i:1}"
        echo $aa
        
        if [[ "$aa" == "1" ]]; then
            #HID
            his+="0\n"
        elif [[ "$aa" == "2" ]]; then
            #HIP
            his+="2\n"
        elif [[ "$aa" == "3" ]]; then
            #GLH
            glu+="1\n"
        elif [[ "$aa" == "4" ]]; then
            #ASH
            asp+="1\n"
        elif [[ "$aa" == "5" ]]; then
            #LYN
            lys+="0\n"
        elif [[ "$aa" == "D" ]]; then 
            #Asp
            asp+="0\n"
        elif [[ "$aa" == "E" ]]; then
            glu+="0\n"
        elif [[ "$aa" == "Q" ]];then
            gln+="0\n"
        elif [[ "$aa" == "H" ]];then
            #HIE 
            his+="1\n"
        elif [[ "$aa" == "K" ]];then 
            lys+="1\n"       
        elif [[ "$aa" == "R" ]]; then
            arg+="1\n"
        fi
    done
    gmx_pipe="1\n1\n"$lys$arg$gln$asp$glu$his
    echo $gmx_pipe
    printf "${gmx_pipe} " | gmx pdb2gmx -f pep.pdb -ignh -inter
    [ -s topol.top ] || rm topol.top
    cd /hits/fast/mbm/hartmaec/workdir/FF99SBILDNPX_OpenMM
done

#!/bin/bash

#run paramaters
root_dir="datasets"
sidechain='sidechain'
nter='nter'
cter='cter'
caps='caps'

eval "$(/hits/fast/mbm/hartmaec/sw/conda/bin/conda shell.bash hook)"
conda activate pepgen

echo "** Starting pepgen **"
bash run_pepgen_single.sh $root_dir $sidechain $nter $cter $caps >> out.log
echo pdbs number of files: sidechain:$(ls -lt $root_dir/$sidechain/*/*pdb | wc -l)/28, nter: $(ls -lt $root_dir/$nter/*/*pdb | wc -l)/28, cter: $(ls -lt $root_dir/$cter/*/*pdb | wc -l)/28, cap: $(ls -lt $root_dir/$caps/*/*pdb | wc -l)/2 >> eval_run.txt

echo "** Starting gmx **"
bash run_gmx_single.sh $root_dir $sidechain $nter $cter $caps >> out.log
echo top number of files: sidechain:$(ls -lt $root_dir/$sidechain/*/*top | wc -l)/28, nter: $(ls -lt $root_dir/$nter/*/*top | wc -l)/28, cter: $(ls -lt $root_dir/$cter/*/*top | wc -l)/28, cap: $(ls -lt $root_dir/$caps/*/*top | wc -l)/2 >> eval_run.txt

echo "** Starting grappa **"
conda activate grappa_cpu
bash run_grappa_single.sh $root_dir $sidechain $nter $cter $caps >> out.log
outputdir_expected=topologies_grappa-1.3.0_amber99
echo grappa top number of files: sidechain:$(ls -lt $outputdir_expected/$sidechain/*top | wc -l)/28, nter: $(ls -lt $outputdir_expected/$nter/*top | wc -l)/28, cter: $(ls -lt $outputdir_expected/$cter/*top | wc -l)/28, cap: $(ls -lt $outputdir_expected/$caps/*top | wc -l)/2 >> eval_run.txt

echo "** Writing dataset configs **"
printf "type: single\nresid: 3\nres: seq\nresname_prefix: None\n" > $outputdir_expected/$sidechain/dataset.yml
printf "type: single\nresid: 1\nres: seq\nresname_prefix: N\n" > $outputdir_expected/$nter/dataset.yml
printf "type: single\nresid: 3\nres: seq\nresname_prefix: C\n" > $outputdir_expected/$cter/dataset.yml
printf "type: single\nresid: cap\nres: seq\nresname_prefix: None\n" > $outputdir_expected/$caps/dataset.yml

echo "Done!"
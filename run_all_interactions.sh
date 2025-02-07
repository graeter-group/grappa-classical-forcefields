#!/bin/bash

#run paramaters
root_dir="datasets"
sidechain='ss-interaction'
nter='ns-interaction'
cter='cs-interaction'

eval "$(/hits/fast/mbm/hartmaec/sw/conda/bin/conda shell.bash hook)"
conda activate pepgen

echo "** Starting pepgen **"
bash run_pepgen_interactions.sh $root_dir $sidechain $nter $cter >> out.log
echo pdbs number of files: sidechain:$(ls -lt $root_dir/$sidechain/*/*pdb | wc -l)/784, nter: $(ls -lt $root_dir/$nter/*/*pdb | wc -l)/28, cter: $(ls -lt $root_dir/$cter/*/*pdb | wc -l)/28 >> eval_interactions.txt

echo "** Starting gmx **"
bash run_gmx_interactions.sh $root_dir $sidechain $nter $cter >> out.log
echo top number of files: sidechain:$(ls -lt $root_dir/$sidechain/*/*top | wc -l)/784, nter: $(ls -lt $root_dir/$nter/*/*top | wc -l)/28, cter: $(ls -lt $root_dir/$cter/*/*top | wc -l)/28 >> eval_interactions.txt

echo "** Starting grappa **"
conda activate grappa_cpu
bash run_grappa_interactions.sh $root_dir $sidechain $nter $cter >> out.log
outputdir_expected=topologies_grappa-1.3.0_amber99
echo grappa top number of files: sidechain:$(ls -lt $outputdir_expected/$sidechain/*top | wc -l)/784, nter: $(ls -lt $outputdir_expected/$nter/*top | wc -l)/28, cter: $(ls -lt $outputdir_expected/$cter/*top | wc -l)/28 >> eval_interactions.txt

echo "** Writing dataset configs **"
printf "type: interaction\nresid: 3,4\nres: None\nresname_prefix: None\n" > $outputdir_expected/$sidechain/dataset.yml
printf "type: interaction\nresid: 1,2\nres: None\nresname_prefix: None\n" > $outputdir_expected/$nter/dataset.yml
printf "type: interaction\nresid: 3,4\nres: None\nresname_prefix: None\n" > $outputdir_expected/$cter/dataset.yml

echo "Done!"
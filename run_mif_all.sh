runlvl=$1
mkdir output/
echo "Running at run level $runlvl :)"
Rscript scripts/pomp_all_dept.R
nohup Rscript scripts/run_mif_all_dept.R $runlvl > output/all &

dept=$1
runlvl=$2
echo "Running for $dept at run level $runlvl :)"
Rscript scripts/pomp_cholera_haitiOCV.R $dept
nohup Rscript scripts/run_mif_haitiOCV.R  $dept $runlvl > out &
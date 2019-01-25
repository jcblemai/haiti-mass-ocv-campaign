runlvl=2

for dept in Centre Ouest Nord-Ouest Nord Sud  Nippes Nord-Est  Sud-Est Grande_Anse; do
    rm -R output/$dept/
    mkdir output/$dept;
    Rscript scripts/pomp_cholera_haitiOCV.R $dept
    nohup Rscript scripts/run_mif_haitiOCV.R  $dept $runlvl > output/$dept/out$dept &
done;


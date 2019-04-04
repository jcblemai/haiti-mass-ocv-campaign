#runlvl=$1
output_dir="output_16-04-init"

for dept in Artibonite; do #Centre Ouest Nord-Ouest Nord Sud  Nippes Nord-Est Sud-Est Grande_Anse; do
    rm -R output/$dept/
    mkdir output/$dept;
    Rscript scripts/pomp_cholera_haitiOCV.R $dept;
    cp -v "output/$dept/sirb_cholera_pomped_$dept.rda" "$output_dir/$dept/sirb_cholera_pomped_$dept.rda"
    #cp -v $output_dir/Grande_Anse/Haiti_OCV-Grande_Anse-param_logliks-10-l$runlvl.csv $output_dir/$dept/Haiti_OCV-$dept-param_logliks-10-l$runlvl.csv;
    # nohup Rscript scripts/run_tiny_mif_haitiOCV.R  $dept $runlvl > output/$dept/out$dept &
done;

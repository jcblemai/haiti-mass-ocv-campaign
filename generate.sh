runlvl=$1
output_dir="output_12-20-2gammaMOD"

for dept in Artibonite Centre Ouest Nord-Ouest Nord Sud  Nippes Nord-Est Sud-Est Grande_Anse; do
    mkdir output/$dept;
    Rscript scripts/pomp_cholera_haitiOCV.R $dept;
    cp -v "output/$dept/sirb_cholera_pomped_$dept.rda" "$output_dir/$dept/sirb_cholera_pomped_$dept.rda"
    #cp -v $output_dir/Grande_Anse/Haiti_OCV-Grande_Anse-param_logliks-10-l$runlvl.csv $output_dir/$dept/Haiti_OCV-$dept-param_logliks-10-l$runlvl.csv;
done;

runlvl=$1
output_dir="output_13-24-mob"

for dept in Artibonite Centre Ouest Nord-Ouest Nord Sud  Nippes Nord-Est  Sud-Est Grande_Anse; do
    mkdir output/$dept;
    Rscript scripts/pomp_cholera_haitiOCV.R $dept;
    cp -v "output/$dept/sirb_cholera_pomped_$dept.rda" "$output_dir/$dept/sirb_cholera_pomped_$dept.rda"
    #cp output/Artibonite/Haiti_OCV-Artibonite-param_logliks-10-l$runlvl.csv output/$dept/Haiti_OCV-###$dept-param_logliks-10-l$runlvl.csv;
done;

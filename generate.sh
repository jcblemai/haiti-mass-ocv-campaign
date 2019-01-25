runlvl=$1
for dept in Centre Ouest Nord-Ouest Nord Sud  Nippes Nord-Est  Sud-Est Grande_Anse; do
    mkdir output/$dept;
    Rscript scripts/pomp_cholera_haitiOCV.R $dept;
    cp output/Artibonite/Haiti_OCV-Artibonite-param_logliks-10-l$runlvl.csv output/$dept/Haiti_OCV-$dept-param_logliks-10-l$runlvl.csv;
done;

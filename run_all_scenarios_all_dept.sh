
mkdir output/
rm -R output/*
mkdir output/Simulations
Rscript scripts/pomp_all_dept.R;
nohup python3 all_scenarios_all_dept.py > out &

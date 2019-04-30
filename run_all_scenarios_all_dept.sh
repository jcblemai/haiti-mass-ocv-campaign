
mkdir output/
rm -R output/*
mkdir output/Simulations
Rscript scripts/pomp_all_dept;
nohup python3 all_scenarios_all_dept.py &

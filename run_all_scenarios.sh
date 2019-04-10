
mkdir output/
rm -R output/*
mkdir output/Simulations
sh generate.sh

nohup python3 all_scenarios.py &

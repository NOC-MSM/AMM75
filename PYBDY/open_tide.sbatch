#!/bin/bash 
#SBATCH --partition=short-serial 
#SBATCH --mem=20000
#SBATCH --job-name=amm75_tides
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1

start=`date +%s`
# executable
module add jaspy
micromamba activate /home/users/anwise/miniconda3/envs/pybdy
export JAVA_HOME=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.412.b08-1.el7_9.x86_64
pybdy -s namelist_AMM75_open_tide.bdy

end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))

echo "$1 Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
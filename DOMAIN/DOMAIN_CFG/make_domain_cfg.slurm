#!/bin/bash
#SBATCH --job-name=nemo_dom
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# Define memory required for this jobs. By default, you would 
# get just under 2 GB, but you can ask for up to 125 GB.
#SBATCH --mem=16G

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically 
#   using threading.
export OMP_NUM_THREADS=1

# Ensure the correct modules are loaded for ARCHER2
# Load modules listed in /work/n01/shared/nemo/setup
# Tested 10Jan22
# module swap craype-network-ofi craype-network-ucx
# module swap cray-mpich cray-mpich-ucx
# module load cray-hdf5-parallel/1.12.0.7
# module load cray-netcdf-hdf5parallel/4.7.4.7
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.2.1
module load cray-netcdf-hdf5parallel/4.9.0.1


# Launch the serial job
/work/n01/n01/anwise/NEMO/nemo_4.0.4/tools/DOMAINcfg/BLD/bin/make_domain_cfg.exe


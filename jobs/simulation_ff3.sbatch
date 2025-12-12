#!/bin/bash
#SBATCH --job-name=sim_ff3
#SBATCH --partition=compute
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --output=logs/sim_ff3_%j.out
#SBATCH --error=logs/sim_ff3_%j.err

# -------------------- environment --------------------
module purge
module load 2023
module load gcc/12.2.0

# R
export PATH="$HOME/R/4.4.1/bin:$PATH"
export R_LIBS_USER="$HOME/R/library/4.4"
unset R_HOME 2>/dev/null

# Gurobi
export GUROBI_HOME="$HOME/gurobi1203/linux64"
export GRB_LICENSE_FILE="$HOME/gurobi.lic"
export LD_LIBRARY_PATH="$GUROBI_HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# Threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export BLIS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Sequential MIQP → allow Gurobi to use all cores
export SPS_GUROBI_THREADS=$SLURM_CPUS_PER_TASK

# -------------------- run --------------------
cd $HOME/work/SparsePortfolioSelection
Rscript inst/simulations/hpc_mv_analysis_ff3.R

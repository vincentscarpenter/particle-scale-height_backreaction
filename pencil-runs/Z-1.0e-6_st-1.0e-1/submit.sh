#!/bin/bash
# Standard output and error:
#SBATCH -o ./Z-1.0e-6_st-1.0e-1.out.%j
#SBATCH -e ./Z-1.0e-6_st-1.0e-1.err.%j
#SBATCH -D ./
#SBATCH -J Z-1.0e-6_st-1.0e-1
#SBATCH --nodes=8
#SBATCH --tasks-per-node=40
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=all
#SBATCH --mail-user=carpenter@mpia.de
# Wall clock limit
#SBATCH --time=48:00:00

pc_run -f $PENCIL_HOME/config/hosts/isaac/isaac1.bc.rzg.mpg.de.conf

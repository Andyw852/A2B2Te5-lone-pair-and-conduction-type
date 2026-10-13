#!/bin/bash
#SBATCH --partition=cpu192
#SBATCH --job-name=Sn2Sb2Te5-ke-dft-cpu-S8.2_mass_soc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --qos=regular
cd $SLURM_SUBMIT_DIR
cp ../step3_uniform/CHGCAR CHGCAR
module purge
module load aocc/aocc-compiler-4.1.0 aocl/4.1.0 openmpi/4.1.5-aocc
ulimit -s unlimited
export OMP_NUM_THREADS=1
export OMPI_MCA_btl=^openib
export UCX_TLS=sm,self
mpirun -np $SLURM_NTASKS /public/home/wangchao/software/vasp.6.4.3/bin/vasp_ncl

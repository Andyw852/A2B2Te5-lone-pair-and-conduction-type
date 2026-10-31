#!/bin/bash
# SOC 机制检验 d1_nosoc_Sn2Bi2Te5
# 与 c2_pbe 同一 vasp 6.4.3；jzzn 默认 QOS 只有 1 分钟墙钟，必须 --qos=regular。
#SBATCH --partition=cpu192
#SBATCH --qos=regular
#SBATCH --job-name=1-0_SOC_d1_nosoc_SnBi
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --time=02:00:00

cd $SLURM_SUBMIT_DIR

# 该 VASP 以 AOCC/AOCL + OpenMPI 4.1.5-aocc 编译，运行库须一致
module purge
module load aocc/aocc-compiler-4.1.0 aocl/4.1.0 openmpi/4.1.5-aocc
ulimit -s unlimited
export OMP_NUM_THREADS=1
export OMPI_MCA_btl=^openib
export UCX_TLS=sm,self

rm -f OUTCAR OSZICAR
echo "[start] $(date) $(hostname)"
mpirun -np $SLURM_NTASKS /public/home/wangchao/software/vasp.6.4.3/bin/vasp_std > out.log 2>&1
echo "[done] $(date)"
grep -q "General timing" OUTCAR || { echo NO_GENERAL_TIMING; exit 1; }
grep -q "aborting loop because EDIFF is reached" OUTCAR || { echo NOT_CONVERGED; exit 1; }
echo OK

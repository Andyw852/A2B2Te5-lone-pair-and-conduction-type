#!/bin/bash
# SOC机制检验_HSE e1_hse_nosoc_Pb2Bi2Te5：PBE 无 SOC 自洽 → HSE06 无 SOC 自洽 → HSE06 无 SOC 的 c2 路径（204 点）
# 集群、模块、vasp_std 路径按 c2_hse 当时用的环境改（这里照抄 SOC机制检验/d1 的 jzzn 设置）。
#SBATCH --partition=cpu192
#SBATCH --qos=regular
#SBATCH --job-name=1-0_SOC_HSE_e1_PbBi
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --time=48:00:00
cd $SLURM_SUBMIT_DIR
module purge
module load aocc/aocc-compiler-4.1.0 aocl/4.1.0 openmpi/4.1.5-aocc
ulimit -s unlimited
export OMP_NUM_THREADS=1 OMPI_MCA_btl=^openib UCX_TLS=sm,self
VASP=/public/home/wangchao/software/vasp.6.4.3/bin/vasp_std
[ -f POTCAR ] || { echo "缺 POTCAR：从 c2_hse 的作业目录复制（元素顺序同 POSCAR 第 6 行）"; exit 1; }
step() {   # $1 = 步号，$2 = INCAR 文件
  cp $2 INCAR; rm -f OUTCAR OSZICAR
  echo "[start $1] $(date) $(hostname)"
  mpirun -np $SLURM_NTASKS $VASP > out.log.$1 2>&1
  grep -q "General timing" OUTCAR || { echo "NO_GENERAL_TIMING step $1"; exit 1; }
  grep -q "aborting loop because EDIFF is reached" OUTCAR || { echo "NOT_CONVERGED step $1"; exit 1; }
  cp OUTCAR OUTCAR.$1; cp OSZICAR OSZICAR.$1; [ -f vasprun.xml ] && cp vasprun.xml vasprun.xml.$1
}
rm -f KPOINTS_OPT
[ -f OUTCAR.1 ] || step 1 INCAR.1_pbe
[ -f OUTCAR.2 ] || step 2 INCAR.2_hse_scf
cp KPOINTS_OPT.path KPOINTS_OPT
step 3 INCAR.3_hse_path
gzip -f PROCAR_OPT && gzip -f vasprun.xml && rm -f WAVECAR vasprun.xml.1 vasprun.xml.2
echo OK

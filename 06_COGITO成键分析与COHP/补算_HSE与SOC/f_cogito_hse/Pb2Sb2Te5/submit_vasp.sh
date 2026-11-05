#!/bin/bash
# 补算 f（Pb2Sb2Te5）：PBE 不加 SOC 自洽 → HSE06 不加 SOC 自洽，两步各留一份 COGITO 所需文件（cogito_pbe/、cogito_hse/）。
#SBATCH --partition=cpu192
#SBATCH --qos=regular
#SBATCH --job-name=cohp_f_PbSb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48   # 144 条带可被 48 整除，HSE 步（无 NCORE）不会被 VASP 加到 160
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --time=72:00:00
cd $SLURM_SUBMIT_DIR
module purge
module load aocc/aocc-compiler-4.1.0 aocl/4.1.0 openmpi/4.1.5-aocc
ulimit -s unlimited
export OMP_NUM_THREADS=1 OMPI_MCA_btl=^openib UCX_TLS=sm,self
VSTD=/public/home/wangchao/software/vasp.6.4.3/bin/vasp_std
VNCL=/public/home/wangchao/software/vasp.6.4.3/bin/vasp_ncl
[ -f POTCAR ] || { echo "缺 POTCAR：从 <材料>/c2_hse 的作业目录复制（元素顺序同 POSCAR 第 6 行：B 位、Te、A 位）"; exit 1; }
step() {   # $1 = 步号，$2 = INCAR 文件，$3 = vasp 可执行文件
  cp $2 INCAR; rm -f OUTCAR OSZICAR
  echo "[start $1] $(date) $(hostname)"
  mpirun -np $SLURM_NTASKS $3 > out.log.$1 2>&1
  grep -q "General timing" OUTCAR || { echo "NO_GENERAL_TIMING step $1"; exit 1; }
  grep -q "aborting loop because EDIFF is reached" OUTCAR || { echo "NOT_CONVERGED step $1"; exit 1; }
  grep -q "NBANDS= *$(grep -o 'NBANDS *= *[0-9]*' $2 | grep -o '[0-9]*$')" OUTCAR || echo "WARN step $1: OUTCAR 中的 NBANDS 与 INCAR 不同（VASP 按核数调整过），记下实际值"
  cp OUTCAR OUTCAR.$1; cp OSZICAR OSZICAR.$1; cp IBZKPT IBZKPT.$1 2>/dev/null; cp INCAR INCAR.used.$1
}
keep() {   # 把 COGITO 读的文件复制到 $1/
  mkdir -p $1; cp POSCAR POTCAR KPOINTS INCAR OUTCAR vasprun.xml WAVECAR IBZKPT $1/
}
[ -f OUTCAR.1 ] || { step 1 INCAR.1_pbe $VSTD; keep cogito_pbe; }
[ -f OUTCAR.2 ] || { step 2 INCAR.2_hse $VSTD; keep cogito_hse; }
cmp -s IBZKPT.1 IBZKPT.2 || echo "WARN: 两步的 IBZKPT 不同（ISYM 2 → 3），记下来交回"
rm -f WAVECAR
echo OK

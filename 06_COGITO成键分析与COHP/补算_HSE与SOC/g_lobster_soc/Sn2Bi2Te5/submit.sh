#!/bin/bash
# 补算 g（Sn2Bi2Te5）：每步 VASP 之后在原地跑 LOBSTER，结果移到 lobster_<步>/。
# 第 3 步（HSE06+SOC）只在 RUN_HSE=1 且第 2 步的 LOBSTER 成功（即 LOBSTER 接受 SOC）时才跑：  sbatch --export=ALL,RUN_HSE=1 submit.sh
# 基组函数数：56（不加 SOC）、112（SOC 旋量）；NBANDS 取 64、128，都不少于它。
#SBATCH --partition=cpu192
#SBATCH --qos=regular
#SBATCH --job-name=cohp_g_SnBi
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --time=96:00:00
cd $SLURM_SUBMIT_DIR
module purge
module load aocc/aocc-compiler-4.1.0 aocl/4.1.0 openmpi/4.1.5-aocc
ulimit -s unlimited
export OMP_NUM_THREADS=1 OMPI_MCA_btl=^openib UCX_TLS=sm,self
VSTD=/public/home/wangchao/software/vasp.6.4.3/bin/vasp_std
VNCL=/public/home/wangchao/software/vasp.6.4.3/bin/vasp_ncl
[ -f POTCAR ] || { echo "缺 POTCAR：从 <材料>/c2_hse 的作业目录复制（元素顺序同 POSCAR 第 6 行：B 位、Te、A 位）"; exit 1; }
LOBSTER=${LOBSTER:-$(command -v lobster || ls -d /public/home/wangchao/software/*obster*/lobster* 2>/dev/null | head -1)}
[ -x "$LOBSTER" ] || { echo "找不到 LOBSTER 可执行文件：export LOBSTER=<路径> 后再交"; exit 1; }
echo "[lobster] $LOBSTER"
step() {   # $1 = 步号，$2 = INCAR 文件，$3 = vasp 可执行文件
  cp $2 INCAR; rm -f OUTCAR OSZICAR
  echo "[start $1] $(date) $(hostname)"
  mpirun -np $SLURM_NTASKS $3 > out.log.$1 2>&1
  grep -q "General timing" OUTCAR || { echo "NO_GENERAL_TIMING step $1"; exit 1; }
  grep -q "aborting loop because EDIFF is reached" OUTCAR || { echo "NOT_CONVERGED step $1"; exit 1; }
  grep -q "NBANDS= *$(grep -o 'NBANDS *= *[0-9]*' $2 | grep -o '[0-9]*$')" OUTCAR || echo "WARN step $1: OUTCAR 中的 NBANDS 与 INCAR 不同（VASP 按核数调整过），记下实际值"
  cp OUTCAR OUTCAR.$1; cp OSZICAR OSZICAR.$1; cp IBZKPT IBZKPT.$1 2>/dev/null; cp INCAR INCAR.used.$1
}
lob() {   # $1 = 结果目录名
  OMP_NUM_THREADS=$SLURM_NTASKS $LOBSTER > lobster.out 2>&1; rc=$?
  mkdir -p $1; mv lobster.out lobsterout *.lobster $1/ 2>/dev/null; cp lobsterin INCAR OUTCAR $1/
  gzip -f $1/COHPCAR.lobster $1/DOSCAR.lobster 2>/dev/null
  grep -qi "finished" $1/lobsterout && [ $rc -eq 0 ] && [ -f $1/ICOHPLIST.lobster ] || { echo "LOBSTER_FAILED $1"; tail -20 $1/lobster.out; return 1; }
  grep -i "spilling" $1/lobsterout
}
[ -f lobster_1_pbe_nosoc/ICOHPLIST.lobster ] || { step 1 INCAR.1_pbe_nosoc $VSTD; lob lobster_1_pbe_nosoc || exit 1; rm -f WAVECAR; }
[ -f lobster_2_pbe_soc/ICOHPLIST.lobster ]   || { step 2 INCAR.2_pbe_soc $VNCL;  lob lobster_2_pbe_soc || { echo "LOBSTER 不接受 SOC 的 WAVECAR：停在这里，不跑第 3 步，把 lobster_2_pbe_soc/ 交回"; exit 1; }; }
if [ "$RUN_HSE" = 1 ]; then
  [ -f WAVECAR ] || { echo "第 3 步要从第 2 步的 WAVECAR 续算，但 WAVECAR 已不在：删掉 lobster_2_pbe_soc/ 重跑第 2 步"; exit 1; }
  step 3 INCAR.3_hse_soc $VNCL; lob lobster_3_hse_soc || exit 1
  rm -f WAVECAR
else
  echo "未跑第 3 步：保留第 2 步的 WAVECAR，之后可直接 sbatch --export=ALL,RUN_HSE=1 submit.sh"
fi
echo OK

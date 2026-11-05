#!/bin/bash
# 补算 f（Pb2Sb2Te5）：在 cogito_pbe/、cogito_hse/ 里各跑一遍 COGITO 三件套（与 06_ 原来的 PBE 运行同一模板：06_COGITO成键分析与COHP/运行日志/<材料>/submit.sh）。
#SBATCH --partition=cpu192
#SBATCH --qos=regular
#SBATCH --job-name=cohp_f_cog_PbSb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --output=queue.out
#SBATCH --error=queue.err
#SBATCH --time=08:00:00
cd $SLURM_SUBMIT_DIR || exit 1
source /public/home/wangchao/miniconda3/etc/profile.d/conda.sh
conda activate atomate2_p_a
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "[cohp] COGITO=$(command -v COGITO) $(pip show cogito-dft 2>/dev/null | grep -i ^version)"
for d in cogito_pbe cogito_hse; do
  [ -f $d/WAVECAR ] || { echo "缺 $d/WAVECAR：先跑 submit_vasp.sh"; exit 1; }
  ( cd $d
    echo "=== $d [1/3] COGITO ===";        COGITO --dir ./ > cogito_run.log 2>&1 || { tail -30 cogito_run.log; exit 1; }
    echo "=== $d [2/3] COGITOanalyze ==="; COGITOanalyze --dir ./ > analyze.log 2>&1 || { tail -30 analyze.log; exit 1; }
    echo "=== $d [3/3] COGITOpost ===";    COGITOpost --dir ./ > post.log 2>&1 || { tail -30 post.log; exit 1; }
    grep -i "charge spilling\|ERROR\|Warning" cogito_run.log analyze.log post.log | head -20
    ls -l bond_info.txt ) || exit 1
done
echo "[cohp] all done $(date '+%F %T')"

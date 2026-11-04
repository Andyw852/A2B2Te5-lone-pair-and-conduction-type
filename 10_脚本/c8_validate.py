#!/usr/bin/env python3
"""HSE06+SOC 密网格 BoltzTraP2（c8_hse_dense，28×28×7/6 KPOINTS_OPT，分支 c8-hse-dense）的验收（2026-11-04）。只读。

对每种材料（结果不在本地时，直接从 origin/c8-hse-dense 用 git show 读）：
1. 插值倍数收敛：两个 mult 的 S，300/450/600 K，±3×10¹⁹–1×10²¹ cm⁻³，最大差。
2. 插值能带：bands_check.txt（主倍数）对直接计算的 HSE 路径 c2_hse/path_eigen.txt：
   路径价带顶之差；带边 0.15、0.3 eV 内与 0.3 eV 外的均方根与最大误差；Sn 系另对 c5 的真实价带谷（检验点 204、205）。
3. 与样品：霍尔浓度下 HSE 与 PBE（bt2_pbe_dense）的 S，与实验 S(T) 的均方根偏差（300–620 K，计算按 25 K 网格内插）。

用法：python 10_脚本/c8_validate.py [--branch origin/c8-hse-dense] [--mats Sn2Sb2Te5 Sn2Bi2Te5 …]
"""
import argparse, json, os, re, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import paper_figures as pf
import c2_valleys as cv

MAIN = 60                                    # 报告里用的主倍数（四种材料一致）


def read(rel, branch):
    p = os.path.join(ROOT, rel)
    if os.path.exists(p): return open(p, encoding='utf-8').read()
    r = subprocess.run(['git', '-C', ROOT, 'show', f'{branch}:{rel}'], capture_output=True)
    return r.stdout.decode('utf-8') if r.returncode == 0 else None


def mults(m, branch):
    base = f'22_能带基准_1-0/{m}/c8_hse_dense/gene'
    r = subprocess.run(['git', '-C', ROOT, 'ls-tree', '-z', '--name-only', f'{branch}:{base}'], capture_output=True)
    names = r.stdout.decode('utf-8').split('\0') if r.returncode == 0 else []
    if os.path.isdir(os.path.join(ROOT, base)): names += os.listdir(os.path.join(ROOT, base))
    return sorted({int(x[len('bt2_mult'):]) for x in names if re.fullmatch(r'bt2_mult\d+', x)})


def s_func(doc):
    T = list(doc['temperatures_K']); C = np.array(doc['carrier_conc_cm-3'], float)
    S = (np.array(doc['seebeck_xx_V_per_K'], float) + np.array(doc['seebeck_yy_V_per_K'], float)) / 2 * 1e6

    def f(N, T0):
        i = T.index(T0); sg = np.sign(N); ok = sg * C[i] > 0; x = np.log10(sg * C[i][ok]); o = np.argsort(x)
        return float(np.interp(np.log10(abs(N)), x[o], S[i][ok][o]))
    return np.array(T, float), f


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--branch', default='origin/c8-hse-dense'); ap.add_argument('--mats', nargs='+', default=list(pf.MATS))
    a = ap.parse_args()
    for m in a.mats:
        base = f'22_能带基准_1-0/{m}/c8_hse_dense/gene'; ms = mults(m, a.branch)
        print(f'\n== {m}  倍数 {ms}')
        if not ms: print('   还没有结果'); continue
        src = json.loads(read(f'{base}/source.json', a.branch))
        print(f"   {src['grid']}；网格带隙 {src['gap_eV']} eV；窗口内 Kramers 最大差 {1000 * src['kramers_max_split_eV']:.1f} meV")
        docs = {k: json.loads(read(f'{base}/bt2_mult{k}/boltztrap_crta.json', a.branch)) for k in ms}
        for k, d in docs.items():
            if d['functional'] != 'HSE06+SOC' or d['kmesh_mult'] != k or 'no scissor' not in d['protocol']:
                print(f"   [!] mult {k}：functional={d['functional']} kmesh_mult={d['kmesh_mult']} protocol={d['protocol']}")
        main_ = MAIN if MAIN in docs else max(ms)
        # 1 倍数收敛
        if len(ms) > 1:
            fs = {k: s_func(d)[1] for k, d in docs.items()}; worst = (0, None)
            for T0 in (300, 450, 600):
                for N in (3e19, 1e20, 3e20, 1e21, -3e19, -1e20, -1e21):
                    v = [fs[k](N, T0) for k in ms]; dv = max(v) - min(v)
                    if dv > worst[0]: worst = (dv, (T0, N))
            print(f'   1. 倍数收敛：{ms} 之间 S 最多差 {worst[0]:.2f} µV/K（{worst[1][0]} K，{worst[1][1]:+.0e} cm⁻³）'
                  + ('  ✓' if worst[0] < 1 else '  ✗ 未收敛'))
        # 2 插值能带
        txt = read(f'{base}/bt2_mult{main_}/bands_check.txt', a.branch)
        nel = pf.NEL[m]; Ep = cv.read_path_eigen(os.path.join(pf.rf.base(m), 'c2_hse', 'path_eigen.txt'))
        if txt:
            B = np.array([[float(x) for x in l.split()] for l in txt.splitlines() if l.strip() and not l.startswith('#')])
            get = lambda i, b: B[(B[:, 0] == i) & (B[:, 1] == b), 2][0]
            vbm, cbm = np.nanmax(Ep[:, nel - 1]), np.nanmin(Ep[:, nel])
            print(f'   2. 插值能带（mult {main_}）：路径价带顶 插值 − 直接 = {1000 * (max(get(i, nel) for i in range(len(Ep))) - vbm):+.1f} meV')
            for lo, hi, lab in ((0, 0.15, '带边 0.15 eV 内'), (0, 0.3, '带边 0.3 eV 内'), (0.3, 9, '带边 0.3 eV 外')):
                d = np.array([1000 * (get(i, b) - Ep[i, b - 1]) for i in range(len(Ep)) for b in range(nel - 3, nel + 5)
                              if lo <= ((vbm - Ep[i, b - 1]) if b <= nel else (Ep[i, b - 1] - cbm)) < hi])
                print(f'      {lab}：{len(d)} 个，均方根 {np.sqrt((d ** 2).mean()):.1f}、最大 {abs(d).max():.1f} meV')
            kc = read(f'22_能带基准_1-0/{m}/bt2_check_kpoints.txt', a.branch) or ''
            for i, e in re.findall(r'(\d+) = c5 v\d（HSE 直接计算 ([\d.]+) eV）', kc.splitlines()[0] if kc else ''):
                print(f'      c5 真实价带谷（检验点 {i}）：插值 − 直接 = {1000 * (get(int(i), nel) - float(e)):+.1f} meV')
        # 3 与样品
        T, f = s_func(docs[main_]); H = pf.HALL[m]['N']; Tp, fp = pf.pbe_dense(m)
        vh = np.array([f(H, int(t)) for t in T]); vp = np.array([fp(H, i) for i in range(len(Tp))])
        dh = [np.interp(t, T, vh) - s for t, s in pf.EXP[m]]; dp = [np.interp(t, Tp, vp) - s for t, s in pf.EXP[m]]
        print(f"   3. 霍尔浓度 {H:+.2e} cm⁻³{'（仪器标红）' if pf.HALL[m]['flagged'] else ''}：")
        print(f"      约 300 K：实验 {pf.EXP[m][0][1]:.1f}，HSE {np.interp(pf.EXP[m][0][0], T, vh):.1f}，PBE {np.interp(pf.EXP[m][0][0], Tp, vp):.1f} µV/K")
        print(f"      约 615 K：实验 {pf.EXP[m][-1][1]:.1f}，HSE {np.interp(pf.EXP[m][-1][0], T, vh):.1f}，PBE {np.interp(pf.EXP[m][-1][0], Tp, vp):.1f} µV/K")
        print(f'      300–620 K 均方根偏差：HSE {np.sqrt(np.mean(np.square(dh))):.1f}，PBE {np.sqrt(np.mean(np.square(dp))):.1f} µV/K')
        print(f'      3×10¹⁹ cm⁻³、300 K：p 型 HSE {f(3e19, 300):.1f} / PBE {fp(3e19, 0):.1f}；n 型 HSE {f(-3e19, 300):.1f} / PBE {fp(-3e19, 0):.1f} µV/K')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Generate UNRELAXED hypothetical Pb2B2Te5 rocksalt models and defect proposals.

No energies, DFT calculations, or annealing simulations are performed. 6.2 A is
an explicit synthetic starting lattice parameter, not a measured/fitted value.
The output POSCAR files do not define a DFT electronic charge state by themselves.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np


VERSION = "1.0"
GHOST = "Vac"
SEEDS = (1729, 2718, 31415, 161803)


def dump_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def fingerprint(lattice, sites, occupancies):
    payload = {
        "lattice": np.round(lattice, 12).tolist(),
        "frac": [s["fractional_coordinates"] for s in sites],
        "occupancy": occupancies,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def lattice_sites(a, repeat):
    # Rows are FCC translations in units a/2. Determinant = 10, so the cell
    # contains 5 FCC cation sites and 5 FCC anion sites before vacancies.
    matrix = np.array([[1, -1, 0], [0, 1, -1], [4, 3, 3]], dtype=float)
    base = 0.5 * a * matrix
    lattice = np.diag(repeat) @ base
    inverse = np.linalg.inv(matrix)
    sites = []
    for sublattice in ("cation", "anion"):
        offset = int(sublattice == "anion")
        for q in range(5):
            # q = (x+y+z)/a for cations. Anions have q+1/2.
            frac0 = np.mod(np.array([2 * q + offset, 0, 0]) @ inverse, 1.0)
            frac0[np.isclose(frac0, 1.0, atol=1e-12)] = 0.0
            for ijk in itertools.product(*(range(n) for n in repeat)):
                frac = (frac0 + np.array(ijk)) / np.array(repeat)
                sites.append({
                    "site_id": f"{'C' if sublattice == 'cation' else 'T'}_q{2*q+offset:02d}_"
                               + "_".join(map(str, ijk)),
                    "sublattice": sublattice,
                    "q_times_2_mod_10": 2 * q + offset,
                    "base_plane_index": q,
                    "replica": list(ijk),
                    "fractional_coordinates": np.round(frac, 14).tolist(),
                    "cartesian_coordinates_A": (frac @ lattice).tolist(),
                })
    multiplicity = int(np.prod(repeat))
    assert len(sites) == 10 * multiplicity
    assert np.isclose(np.linalg.det(lattice), 1.25 * a ** 3 * multiplicity)
    # Check stacking labels from Cartesian coordinates, modulo the 5-plane cell.
    for site in sites:
        q2 = 2 * sum(site["cartesian_coordinates_A"]) / a
        assert abs((q2 - site["q_times_2_mod_10"] + 5) % 10 - 5) < 1e-9
    return lattice, sites


def periodic_distances(lattice, sites):
    coords = np.array([s["fractional_coordinates"] for s in sites])
    differences = coords[:, None, :] - coords[None, :, :]
    differences -= np.rint(differences)
    dist2 = np.full(differences.shape[:2], np.inf)
    # All 27 adjacent images are sufficient for this fixed, modestly skewed basis.
    for shift in itertools.product((-1, 0, 1), repeat=3):
        cart = (differences + np.array(shift)) @ lattice
        dist2 = np.minimum(dist2, np.sum(cart * cart, axis=-1))
    np.fill_diagonal(dist2, np.inf)
    return np.sqrt(dist2)


def ordered_occupancies(sites, b):
    sequence = [b, "Pb", "Pb", b, GHOST]
    return ["Te" if s["sublattice"] == "anion" else sequence[s["base_plane_index"]]
            for s in sites]


def make_parent_occupancies(sites, b, family, seed):
    result = ordered_occupancies(sites, b)
    if family == "O":
        return result
    rng = np.random.default_rng(seed)
    chosen = [i for i, s in enumerate(sites)
              if s["sublattice"] == "cation" and (family == "D" or result[i] != GHOST)]
    values = [result[i] for i in chosen]
    rng.shuffle(values)
    for i, val in zip(chosen, values):
        result[i] = val
    return result


def validate_parent(sites, occupancy, b, family, multiplicity):
    assert Counter(occupancy) == Counter({"Pb": 2 * multiplicity, b: 2 * multiplicity,
                                         "Te": 5 * multiplicity, GHOST: multiplicity})
    assert all(occupancy[i] == "Te" for i, s in enumerate(sites) if s["sublattice"] == "anion")
    if family in ("O", "V"):
        assert all((occupancy[i] == GHOST) == (s["base_plane_index"] == 4)
                   for i, s in enumerate(sites) if s["sublattice"] == "cation")


def write_poscar(path, title, lattice, sites, occupancy, b):
    species = ["Pb", b, "Te"]
    groups = [[i for i, symbol in enumerate(occupancy) if symbol == element] for element in species]
    assert all(groups)
    indices = sum(groups, [])
    lines = [title, "1.0"]
    lines += ["  " + " ".join(f"{v:.14f}" for v in vec) for vec in lattice]
    lines += [" ".join(species), " ".join(str(len(group)) for group in groups), "Direct"]
    lines += ["  " + " ".join(f"{x:.14f}" for x in sites[i]["fractional_coordinates"]) for i in indices]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return {i: j + 1 for j, i in enumerate(indices)}


def ordered_class(site, occupancy):
    if occupancy != "Te":
        return occupancy
    q2 = site["q_times_2_mod_10"]
    # Inversion q -> 3-q maps q2 -> 6-q2 (mod 10).
    if q2 in (1, 5):
        return "Te_between_B_Pb"
    if q2 == 3:
        return "Te_between_Pb_Pb"
    if q2 in (7, 9):
        return "Te_between_B_Vac"
    raise AssertionError(q2)


def defect_definitions(b):
    # Occupied native cation vacancies are distinct from the structural empty
    # cation sublattice sites, denoted *_on_Vac (not a generic interstitial).
    return [
        (f"{b}_Pb", "Pb", b, "cation", True),
        (f"Pb_{b}", b, "Pb", "cation", True),
        ("Pb_on_Vac", GHOST, "Pb", "cation", True),
        (f"{b}_on_Vac", GHOST, b, "cation", True),
        ("V_Te", "Te", GHOST, "anion", True),
        ("V_Pb", "Pb", GHOST, "cation", False),
        (f"V_{b}", b, GHOST, "cation", False),
        ("Te_on_Vac", GHOST, "Te", "cation", False),
        ("Pb_Te", "Te", "Pb", "anion", False),
        (f"{b}_Te", "Te", b, "anion", False),
        ("Te_Pb", "Pb", "Te", "cation", False),
        (f"Te_{b}", b, "Te", "cation", False),
    ]


def choose_sites(sites, occupancy, source, sublattice, family, seed, focus, key):
    eligible = [i for i, s in enumerate(sites) if occupancy[i] == source and s["sublattice"] == sublattice]
    if family == "O":
        representatives = {}
        for i in eligible:
            representatives.setdefault(ordered_class(sites[i], source), i)
        return [(i, group) for group, i in sorted(representatives.items())]
    # Deterministic sample, not an assertion that random-site environments are
    # equivalent. Full site screening is still needed if ranking is close.
    local_seed = int(hashlib.sha256(f"{seed}:{key}".encode()).hexdigest()[:16], 16)
    rng = np.random.default_rng(local_seed)
    count = min(len(eligible), 4 if focus else 1)
    selected = sorted(rng.choice(eligible, size=count, replace=False).tolist())
    return [(i, "sampled_site_not_symmetry_reduced") for i in selected]


def geometry_summary(occupancy, distances):
    included = [i for i, val in enumerate(occupancy) if val != GHOST]
    minimum = float(distances[np.ix_(included, included)].min())
    assert minimum > 1.0, "overlap or unit error"
    return {"atom_count": len(included),
            "species_counts": dict(sorted(Counter(occupancy[i] for i in included).items())),
            "ghost_vacancy_count": occupancy.count(GHOST),
            "minimum_periodic_pair_distance_A": minimum}


def readme_text(material, args, counts):
    return f"""理论初始模型与缺陷候选（全部未弛豫，未开展 DFT）

材料：{material}。初始立方岩盐母格常数 a={args.a:g} Å；这是人为起始值，不是实验值、
不是 DFT 优化值。共同晶胞重复数：{args.repeat}。所有坐标来自同一个 FCC 岩盐母格，
不代表退火前后真实结构，也不含已经计算出来的“平整化”。

基础晶胞三个行向量为 a/2 × (1,-1,0)、(0,1,-1)、(4,3,3)，有 5 个阳离子母格位点
及 5 个 Te 位点。保持 Pb:B:Te=2:2:5 时，每 5 个阳离子母格位点中有 1 个结构空位。
默认 4×4×1 为 144 原子、16 个 ghost vacancy，完整母格共 160 位点。

O：阳离子 (111) 层 q mod 5 为 [B,Pb,Pb,B,Vac]，一个有序假设模型。
V：结构空位固定在 q=4 层，其余阳离子 Pb/B 随机混排，4 个独立种子。
D：Pb/B/结构空位在阳离子子晶格完全随机混排，4 个独立种子。
Te 子晶格先保持完整。D 是有限随机构型集合，不是 SQS，也不是平衡无序系综。
O/D/V 是三个静态占位假设，不能未经验证直接称为“退火后/退火前”状态。

父结构数：{counts['parents']}；缺陷候选数：{counts['defects']}；作业清单行数：{counts['jobs']}。
每个 D/V 父胞五类重点缺陷各抽取最多 4 个位点，另外七类各取 1 个位点。
重点缺陷：B_Pb、Pb_B、Pb_on_Vac、B_on_Vac、V_Te。
补充缺陷：V_Pb、V_B、Te_on_Vac、Pb_Te、B_Te、Te_Pb、Te_B。
有序 O 的同层平移及 q→3−q 反演等价位点已合并；Pb、B、Vac 各一类，Te 三类。
这里只利用理想母胞中显式成立的等价关系，未假设 D/V 中抽样位点等价。

这十二类覆盖该岩盐位点模型的文献缺陷集合，但并非穷尽一切间隙、复合缺陷。
*_on_Vac 表示填充原有结构空位，不是一般间隙原子；V_Pb/V_B 是新增的原子空位。
每个候选只改变一个母格位点，atom_count 随增/减原子变化；并非每个缺陷仍为2:2:5。
POSCAR 本身没有设置电荷态。defect.json / jobs.csv 的 intended_charge_state=0 仅表示这一批的
中性预筛用途；不能由中性形成能判断施主/受主电离和 P/N。带电态与修正需另行计算。

建议执行顺序：
1. 收敛同一赝势、截断能、k 点密度与 SOC 策略，在合理晶格常数附近优化父结构。
2. 记录放松前后占位有序度、空位聚集、层起伏、局部键长及能量；检查是否留在目标相。
3. 从已弛豫父结构重新生成/映射单点缺陷，比较多位点/多构型能量分布。
   本目录缺陷 POSCAR 是未弛豫几何草案，不能将它们与后续已弛豫父胞直接混用算形成能。
   初始生成脚本只生成理想 FCC 初始坐标；另附 regenerate_defects_from_relaxed_parent.py
   可从保持原子身份/顺序的父胞 CONTCAR 再生缺陷。需独立确认弛豫收敛，不能仅凭文件名。
   弛豫后须在实际计算流程中按 site_map 保持原子身份、核对占位是否迁移，再重新施加缺陷；
   对填充空位的候选，应结合弛豫后局部邻居重新定位空腔，不能直接沿用旧空位坐标。
4. 中性十二类粗筛后，对有竞争力的缺陷扩大位点采样，再计算多个电荷态与有限尺寸修正。
5. 使用明确的共同外部化学势库，并说明各分支的亚稳/部分平衡假设，计算缺陷形成能
   及热平衡/冻结缺陷电荷中性解；不得直接照搬有序稳定相的化学势顶点用于亚稳无序相。
   电荷中性解给出 n、p；Seebeck 与 Hall 符号还取决于两类载流子的迁移率及局域化。
   多数载流子数量变化不能直接等同实验 S/Hall 符号翻转，也不保证发生任何 P↔N 转变。

文件：
- parents/*/POSCAR：未弛豫父结构。
- parents/*/site_map.json：包含 ghost vacancy 的完整母格索引及 POSCAR 原子映射。
- defects/*/*/POSCAR：单点缺陷几何候选。
- defects/*/*/defect.json：发生改变的母格位点、组成变化及父结构指纹。
- jobs.csv：全部状态为 NOT_RUN；不包含任何计算成功或能量数据。
- manifest.json：几何、种子、数量、限制及指纹。
- validation.json：整数化学计量、(111) 层标签、重复结构、最近距离、无重叠检查。

不提供 POTCAR、INCAR 或运行提交脚本，不会自动运行或提交 VASP。
所有文件由 scripts/build_rocksalt_models.py 可重现；再次运行请指定新的 --out 目录。

在“补算实施方案_20261007”目录中的命令示例：
python scripts/build_rocksalt_models.py --material Pb2Sb2Te5 --a 6.2 --repeat 4 4 1 --out new_PbSb_144
替换 --material Pb2Bi2Te5 可生成 Bi 对照；相同 a 仍只是人工猜测，不代表两者平衡晶格相同。
脚本依赖 Python 与 numpy。--out 只允许不存在或空目录，遇到非空目录会报错保护已有文件。
脚本不会覆盖当前模型，也不会自动读取、弛豫、转换或提交计算结果。

弛豫后映射工具示例（请将 PATH_TO_CONTCAR 替换为实际且已核验的父胞弛豫文件）：
python scripts/regenerate_defects_from_relaxed_parent.py --site-map rocksalt_Pb2Sb2Te5_144atoms_UNRELAXED/parents/O_ordered/site_map.json --relaxed-poscar PATH_TO_CONTCAR --templates rocksalt_Pb2Sb2Te5_144atoms_UNRELAXED/defects/O_ordered --ghost-method neighbor --out mapped_O
工具检查元素计数和元素块顺序，但不能从 CONTCAR 单独证明同种元素原子没有互换身份。
默认非仿射周期位移超过 1.5 Å 或原子间距小于 1.5 Å 即拒绝映射；两个阈值可显式调整，
但调整前必须检查重构和原子身份。--check-only 只校验，不写候选。
填充空位必须显式选择 --ghost-method：affine 保留原空位分数坐标；neighbor 用原六个
最近 Te 邻居的周期非仿射位移平均估计位置。两者都是缺陷弛豫起始猜测，不是真实空位测量。
输出记录来源文件 SHA256、位移诊断和空位定位假设，状态仍为 NOT_RUN；非空输出目录拒绝覆盖。
映射工具不继承 Selective dynamics 约束，不设置电荷或生成 DFT 输入参数。
目前只通过恒等映射、仿射晶格变化及周期边界位移的软件验证，没有实际弛豫结果可供映射。
"""


def build(args):
    if args.a <= 0 or any(n < 1 for n in args.repeat):
        raise ValueError("a and repeat must be positive")
    out = args.out.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Refusing to overwrite nonempty output: {out}")
    out.mkdir(parents=True, exist_ok=True)
    b = "Sb" if args.material == "Pb2Sb2Te5" else "Bi"
    lattice, sites = lattice_sites(args.a, args.repeat)
    distances = periodic_distances(lattice, sites)
    assert np.isclose(distances.min(), args.a / 2, atol=1e-9)
    multiplicity = int(np.prod(args.repeat))
    parents, candidates, jobs, seen = [], [], [], set()
    setups = [("O", None)] + [(family, seed) for family in ("V", "D") for seed in SEEDS]
    for family, seed in setups:
        parent_id = "O_ordered" if family == "O" else f"{family}_seed{seed}"
        occupancy = make_parent_occupancies(sites, b, family, seed)
        validate_parent(sites, occupancy, b, family, multiplicity)
        parent_dir = out / "parents" / parent_id
        parent_dir.mkdir(parents=True)
        digest = fingerprint(lattice, sites, occupancy)
        if digest in seen:
            raise AssertionError("Repeated parent occupancy; enlarge small supercell")
        seen.add(digest)
        geometry = geometry_summary(occupancy, distances)
        atom_mapping = write_poscar(parent_dir / "POSCAR", f"UNRELAXED synthetic a={args.a:g} {parent_id}",
                                    lattice, sites, occupancy, b)
        mapped = [{**s, "occupancy": occupancy[i],
                   "poscar_atom_index_1based": atom_mapping.get(i),
                   "is_structural_ghost_vacancy": occupancy[i] == GHOST}
                  for i, s in enumerate(sites)]
        dump_json(parent_dir / "site_map.json", {
            "status": "UNRELAXED_NOT_RUN", "parent_id": parent_id,
            "lattice_vectors_A_rows": lattice.tolist(), "sites": mapped,
        })
        metadata = {"parent_id": parent_id, "family": family, "random_seed": seed,
                    "status": "UNRELAXED_NOT_RUN", "fingerprint_sha256": digest,
                    "poscar": str((parent_dir / "POSCAR").relative_to(out)), **geometry}
        dump_json(parent_dir / "parent.json", metadata)
        parents.append(metadata)
        jobs.append({"job_id": parent_id, "stage": "PARENT_RELAX_REQUIRED", "parent_id": parent_id,
                     "defect_type": "", "site_id": "", "intended_charge_state": "",
                     "atom_count": geometry["atom_count"], "status": "NOT_RUN",
                     "poscar": metadata["poscar"]})
        for name, source, target, sublattice, focus in defect_definitions(b):
            selected = choose_sites(sites, occupancy, source, sublattice, family, seed, focus, name)
            for nth, (site_index, site_class) in enumerate(selected, 1):
                defect_id = f"{name}__s{nth:02d}"
                defect_dir = out / "defects" / parent_id / defect_id
                defect_dir.mkdir(parents=True)
                changed = occupancy.copy()
                changed[site_index] = target
                assert sum(a != c for a, c in zip(occupancy, changed)) == 1
                expected = Counter(occupancy)
                expected[source] -= 1
                expected[target] += 1
                assert Counter(changed) == +expected
                write_poscar(defect_dir / "POSCAR", f"UNRELAXED {parent_id} {defect_id} neutral-preview",
                             lattice, sites, changed, b)
                diff = {element: changed.count(element) - occupancy.count(element)
                        for element in ("Pb", b, "Te")}
                info = {"defect_id": defect_id, "defect_type": name, "parent_id": parent_id,
                        "parent_fingerprint_sha256": digest, "status": "UNRELAXED_NOT_RUN",
                        "intended_charge_state": 0,
                        "charge_note": "Geometry only. POSCAR does not set electronic charge. Charged screening remains required.",
                        "source_occupancy": source, "target_occupancy": target,
                        "composition_change_vs_parent": diff,
                        "changed_site": sites[site_index], "site_class": site_class,
                        "sampling_priority": "mechanism_focused" if focus else "remaining_twelve_type_screen",
                        "fingerprint_sha256": fingerprint(lattice, sites, changed),
                        "poscar": str((defect_dir / "POSCAR").relative_to(out)),
                        **geometry_summary(changed, distances)}
                dump_json(defect_dir / "defect.json", info)
                candidates.append(info)
                jobs.append({"job_id": parent_id + "__" + defect_id,
                             "stage": "REGENERATE_FROM_RELAXED_PARENT_THEN_NEUTRAL_SCREEN",
                             "parent_id": parent_id, "defect_type": name,
                             "site_id": sites[site_index]["site_id"], "intended_charge_state": 0,
                             "atom_count": info["atom_count"], "status": "NOT_RUN", "poscar": info["poscar"]})
    counts = {"parents": len(parents), "defects": len(candidates), "jobs": len(jobs)}
    manifest = {"generator_version": VERSION, "material": args.material,
                "status": "ALL_UNRELAXED_NO_DFT_PERFORMED", "initial_cubic_a_A": args.a,
                "a_provenance": "synthetic user-independent starting guess; not experimental or optimized",
                "repeat": args.repeat, "lattice_vectors_A_rows": lattice.tolist(),
                "base_cation_plane_sequence_q_mod_5": [b, "Pb", "Pb", b, GHOST],
                "initial_lattice_parameter_is_same_for_all_families": True,
                "random_seeds_per_disordered_family": list(SEEDS),
                "lattice_site_count_including_ghosts": len(sites),
                "counts": counts,
                "defect_counts_by_parent": dict(Counter(c["parent_id"] for c in candidates)),
                "defect_counts_by_type": dict(Counter(c["defect_type"] for c in candidates)),
                "ordered_symmetry_reduction": "in-plane translations and inversion q -> 3-q; 1 Pb,1 B,1 Vac,3 Te classes",
                "limitations": ["O,V,D denote hypothetical occupational order; no assignment to measured annealing states",
                                "No atomic distortions, cell relaxation, energies, or finite-temperature entropy computed",
                                "D random cells are not optimized SQS structures",
                                "Finite-size and random-configuration convergence remain required",
                                "Twelve site-defect types do not exhaust interstitials or complexes",
                                "Candidate POSCARs must be regenerated from consistently relaxed parents before energy comparisons",
                                "Neutral geometry proposals do not establish charge states or p/n conduction"],
                "parents": parents}
    dump_json(out / "manifest.json", manifest)
    with (out / "jobs.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(jobs[0]))
        writer.writeheader()
        writer.writerows(jobs)
    # Verify the actual serialized POSCAR counts and all row/job paths.
    checked = 0
    for job in jobs:
        lines = (out / job["poscar"]).read_text(encoding="ascii").splitlines()
        n = sum(map(int, lines[6].split()))
        assert n == job["atom_count"] == len(lines) - 8
        xyz = np.array([[float(v) for v in line.split()] for line in lines[8:]])
        assert np.all(xyz >= -1e-12) and np.all(xyz < 1.0 + 1e-12)
        assert len(np.unique(np.round(xyz, 10), axis=0)) == n
        checked += 1
    validation = {"result": "PASS_GEOMETRY_AND_BOOKKEEPING_ONLY", "df_t_calculations_performed": 0,
                  "checks": {"parent_exact_stoichiometry": True, "fcc_cell_volume_and_plane_labels": True,
                             "unique_parent_fingerprints": True, "ordered_vacancy_plane_in_O_and_V": True,
                             "one_site_change_per_candidate": True, "candidate_composition_deltas": True,
                             "serialized_poscar_atom_counts": True, "no_coincident_sites": True,
                             "all_jobs_NOT_RUN": all(j["status"] == "NOT_RUN" for j in jobs)},
                  "serialized_poscar_files_checked": checked,
                  "minimum_distance_full_parent_lattice_A": float(distances.min()),
                  "expected_rocksalt_nearest_neighbor_A": args.a / 2,
                  "minimum_occupied_distance_A": min(c["minimum_periodic_pair_distance_A"] for c in candidates),
                  "note": "Checks only verify constructed geometries, not physical stability or electronic properties."}
    dump_json(out / "validation.json", validation)
    (out / "README_模型说明.txt").write_text(readme_text(args.material, args, counts), encoding="utf-8")
    print(json.dumps({"output": str(out), **counts, "validation": validation["result"]}, ensure_ascii=False))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--material", choices=("Pb2Sb2Te5", "Pb2Bi2Te5"), default="Pb2Sb2Te5")
    parser.add_argument("--a", type=float, default=6.2, help="synthetic initial conventional cubic parameter in A")
    parser.add_argument("--repeat", type=int, nargs=3, default=[4, 4, 1], metavar=("N1", "N2", "N3"))
    parser.add_argument("--out", type=Path, required=True, help="new or empty output directory")
    return parser.parse_args()


if __name__ == "__main__":
    build(parse_args())

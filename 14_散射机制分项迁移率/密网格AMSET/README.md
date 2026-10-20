# 密网格 AMSET 结果入库约定

本目录存放 hanhai25 上 128 核全量 AMSET（`interpolation_factor=1`，`energy_cutoff` 用默认 1.5）的结果。

- `14_/散射机制分项迁移率/` 的 `scattering_mobility.py` 与
  `18_ZT参数与孤对电子/` 的 `zt_parameters.py` 都只读本目录下的 `transport.json`；
  目录缺失或为空时，两个脚本会跳过该材料并打印提示。
- 12x12x2 粗网格结果已作废（低估 3–4 倍），不要再放回来。

## 每个材料目录需要的文件

| 文件 | 来源 | 用途 |
|---|---|---|
| `transport.json` | `<材料>/ke-dft-cpu/step8_amset_full_f1/transport_<mesh>.json` | 迁移率、电导、Seebeck、κ_e |
| `settings.yaml` | 同目录 | 溯源（温度/浓度网格、interpolation_factor、bandgap） |
| `meta.txt` | 手工写 | 作业号、墙钟、峰值内存、cluster |

## 入库步骤（结果跑完后）

```bash
# hanhai25 上（复用 ControlMaster）
ssh -S <socket> hanhai25 'cat ~/work/ke_soc_dense_20261019/<材料>/ke-dft-cpu/step8_amset_full_f1/transport_*.json' \
  > 14_散射机制分项迁移率/密网格AMSET/<材料>/transport.json
```

然后在仓库根目录运行：

```bash
python 10_脚本/scattering_mobility.py .      # 14_ 的迁移率表与图
python 10_脚本/zt_parameters.py .            # 18_ 的 ZT 参数与 PF
```

## 当前状态（2026-10-20 深夜）

| 材料 | 作业 | 状态 |
|---|---|---|
| Pb2Sb2Te5 | 296539 | 全量在跑（inelastic 82%+） |
| Sn2Sb2Te5 | 296540 | 全量在跑 |
| Pb2Bi2Te5 | 297159 | 精简版（6 浓度 x 5 温度，700G）在跑 |
| Sn2Bi2Te5 | 297160 | 精简版在跑 |

注意：Bi 系的精简版是 6 浓度 x 5 温度，与 Sb 系的 13 浓度 x 17 温度**口径不同**，
入库时要在 `meta.txt` 里写明，图与表中不能直接逐点对比。
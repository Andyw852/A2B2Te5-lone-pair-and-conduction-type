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

## 当前状态（2026-10-22）

| 材料 | 已入库 | 作业 | 说明 |
|---|---|---|---|
| Pb2Sb2Te5 | `transport_300K.json` | 297926 | 300 K × 6 浓度；400–700 K 未算 |
| Sn2Sb2Te5 | `transport_300K.json` | 297927 | 同上 |
| Pb2Bi2Te5 | — | 297928（300 K，在跑） | 300–700 K 的精简版 297159 在 24 h 墙钟内超时，没有结果 |
| Sn2Bi2Te5 | `transport_300-700K.json` | 297160 | 300–700 K × 6 浓度 |

- 四个材料用同一组浓度：−1e20、−1e19、−1e18、1e18、1e19、1e20，可以逐点对比。
- 同一材料的多个 `transport*.json` 会按温度合并，**温度不能重复**。Sn2Bi2Te5 的 300 K 已经在 `transport_300-700K.json` 里，单独算 300 K 的作业 297929 结果不要放进这个目录（否则合并时报错）；它只能用来核对可复现性。
- 实测耗时（128 进程）：Sb 系 300 K × 6 浓度 3.3–3.9 h；SnBi 5 个温度 × 6 浓度 23.1 h；PbBi 5 个温度超过 24 h。

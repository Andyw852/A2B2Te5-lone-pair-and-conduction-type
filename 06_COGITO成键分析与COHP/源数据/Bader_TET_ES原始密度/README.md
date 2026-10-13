# Bader_TET/ES 原始密度

来源：工作目录 `band_and_elf/Bader_TET/<材料>/ES/`。四个材料的
`CONTCAR` 与当前 COGITO POSCAR 的晶胞体积、元素和分数坐标逐项一致。
这些是 PBE-D3 静态计算的原始电荷/ELF 输入（12×12×2，`CHGCAR`、
`AECCAR0/1/2`、`ELFCAR` 及设置和 OUTCAR），可用标准 Bader 工具重新生成
`ACF.dat`。

原目录没有 `ACF.dat`，也没有保存旧 TET/ET 的定义、参考价电子或归一化公式。
因此本报告将此目录作为 Bader 后处理的可复算源，不把尚未定义的 Bader 电荷
直接命名为 TET，也不用于恢复旧 ES=2×ICOBI 表。旧 ES 所需的 LOBSTER
`ICOBILIST` 仍缺失；唯一找到的 Pb2Sb2Te5 LOBSTER 输出因投影失败而无效。

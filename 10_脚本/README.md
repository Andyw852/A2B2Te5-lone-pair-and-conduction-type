# 后处理复现

安装 requirements.txt 的依赖，以交付目录为当前目录：

```bash
python 10_脚本/analyze_sources.py .
python 10_脚本/extra_audit.py .
python 10_脚本/build_report.py .
python 10_脚本/final_audit.py .
```

build_report.py 与 paper_revision.py 共同生成主报告、全元素对图及共享/ET/ELF 图。final_audit.py 调用 audit_sharing.py，核对表格、原矩阵算术、图片和链接，再生成 SHA256。

共享与局域化/ 内脚本保存计算来源，使用原工作区绝对路径；不属于无需原工作区即可重跑的 DFT 流程。DI 的对称展开与分摊近似尚未验证，不能据算术一致宣称严格 ES。全套 WAVECAR、POTCAR 不放入本交付。

原始密度以 gzip 保存于 06_COGITO成键分析与COHP/源数据/Bader_TET_ES原始密度。解压后可按 ELF 阈值复核密度积分。原网格积分检查见 09_审查/ELF原网格积分校验.json。

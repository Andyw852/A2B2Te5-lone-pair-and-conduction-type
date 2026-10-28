"""把论文报告的 HTML 模板中的 {{FIG:名称}} 换成嵌入的图片（缩到 1800 px 宽的 PNG，base64），输出单个自包含 HTML。

用法：python 10_脚本/build_paper_report.py
输入：21_阶段报告与交接/论文报告_20261028.src.html、21_阶段报告与交接/论文图_20261028/*.png
输出：21_阶段报告与交接/论文报告_20261028.html
"""
import base64, io, os, re
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, '21_阶段报告与交接')
src = open(os.path.join(D, '论文报告_20261028.src.html'), encoding='utf-8').read()


def embed(m):
    im = Image.open(os.path.join(D, '论文图_20261028', m.group(1) + '.png')).convert('RGB')
    if im.width > 1800:
        im = im.resize((1800, round(im.height * 1800 / im.width)), Image.LANCZOS)
    buf = io.BytesIO(); im.save(buf, 'PNG', optimize=True)
    return 'data:image/png;base64,' + base64.b64encode(buf.getvalue()).decode()


html = re.sub(r'\{\{FIG:([^}]+)\}\}', embed, src)
out = os.path.join(D, '论文报告_20261028.html')
open(out, 'w', encoding='utf-8').write(html)
print(out, f'{len(html) / 1e6:.1f} MB')

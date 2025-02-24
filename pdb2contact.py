# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 11:23:03 2025

@author: ijinh
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from scipy.spatial.distance import pdist, squareform

# 获取PDB文件路径和参数
pdb_file = sys.argv[1]  # PDB 文件路径
distType = sys.argv[2]  # "T" 表示真实距离矩阵, 其他表示 Contact Map

# 解析PDB文件
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", pdb_file)
model = structure[0]

# 获取 C-alpha 原子坐标
residues = [res for res in model.get_residues() if res.has_id("CA")]
coords = np.array([res["CA"].coord for res in residues])

# 计算真实距离矩阵
dist_matrix = squareform(pdist(coords))

# 计算 Contact Map（阈值 6Å）
contact_map = (dist_matrix < 8).astype(int)

# 选择输出矩阵
if distType == "T":
    mt = dist_matrix  # 真实距离矩阵
    fmt = "%.3f"  # 保留 3 位小数
else:
    mt = contact_map  # Contact Map
    fmt = "%d"  # 二值化

# 生成输出文件名（去掉 .pdb 后缀）
output_prefix = pdb_file.rsplit(".", 1)[0]

# 保存矩阵
np.savetxt(output_prefix+ "."+distType + ".cm.txt", mt, fmt=fmt)

# 绘制矩阵（无刻度，全屏）
plt.figure(figsize=(12, 12))
plt.imshow(mt, cmap="gray_r", origin="lower")
plt.xticks([])
plt.yticks([])
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)  # 让图像铺满
plt.savefig(output_prefix + "."+distType+".cm.png", bbox_inches="tight", dpi=300)
plt.show()

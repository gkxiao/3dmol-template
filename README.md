## 3dmol.js 打包展示蛋白配体相互作用

建议一个项目一个独立文件夹（避免多个 .pdb / .sdf 时脚本只取排序第一个）：

1. 复制 3 个程序文件过去：index.html、3Dmol-min.js、build_data.py
2. 放入新的 *.pdb 和 *.sdf（各一个）
3. 在该目录执行：python build_data.py
   多个文件时显式指定：python build_data.py receptor.pdb hits.sdf
   → 重新生成 data.js，刷新页面即可；最后把文件夹打包分发给同事
4. 分发：将index.html、3Dmol-min.js、data.js打包，其它文件不需要

比如：
```bash
python build_data.py
已生成 data.js (790.2 KB)
  蛋白: 4ivc_receptor.pdb | 原子 4867 | 链 A | 分辨率 2.35 Å
  配体: 51 个 pose / 36 个分子 (来自 test.sdf)
  现在可直接双击 index.html（离线、无需服务器）
```

build_data.py 只用 Python 标准库，且只在制作者换数据时跑一次；同事端仍只需浏览器，不需要 Python 或任何包。

数据准备的两个硬要求：

- 坐标系必须一致——新配体的 pose 必须已经对接/叠合到这个新蛋白的坐标框架里，否则会看到配体飘在蛋白外。
- PDB 建议单模型、删掉不想显示的原有共晶配体/水（HETATM 会一起被读入并显示）；SDF 必须是 3D 坐标、每条以 $$$$ 结尾，V2000/V3000 都行。第一行是分子名，同名记录自动归为同一分子的不同 pose；想每个分子独立成组就用互不相同的名字。

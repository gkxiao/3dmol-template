## 3dmol.js 打包展示蛋白配体相互作用

建议一个项目一个独立文件夹（避免多个 .pdb / .sdf 时脚本只取排序第一个）：

1. 复制 3 个程序文件过去：index.html、3Dmol-min.js、build_data.py
2. 放入新的 *.pdb 和 *.sdf（各一个）
3. 在该目录执行：python build_data.py
   多个文件时显式指定：python build_data.py receptor.pdb hits.sdf
   → 重新生成 data.js，刷新页面即可；最后把文件夹打包分发给同事

比如：
```bash
python build_data.py
已生成 data.js (790.2 KB)
  蛋白: 4ivc_receptor.pdb | 原子 4867 | 链 A | 分辨率 2.35 Å
  配体: 51 个 pose / 36 个分子 (来自 test.sdf)
  现在可直接双击 index.html（离线、无需服务器）
```

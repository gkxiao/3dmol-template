#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_data.py - 把同目录的 PDB + SDF 打包成浏览器可直接读取的 data.js

用法:
    python build_data.py                       # 自动查找同目录的 *.pdb / *.sdf
    python build_data.py 7s8l.pdb my.sdf       # 显式指定文件
    python build_data.py --v2000               # 额外输出一份保留全氢的 V2000 SDF

输出:
    data.js  (window.MOL_DATA = {...})

为什么要把结构文本内嵌进 .js:
    浏览器用 file:// 双击打开 index.html 时, JS 的 fetch/XHR 无法读取同目录的
    本地文件 (CORS 限制), 但 <script src="data.js"></script> 可以。这样整个文件夹
    压缩后发给同事, 双击即可用, 不需要联网、不需要启动本地服务器。

数据是纯文本, 不改动原始 PDB/SDF; 换一批分子后重跑本脚本即可。
"""

import json
import re
import sys
import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent

# 这些字段是 Cresset 的场点坐标块, 体积大且对展示无意义, 不写入 data.js
SKIP_PROPS = {
    "_cresset_fieldpoints",
    "_cresset_fieldpoints2",
    "_cresset_fieldconstraints",
}


def pick_file(args, ext):
    """按命令行参数或目录扫描挑一个输入文件"""
    for a in args:
        if a.lower().endswith(ext):
            p = Path(a)
            if not p.is_absolute():
                p = HERE / p
            if p.exists():
                return p
    found = sorted(HERE.glob("*" + ext))
    if not found:
        raise SystemExit("找不到 *%s 文件，请把它放在本脚本同目录" % ext)
    if len(found) > 1:
        print("提示: 发现多个 %s，使用 %s" % (ext, found[0].name))
    return found[0]


def split_sdf(text):
    """按 $$$$ 切成单条记录，保留每条的 '$$$$' 结束符（3Dmol 的 SDF 解析器需要）"""
    out = []
    for part in text.split("$$$$"):
        part = part.strip("\r\n")
        if part.strip():
            out.append(part + "\n$$$$\n")
    return out


def parse_props(lines):
    """解析 SDF 的属性块 '>  <Field>' -> 值（可能多行）"""
    props = {}
    i = 0
    while i < len(lines):
        m = re.match(r">\s*<([^>]+)>", lines[i])
        if not m:
            i += 1
            continue
        key = m.group(1).strip()
        vals = []
        i += 1
        while i < len(lines) and not lines[i].startswith(">"):
            if lines[i].strip() != "$$$$":
                vals.append(lines[i])
            i += 1
        while vals and not vals[-1].strip():
            vals.pop()
        while vals and not vals[0].strip():
            vals.pop(0)
        if key not in SKIP_PROPS:
            props[key] = "\n".join(vals).strip()
    return props


def geometry(lines):
    """返回 (重原子数, 重原子质心 [x,y,z])，同时兼容 V3000 与 V2000"""
    v3000 = any("V3000" in l for l in lines[:6])
    n = 0
    sx = sy = sz = 0.0

    def add(elem, x, y, z):
        nonlocal n, sx, sy, sz
        if elem.upper() == "H":
            return
        n += 1
        sx += x
        sy += y
        sz += z

    if v3000:
        in_atom = False
        for l in lines:
            s = l.strip()
            if s.startswith("M  V30 BEGIN ATOM"):
                in_atom = True
                continue
            if s.startswith("M  V30 END ATOM"):
                in_atom = False
                continue
            if in_atom and s.startswith("M  V30 "):
                t = s[7:].split()
                if len(t) >= 5:
                    try:
                        add(t[1], float(t[2]), float(t[3]), float(t[4]))
                    except ValueError:
                        pass
    else:
        try:
            counts = next(l for l in lines if re.match(r"^\s*\d+\s+\d+", l))
            na = int(counts[0:3])
            start = lines.index(counts) + 1
            for l in lines[start:start + na]:
                try:
                    add(l[31:34].strip(), float(l[0:10]), float(l[10:20]), float(l[20:30]))
                except (ValueError, IndexError):
                    pass
        except StopIteration:
            pass

    center = [round(sx / n, 3), round(sy / n, 3), round(sz / n, 3)] if n else [0, 0, 0]
    return n, center


def v3000_to_v2000(rec):
    """把一条 V3000 记录转成 V2000（保留全部原子，包括所有显式氢、键级与电荷）。
    仅用于生成兼容副本；本查看器直接使用原始 V3000 即可。"""
    lines = rec.splitlines()
    title = lines[0] if lines else "MOL"
    na = nb = 0
    try:
        counts = next(l for l in lines if l.startswith("M  V30 COUNTS"))
        t = counts[13:].split()
        na, nb = int(t[0]), int(t[1])
    except (StopIteration, ValueError, IndexError):
        return None

    atoms, bonds, charges, imap = [], [], [], {}
    section = None
    for l in lines:
        s = l.strip()
        if s.startswith("M  V30 BEGIN ATOM"): section = "ATOM"; continue
        if s.startswith("M  V30 END ATOM"): section = None; continue
        if s.startswith("M  V30 BEGIN BOND"): section = "BOND"; continue
        if s.startswith("M  V30 END BOND"): section = None; continue
        if not s.startswith("M  V30 "): continue
        t = s[7:].split()
        if section == "ATOM" and len(t) >= 5:
            try:
                oid, elem = int(t[0]), t[1]
                x, y, z = float(t[2]), float(t[3]), float(t[4])
            except ValueError:
                continue
            imap[oid] = len(atoms) + 1
            atoms.append((elem, x, y, z))
            for kv in t[5:]:
                if kv.startswith("CHG="):
                    try: charges.append((imap[oid], int(kv.split("=")[1])))
                    except ValueError: pass
        elif section == "BOND" and len(t) >= 4:
            try:
                a1, a2, order = imap[int(t[2])], imap[int(t[3])], int(float(t[1]))
                bonds.append((a1, a2, order))
            except (KeyError, ValueError):
                pass

    if len(atoms) != na or len(bonds) != nb:
        return None                      # 结构不完整则不输出，避免生成错误文件

    out = [title, "  WorkBuddy 3D", ""]
    out.append("%3d%3d  0  0  0  0  0  0  0  0999 V2000" % (na, nb))
    for elem, x, y, z in atoms:
        out.append("%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0" % (x, y, z, elem))
    for a1, a2, o in bonds:
        out.append("%3d%3d%3d  0  0  0  0" % (a1, a2, o))
    for i in range(0, len(charges), 8):
        chunk = charges[i:i + 8]
        out.append("M  CHG%3d" % len(chunk) + "".join("%3d%3d" % (a, c) for a, c in chunk))
    out.append("M  END")
    pstart = next((i for i, l in enumerate(lines) if l.startswith(">")), None)
    if pstart is not None:
        for l in lines[pstart:]:
            if l.strip() == "$$$$": break
            out.append(l)
    out.append("$$$$")
    return "\n".join(out) + "\n"


def write_v2000(sdf_file, recs):
    outdir = HERE / "v2000"
    outdir.mkdir(exist_ok=True)
    outfile = outdir / (sdf_file.stem + "_v2000.sdf")
    converted, skipped = [], 0
    for rec in recs:
        v = v3000_to_v2000(rec)
        if v is None: skipped += 1
        else: converted.append(v)
    outfile.write_text("".join(converted), encoding="utf-8", newline="\n")
    print("  V2000 副本（保留全氢）: %s  %d 条记录%s  %.1f KB" % (
        outfile.relative_to(HERE), len(converted),
        ("，跳过 %d 条非 V3000" % skipped) if skipped else "",
        outfile.stat().st_size / 1024))


def pdb_meta(text):
    lines = text.splitlines()
    title = " ".join(l[10:].strip() for l in lines if l.startswith("TITLE")).strip()
    res = ""
    for l in lines:
        if l.startswith("REMARK   2 RESOLUTION."):
            res = l.split("RESOLUTION.")[1].strip().split()[0]
            break
    chains = sorted({l[21] for l in lines if l.startswith(("ATOM  ", "HETATM")) and l[21].strip()})
    atoms = [l for l in lines if l.startswith(("ATOM  ", "HETATM"))]
    rng = ""
    if atoms:
        rng = "%s %s - %s %s" % (
            atoms[0][17:20].strip(), atoms[0][22:26].strip(),
            atoms[-1][17:20].strip(), atoms[-1][22:26].strip(),
        )
    return title, res, chains, rng, len(atoms)


def js_safe(obj):
    """JSON -> JS 字面量；转义 '<' 防止数据里出现 </script> 提前闭合脚本标签"""
    return (json.dumps(obj, ensure_ascii=False, separators=(",", ":"))
            .replace("<", "\\u003c")
            .replace("\u2028", "\\u2028")
            .replace("\u2029", "\\u2029"))


def main():
    args = sys.argv[1:]
    want_v2000 = "--v2000" in args
    args = [a for a in args if a != "--v2000"]
    pdb_file = pick_file(args, ".pdb")
    sdf_file = pick_file(args, ".sdf")

    protein = pdb_file.read_text(encoding="utf-8", errors="replace")
    sdf_text = sdf_file.read_text(encoding="utf-8", errors="replace")
    recs = split_sdf(sdf_text)
    if not recs:
        raise SystemExit("SDF 中没有读到任何记录")

    title, resolution, chains, res_range, n_atoms = pdb_meta(protein)
    pose_counter = {}
    ligands = []
    for i, rec in enumerate(recs):
        lines = rec.splitlines()
        name = lines[0].strip() if lines and lines[0].strip() else "MOL-%d" % (i + 1)
        group = name  # 同名 = 同一分子的不同 pose
        pose_counter[group] = pose_counter.get(group, 0) + 1
        n_heavy, center = geometry(lines)
        ligands.append({
            "id": i,
            "name": name,
            "group": group,
            "pose": pose_counter[group],
            "mol": rec,
            "props": parse_props(lines),
            "nHeavy": n_heavy,
            "center": center,
        })

    meta = {
        "pdbFile": pdb_file.name,
        "sdfFile": sdf_file.name,
        "pdbId": re.sub(r"\.pdb$", "", pdb_file.name, flags=re.I).upper(),
        "title": title,
        "resolution": resolution,
        "chains": chains,
        "resRange": res_range,
        "nProteinAtoms": n_atoms,
        "nLigands": len(ligands),
        "nGroups": len(pose_counter),
        "generated": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "viewer": "3Dmol.js 2.5.5 (offline)",
    }

    out = HERE / "data.js"
    with out.open("w", encoding="utf-8", newline="\n") as f:
        f.write("/* 自动生成于 %s by build_data.py —— 请勿手工编辑 */\n" % meta["generated"])
        f.write("window.MOL_DATA = {\n")
        f.write('  "meta": %s,\n' % js_safe(meta))
        f.write('  "protein": %s,\n' % js_safe(protein))
        f.write('  "ligands": %s\n' % js_safe(ligands))
        f.write("};\n")

    print("已生成 %s (%.1f KB)" % (out.name, out.stat().st_size / 1024))
    print("  蛋白: %s | 原子 %d | 链 %s | 分辨率 %s Å" % (pdb_file.name, n_atoms, ",".join(chains), resolution or "?"))
    print("  配体: %d 个 pose / %d 个分子 (来自 %s)" % (len(ligands), len(pose_counter), sdf_file.name))
    if want_v2000:
        write_v2000(sdf_file, recs)
    print("  现在可直接双击 index.html（离线、无需服务器）")


if __name__ == "__main__":
    main()

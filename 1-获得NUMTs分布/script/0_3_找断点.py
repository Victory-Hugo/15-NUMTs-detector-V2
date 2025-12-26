#!/usr/bin/env python3
# 2025年5月
# -*- coding: utf-8 -*-
"""
find_breakpoints.py
===================

Detect nuclear–mtDNA (NUMT) breakpoints from a single PSL slice，
可一次处理多条染色体：输入 chr 参数时用逗号分隔，
脚本遇到某条染色体无比对时仅打印警告并跳过，不会停止整个运行。

用法
----
python3 find_breakpoints.py <psl_part> <sampleID> <chrom_list> <start> <end> <output_prefix>

产生多个文件，每条染色体各自输出：
    <output_prefix>.<CHR>.AllBreakpoints.tsv
    <output_prefix>.<CHR>.ConfidentBreakpoints.tsv
    ……
"""

import sys
import json
import os
import pandas as pd

###############################################################################
# ------------------------------- 参数解析 ---------------------------------- #
###############################################################################
try:
    INPUT_PSL, SAMPLEID, CHR_ARG, START, END, OUTPUT_PREFIX = sys.argv[1:]
    START, END = int(START), int(END)
except ValueError:
    sys.exit(
        "❌ 参数不足！\n"
        "    python3 find_breakpoints.py <psl_part> <sampleID> <chrom_list> <start> <end> <output_prefix>"
    )

def load_config():
    cfg_path = os.environ.get("NUMTS_CONFIG")
    if not cfg_path:
        conf_dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "confs")
        )
        preferred = os.path.join(conf_dir, "GRCH38.json")
        fallback = os.path.join(conf_dir, "numts_config.json")
        cfg_path = preferred if os.path.exists(preferred) else fallback
    try:
        with open(cfg_path, "r") as f:
            return json.load(f)
    except FileNotFoundError:
        sys.exit(f"❌ 找不到配置文件: {cfg_path}")
    except json.JSONDecodeError as exc:
        sys.exit(f"❌ 配置文件 JSON 解析失败: {cfg_path} ({exc})")


bp_cfg = load_config().get("breakpoints", {})
TNAME_MAP = bp_cfg.get("tname_map", {})
MT_CHR = bp_cfg.get("mt_chrom")
CHR_PREFIX = bp_cfg.get("chrom_prefix", "")
CHR_ALIASES = bp_cfg.get("chrom_aliases", {})

if not MT_CHR:
    sys.exit("❌ 配置缺少 breakpoints.mt_chrom")
if not isinstance(TNAME_MAP, dict):
    sys.exit("❌ 配置错误：breakpoints.tname_map 需要是字典")

def normalize_chr(raw: str) -> str:
    raw = raw.strip()
    if not raw:
        return ""
    if raw in CHR_ALIASES:
        return CHR_ALIASES[raw]
    if CHR_PREFIX and not raw.startswith(CHR_PREFIX):
        return CHR_PREFIX + raw
    return raw

# 解析多条染色体，逗号分隔
raw_chrs = CHR_ARG.split(',')
chr_list = []
for raw in raw_chrs:
    c = normalize_chr(raw)
    if c:
        chr_list.append(c)

###############################################################################
# ---------------------------- 辅助函数与常量 ------------------------------- #
###############################################################################
PSL_COLS = [
    "match", "misMatch", "repMatch", "Ns", "QgapCount", "QgapBases",
    "TgapCount", "TgapBases", "strand", "Qname", "Qsize", "Qstart",
    "Qend", "Tname", "Tsize", "Tstart", "Tend", "blockCount",
    "blockSizes", "qStarts", "tStarts"
]

def classify_breakpoint(row, which='nu', read_len=150, mismatch=5):
    """给每条 PSL 命中分配断点类别。"""
    if which == 'nu':
        if   row.strand == "+" and row.Qend   >= read_len - mismatch: return "nu_Tstart_Bright"
        elif row.strand == "+" and row.Qstart <= mismatch           : return "nu_Tend_Bleft"
        elif row.strand == "-"                                          : return "nu_NegStrand"
        else                                                            : return "nu_useLess"
    else:
        if   row.strand == "+" and row.Qend   >= read_len - mismatch: return "mt_Tstart"
        elif row.strand == "+" and row.Qstart <= mismatch           : return "mt_Tend"
        elif row.strand == "-" and row.Qend   >= read_len - mismatch: return "mt_Tend"
        elif row.strand == "-" and row.Qstart <= mismatch           : return "mt_Tstart"
        else                                                        : return "mt_useLess"

def group_and_count(df, cols):
    """按给定列分组并统计 reads 数量。"""
    return df.groupby(cols, dropna=False).size().reset_index(name="readsCount")

def unify_pos(x):
    """统一取 Tend 或 Tstart 作为断点位置。"""
    return x.Tend if pd.notna(x.Tend) else x.Tstart

###############################################################################
# ------------------------------- 数据读取 ---------------------------------- #
###############################################################################
# 读入 PSL
df = pd.read_csv(INPUT_PSL, skiprows=5, sep="\t", names=PSL_COLS)
df['matchLEN'] = df['Tend'] - df['Tstart']

# 全局基本过滤
df = df[(df.matchLEN <= 150) & (df.misMatch <= 5)]
df = df[(df.Tend     >= 150) | (df.Tend     <= 5)]

# 染色体名称映射
df['Tname_mapped'] = df.Tname.map(TNAME_MAP).fillna(df.Tname)

# 先提取线粒体数据，后面循环里可用
mt_df = df[df.Tname_mapped == MT_CHR].copy()
mt_df['pointGroup'] = mt_df.apply(classify_breakpoint, axis=1, which='mt')
mt_df['pointGroup'] = mt_df['pointGroup'].fillna('').astype(str)
mt_df['Group']      = mt_df['pointGroup'].str.replace(r'_T.*B', '', regex=True)
mt_df['chr']        = MT_CHR

# 对每条指定的染色体依次处理
for CHR in chr_list:
    nu_df = df[df.Tname_mapped == CHR].copy()
    if nu_df.empty:
        print(f"⚠️  {SAMPLEID}: 在 {CHR}:{START}-{END} 未检测到任何核序列比对，已跳过。")
        continue

    # 分类
    nu_df['pointGroup'] = nu_df.apply(classify_breakpoint, axis=1, which='nu')
    nu_df['Group']      = nu_df.pointGroup.str.replace(r'_T.*B', '', regex=True)
    nu_df['chr']        = CHR

    # 1) 核 DNA 左/右断点
    nu_left  = group_and_count(nu_df[nu_df.pointGroup=="nu_Tend_Bleft"  ],
                               ['pointGroup','Group','chr','Tend' , 'strand'])
    nu_right = group_and_count(nu_df[nu_df.pointGroup=="nu_Tstart_Bright"],
                               ['pointGroup','Group','chr','Tstart','strand'])
    nu_both  = pd.concat([nu_left, nu_right], ignore_index=True)

    # 2) 跨映射核断点（nu_mt_both）
    nu_mt   = nu_df[nu_df.Qname.isin(mt_df.Qname)]
    nu_mt_l = group_and_count(nu_mt[nu_mt.pointGroup=="nu_Tend_Bleft"  ],
                              ['pointGroup','Group','chr','Tend' , 'strand'])
    nu_mt_r = group_and_count(nu_mt[nu_mt.pointGroup=="nu_Tstart_Bright"],
                              ['pointGroup','Group','chr','Tstart','strand'])
    nu_mt_both = pd.concat([nu_mt_l, nu_mt_r], ignore_index=True)

    # 3) 线粒体自身断点 (mt_both 已在循环外算过，也可在此重复)
    mt_tend   = group_and_count(mt_df[mt_df.pointGroup=="mt_Tend"  ],
                                ['pointGroup','Group','chr','Tend'  ,'strand'])
    mt_tstart = group_and_count(mt_df[mt_df.pointGroup=="mt_Tstart"],
                                ['pointGroup','Group','chr','Tstart','strand'])
    mt_both   = pd.concat([mt_tend, mt_tstart], ignore_index=True)

    # 4) 跨映射线粒体断点 (mt_confident)
    mt_conf   = mt_df[mt_df.Qname.isin(nu_df.Qname)]
    mt_conf_l = group_and_count(mt_conf[mt_conf.pointGroup=="mt_Tstart"],
                                ['pointGroup','chr','Tstart','strand'])
    mt_conf_r = group_and_count(mt_conf[mt_conf.pointGroup=="mt_Tend"  ],
                                ['pointGroup','chr','Tend'  ,'strand'])
    mt_conf   = pd.concat([mt_conf_l, mt_conf_r], ignore_index=True) \
                 .assign(Group='mt_confident')

    # 合并 & 输出
    all_breaks = (
        pd.concat([nu_both, mt_both], ignore_index=True)
          .assign(pos=lambda d: d.apply(unify_pos, axis=1), type='all')
    )
    conf_breaks = (
        pd.concat([nu_mt_both, mt_conf], ignore_index=True)
          .assign(pos=lambda d: d.apply(unify_pos, axis=1), type='confident')
    )

    cols = ['sampleID','chr','strand','pointGroup','Group','readsCount','pos','type']
    all_breaks  = all_breaks .assign(sampleID=SAMPLEID)[cols]
    conf_breaks = conf_breaks.assign(sampleID=SAMPLEID)[cols]

    # 文件名包含染色体
    all_path  = f"{OUTPUT_PREFIX}.{CHR}.AllBreakpoints.tsv"
    conf_path = f"{OUTPUT_PREFIX}.{CHR}.ConfidentBreakpoints.tsv"

    all_breaks .to_csv(all_path , sep='\t', index=False)
    conf_breaks.to_csv(conf_path, sep='\t', index=False)

    print(f"✅ 完成 {SAMPLEID} {CHR}:")
    print(f"   → {all_path}")
    print(f"   → {conf_path}")

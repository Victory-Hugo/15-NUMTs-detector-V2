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
python3 find_breakpoints.py <psl_part> <sampleID> <chr1,chr2,…> <start> <end> <output_prefix>

产生多个文件，每条染色体各自输出：
    <output_prefix>.chr1.AllBreakpoints.tsv
    <output_prefix>.chr1.ConfidentBreakpoints.tsv
    <output_prefix>.chr2.AllBreakpoints.tsv
    ……
"""

import sys
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
        "    python3 find_breakpoints.py <psl_part> <sampleID> <chr1,chr2,…> <start> <end> <output_prefix>"
    )

# 解析多条染色体，逗号分隔
raw_chrs = CHR_ARG.split(',')
chr_list = []
for raw in raw_chrs:
    c = raw if raw.startswith("chr") else "chr" + raw
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

TNAME_MAP = {
    'NC_012920.1':'chrM','NC_000001.11':'chr1','NC_000002.12':'chr2',
    'NC_000003.12':'chr3','NC_000004.12':'chr4','NC_000005.10':'chr5',
    'NC_000006.12':'chr6','NC_000007.14':'chr7','NC_000008.11':'chr8',
    'NC_000009.12':'chr9','NC_000010.11':'chr10','NC_000011.10':'chr11',
    'NC_000012.12':'chr12','NC_000013.11':'chr13','NC_000014.9':'chr14',
    'NC_000015.10':'chr15','NC_000016.10':'chr16','NC_000017.11':'chr17',
    'NC_000018.10':'chr18','NC_000019.10':'chr19','NC_000020.11':'chr20',
    'NC_000021.9':'chr21','NC_000022.11':'chr22','NC_000023.11':'chrX',
    'NC_000024.10':'chrY'
}

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
mt_df = df[df.Tname_mapped == 'chrM'].copy()
mt_df['pointGroup'] = mt_df.apply(classify_breakpoint, axis=1, which='mt')
mt_df['pointGroup'] = mt_df['pointGroup'].fillna('').astype(str)
mt_df['Group']      = mt_df['pointGroup'].str.replace(r'_T.*B', '', regex=True)
mt_df['chr']        = 'chrM'

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

# #!/usr/bin/env python
# 2025年2月
# ################################################################################
# ## This script look for NUMT breakpoints
# ################################################################################

# import sys, os
# import pandas as pd

# def classify_breakpoint(row, type='nu'):
#     mismatchLEN = 5 # 默认是3 请在这里选择过滤措施，选择更大的数量会得到更多的结果
#     readLEN = 150 # 默认是150 请在这里选择过滤措施，选择更长的长度会得到更多的结果
#     if type == 'nu':
#         if row['strand'] == "+" and row['Qend'] >= readLEN - mismatchLEN:
#             return "nu_Tstart_Bright"
#         elif row['strand'] == "+" and row['Qstart'] <= mismatchLEN:
#             return "nu_Tend_Bleft"
#         elif row['strand'] == "-":
#             return "nu_NegStrand"
#         else:
#             return "nu_useLess"
#     else:
#         if row['strand'] == "+" and row['Qend'] >= readLEN - mismatchLEN:
#             return "mt_Tstart"
#         elif row['strand'] == "+" and row['Qstart'] <= mismatchLEN:
#             return "mt_Tend"
#         elif row['strand'] == "-" and row['Qend'] >= readLEN - mismatchLEN:
#             return "mt_Tend"
#         elif row['strand'] == "-" and row['Qstart'] <= mismatchLEN:
#             return "mt_Tstart"
#         else:
#             return "mt_useLess"

# # 从命令行读取输入参数
# INPUT_PSL, SAMPLEID, CHR, START, END, OUTPUT = sys.argv[1:]
# START, END = int(START), int(END)

# # 读取数据并预处理
# psl_columns = [
#     "match", "misMatch", "repMatch", "Ns", "QgapCount", "QgapBases", "TgapCount", "TgapBases", "strand",
#     "Qname", "Qsize", "Qstart", "Qend", "Tname", "Tsize", "Tstart", "Tend", "blockCount", "blockSizes",
#     "qStarts", "tStarts"
# ]
# df = pd.read_csv(INPUT_PSL, skiprows=5, sep="\t", names=psl_columns)
# df['matchLEN'] = df['Tend'] - df['Tstart']

# # 过滤数据，过滤措施
# # ! 默认是df[(df['matchLEN'] < 140) & (df['misMatch'] <= 3)]请在这里选择过滤措施，选择更大的数量会得到更多的结果
# filtered_df = df[(df['matchLEN'] < 150) & (df['misMatch'] <= 5)] 
# # ! 默认是[(filtered_df['Tend'] >= 147) | (filtered_df['Tend'] <= 3)]请在这里选择过滤措施，选择更长的长度会得到更多的结果
# filtered_df = filtered_df[(filtered_df['Tend'] >= 147) | (filtered_df['Tend'] <= 3)]
# print(filtered_df.columns.tolist())
# print("所有染色体:", filtered_df['Tname'].unique())

# # Tname映射
# #! 原始脚本并不存在这一步，但是由于GRCH38的参考序列对染色体名称进行了更换，所以在这里将名称替换回来。
# tname_mapping = {
#     'NC_012920.1': 'MT',
#     'NC_000001.11': 'chr1',
#     'NC_000002.12': 'chr2',
#     'NC_000003.12': 'chr3',
#     'NC_000004.12': 'chr4',
#     'NC_000005.10': 'chr5',
#     'NC_000006.12': 'chr6',
#     'NC_000007.14': 'chr7',
#     'NC_000008.11': 'chr8',
#     'NC_000009.12': 'chr9',
#     'NC_000010.11': 'chr10',
#     'NC_000011.10': 'chr11',
#     'NC_000012.12': 'chr12',
#     'NC_000013.11': 'chr13',
#     'NC_000014.9': 'chr14',
#     'NC_000015.10': 'chr15',
#     'NC_000016.10': 'chr16',
#     'NC_000017.11': 'chr17',
#     'NC_000018.10': 'chr18',
#     'NC_000019.10': 'chr19',
#     'NC_000020.11': 'chr20',
#     'NC_000021.9': 'chr21',
#     'NC_000022.11': 'chr22',
#     'NC_000023.11': 'chrX',
#     'NC_000024.10': 'chrY'
# }
# filtered_df = filtered_df.copy()
# filtered_df.loc[:, 'Tname_mapped'] = filtered_df['Tname'].map(tname_mapping).fillna(filtered_df['Tname'])

# # 分离线粒体和核序列
# mt_df = filtered_df.loc[filtered_df['Tname_mapped'] == 'MT'].copy()
# #! 原始代码如下：
# # nu_df = filtered_df.loc[(filtered_df['Tname_mapped'] == CHR) & 
# #                         (filtered_df['Tstart'] >= START) & 
# #                         (filtered_df['Tend'] <= END)].copy()
# #! 放宽要求：给定一个pad,弹性窗口。设定在50至200之间。
# # pad = 2000
# # nu_df = filtered_df[
# #     (filtered_df['Tname_mapped'] == CHR) &
# #     (filtered_df['Tstart'] >= START - pad) &
# #     (filtered_df['Tend']   <= END   + pad)
# # ].copy()

# #! 最宽要求：只要存在于染色体就算
# nu_df = filtered_df.loc[
#      (filtered_df['Tname_mapped'] == CHR)].copy()

# if not nu_df.empty:
#     # 应用分类函数
#     nu_df.loc[:, 'pointGroup'] = nu_df.apply(classify_breakpoint, axis=1, type='nu')
#     if not mt_df.empty:
#         mt_df.loc[:, 'pointGroup'] = mt_df.apply(classify_breakpoint, axis=1, type='mt')

#     nu_df.loc[:, 'Group'] = nu_df['pointGroup'].str.replace(r'_T.*B', '', regex=True)
#     nu_df.loc[:, 'chr'] = CHR
#     if not mt_df.empty:
#         mt_df.loc[:, 'chr'] = 'MT'

#     def group_and_count(df, by_cols):
#         return df.groupby(by_cols).size().reset_index(name="readsCount")

#     # 核DNA断点
#     nu_left = group_and_count(nu_df.loc[nu_df['pointGroup'] == 'nu_Tend_Bleft'], ['pointGroup', 'Group', 'chr', 'Tend', 'strand'])
#     nu_right = group_and_count(nu_df.loc[nu_df['pointGroup'] == 'nu_Tstart_Bright'], ['pointGroup', 'Group', 'chr', 'Tstart', 'strand'])
#     nu_both = pd.concat([nu_left, nu_right])

#     # 同时映射到线粒体的核DNA断点
#     nu_mt = nu_df.loc[nu_df['Qname'].isin(mt_df['Qname'])]
#     nu_mt_left = group_and_count(nu_mt.loc[nu_mt['pointGroup'] == 'nu_Tend_Bleft'], ['pointGroup', 'Group', 'chr', 'Tend', 'strand'])
#     nu_mt_right = group_and_count(nu_mt.loc[nu_mt['pointGroup'] == 'nu_Tstart_Bright'], ['pointGroup', 'Group', 'chr', 'Tstart', 'strand'])
#     nu_mt_both = pd.concat([nu_mt_left, nu_mt_right])

#     # 线粒体断点
#     if not mt_df.empty:
#         mt_tend = group_and_count(mt_df.loc[mt_df['pointGroup'] == 'mt_Tend'], ['pointGroup', 'chr', 'Tend', 'strand'])
#         mt_tstart = group_and_count(mt_df.loc[mt_df['pointGroup'] == 'mt_Tstart'], ['pointGroup', 'chr', 'Tstart', 'strand'])
#         mt_both = pd.concat([mt_tend, mt_tstart])
#         mt_both['Group'] = 'UKn'

#         # 同时映射到核DNA的线粒体断点
#         mt_conf = pd.concat([
#             group_and_count(mt_df.loc[(mt_df['Qname'].isin(nu_df.loc[nu_df['pointGroup'] == 'nu_Tend_Bleft', 'Qname'])) & (mt_df['pointGroup'] == 'mt_Tstart')], ['pointGroup', 'chr', 'Tstart', 'strand']),
#             group_and_count(mt_df.loc[(mt_df['Qname'].isin(nu_df.loc[nu_df['pointGroup'] == 'nu_Tend_Bleft', 'Qname'])) & (mt_df['pointGroup'] == 'mt_Tend')], ['pointGroup', 'chr', 'Tend', 'strand']),
#             group_and_count(mt_df.loc[(mt_df['Qname'].isin(nu_df.loc[nu_df['pointGroup'] == 'nu_Tstart_Bright', 'Qname'])) & (mt_df['pointGroup'] == 'mt_Tstart')], ['pointGroup', 'chr', 'Tstart', 'strand']),
#             group_and_count(mt_df.loc[(mt_df['Qname'].isin(nu_df.loc[nu_df['pointGroup'] == 'nu_Tstart_Bright', 'Qname'])) & (mt_df['pointGroup'] == 'mt_Tend')], ['pointGroup', 'chr', 'Tend', 'strand'])
#         ])

#         if not mt_conf.empty:
#             # 确保 'pointGroup' 列的数据类型为字符串，并且不存在缺失值
#             mt_conf = mt_conf.dropna(subset=['pointGroup'])
#             mt_conf['pointGroup'] = mt_conf['pointGroup'].astype(str)

#             # 使用矢量化方法进行赋值
#             mt_conf['Group'] = mt_conf['pointGroup'].apply(lambda x: 'mtLeft' if 'left' in x else 'mtRight')
#         else:
#             print("mt_conf 数据框为空，跳过 'Group' 列的赋值。")

#         all_breakpoints = pd.concat([nu_both, mt_both])
#         mt_conf['sampleID'] = SAMPLEID
#         mt_conf['Tstart'] = mt_conf['Tstart'].fillna(-1).astype(int)
#         mt_conf['Tend'] = mt_conf['Tend'].fillna(-1).astype(int)
#     else:
#         all_breakpoints = nu_both
#         mt_conf = pd.DataFrame()

#     # 输出结果
#     all_breakpoints['sampleID'] = SAMPLEID
#     all_breakpoints['Tstart'] = all_breakpoints['Tstart'].fillna(-1).astype(int)
#     all_breakpoints['Tend'] = all_breakpoints['Tend'].fillna(-1).astype(int)

#     confident_breakpoints = pd.concat([nu_mt_both, mt_conf])
#     confident_breakpoints['sampleID'] = SAMPLEID
#     confident_breakpoints['Tstart'] = confident_breakpoints['Tstart'].fillna(-1).astype(int)
#     confident_breakpoints['Tend'] = confident_breakpoints['Tend'].fillna(-1).astype(int)

#     # 写入文件
#     confident_breakpoints.to_csv(OUTPUT + '.Breakpoints.tsv', sep='\t', header=False, index=False)
# else:
#     print("过滤后没有发现核序列。")

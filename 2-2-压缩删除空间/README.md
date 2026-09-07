# 2-2-压缩删除空间

清理 NUMTs 检测流程的两个输出目录，删除全部可再生的中间文件，并把 SAM 压缩为 BAM。

## 保留原则

### 步骤 1 —— `1-1-获得NUMTs分布/output`（服务器：`Result_NUMT/1-NUMTs/`）

`*_all_regions.bwa.sam` 是链条上**唯一同时含「序列 + read 名 + 比对坐标」**的文件，
其余下游文件全部可由它重算，因此它被压缩为 BAM 后保留，其余可再生文件删除。

| 处置 | 文件 | 依据 |
|---|---|---|
| 压缩为 BAM | `*_all_regions.bwa.sam` | 序列唯一载体；校验记录数一致后删除原 SAM |
| 保留 | `*.mt.disc.sam`、`*.mt.split.sam` | mt 相关 discordant/split read 的唯一存档 |
| 保留 | `*.breakpointINPUT.tsv`、`*.cluster.tsv`、`*.cluster.summary.tsv` | 区域清单与聚类证据 |
| 保留 | `*.AllBreakpoints.tsv`、`*.ConfidentBreakpoints.tsv` | 最终断点结果 |
| 删除 | `*_all_regions.fasta.cap.*` | CAP3 中间产物，占该目录约 95% 体积 |
| 删除 | `*_all_regions.fasta` | 可由 BAM 的 primary 记录无损还原（已实测逐条一致） |
| 删除 | `*_all_regions.psl`、`*.psl` | 可由 BAM 经 `0_2_sam2psl.py` 重算 |
| 删除 | `*.Breakpoints.old.tsv`、`*.cluster.old.tsv` | legacy 对照文件 |

从 BAM 还原 fasta：

```bash
samtools view sample_all_regions.bwa.bam \
  | awk 'BEGIN{OFS=""} {f=$2; if(and(f,256)||and(f,2048))next; s=$10;
      if(and(f,16)){r="";for(i=length(s);i>0;i--){c=substr(s,i,1);
        r=r (c=="A"?"T":c=="T"?"A":c=="C"?"G":c=="G"?"C":"N")}s=r}
      print ">",$1,"\n",s}' > restored.fasta
```

### 步骤 2 —— `1-5-Vardetection/output`（服务器：`Result_NUMT/1-5-Vardetect/output/`）

两个 `*aln.fa`（clustalo 多序列比对）是序列本体，配合 `*_numts.bed`
即可用 `generateVariantTable.Human.py` / `generateVariantTable.HumanChimp.py`
再生全部逐碱基 TSV，而这些明细表占该目录约 90% 体积。

| 处置 | 文件 |
|---|---|
| 保留 | `alnHuman/*.humanMTaln.fa`、`alnHumanChimp/*.humanchimpMTaln.fa` |
| 保留 | `*_numts.bed`、`*_all_regions.filtered.fasta` |
| 保留 | `*.numtVarFilterPos.tsv`、`*.numtDhumanChimp.tsv`、`*.numtDhumanChimp.sum.tsv` |
| 保留 | `pslHuman/*.psl`、`pslChimp/*.psl`、`logs/success.log` |
| 删除 | `*.humanMTaln.fa.numt.tsv`、`*.numtVar.tsv`、`*.full.tsv`、`*.humanchimpMTaln.fa.numt.tsv` |
| 删除 | `*_all_regions.humanMT.fasta`、`*_all_regions.humanchimpMT.fasta`（= `cat` 参考序列 + `filtered.fasta`） |
| 删除 | `*.fasta.cap*`、`fastaFiles.list` |

再生逐碱基 TSV：

```bash
python3 generateVariantTable.Human.py      alnHuman/X.humanMTaln.fa          X_numts.bed
python3 generateVariantTable.HumanChimp.py alnHumanChimp/X.humanchimpMTaln.fa X_numts.bed
```

## 安全设计

1. **keep 优先于 delete** —— 任何匹配保留清单的文件绝不会被删除，即使它同时匹配删除模式。
2. **未知文件原样保留** —— 既不匹配 keep 也不匹配 delete 的文件记为 `unknown` 并保留，
   在汇总报告中可见，便于发现模式清单未覆盖的新文件类型。
3. **默认 dry-run** —— `conf` 中 `runtime.dry_run: "true"`，只统计不改动。
4. **压缩先校验后删除** —— `samtools view -c` 比对 SAM/BAM 记录数，一致才删除原 SAM。
5. **先删后压** —— 单样本内先删除可再生文件腾出空间，再写 BAM，适应磁盘吃紧的场景。

## 用法

```bash
# 1) 生成 LSF 作业脚本（默认 dry-run，不提交），检查 temp/*/jobs/*.lsf
bash pipe/1-clean-numts-output.sh
bash pipe/2-clean-vardetect-output.sh

# 2) dry-run 全量统计，确认将删除的体积
bash pipe/1-clean-numts-output.sh --limit 0 --submit
bash pipe/1-clean-numts-output.sh --collect     # 作业结束后汇总

# 3) 小批量试删 50 个样本
bash pipe/1-clean-numts-output.sh --dry-run false --limit 50 --submit
bash pipe/1-clean-numts-output.sh --collect

# 4) 全量执行
bash pipe/1-clean-numts-output.sh --dry-run false --limit 0 --submit
bash pipe/1-clean-numts-output.sh --collect
```

`--limit 0` 表示全部样本。作业按 `runtime.samples_per_job` 切块，
每块一个 LSF 作业，块内用 GNU parallel 起 `runtime.jobs` 路并发。

## 产出

- `output/result/<步骤>/*-report.tsv` —— tidy long 长表，每行
  `sample_id / action / category / files / bytes`
- `output/result/<步骤>/*-summary.tsv` —— 按 `action × category` 聚合，含 GiB 列

`action` 取值：`kept` / `unknown` / `would_delete` / `deleted` / `delete_failed` /
`would_compress` / `compressed_from_sam` / `compressed_to_bam` /
`compress_failed` / `verify_failed` / `compress_skipped_exists`

## 目录

```
conf/    1-clean-numts-output.yaml, 2-clean-vardetect-output.yaml
pipe/    1-clean-numts-output.sh,   2-clean-vardetect-output.sh
script/  load_config.sh, clean_numts_sample.sh, clean_vardetect_sample.sh
temp/    <步骤>/{chunks,jobs,reports}/ 与模式清单
log/     LSF stdout/stderr
output/result/<步骤>/
```

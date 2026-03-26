# NUMTs 变异检测结果解读与下游分析建议

## 1. 主要输出目录与文件

每个样本输出位于：

```
/output/<SAMPLE_ID>/
```

常见子目录与文件：

- `assembly/`：CAP3 拼接相关产物（如 `.cap`、`.contigs`）。若 CAP3 不可用或 contigs 未生成，可能为空或只包含日志。
- `pslHuman/`：NUMT 序列对人类 mtDNA 参考的 BLAT 比对结果（`.human.psl`）。
- `pslChimp/`：NUMT 序列对黑猩猩 mtDNA 参考的 BLAT 比对结果（`.chimp.psl`）。
- `alnHuman/`：人类 mtDNA + NUMT 的多序列比对结果与派生表格。
- `alnHumanChimp/`：人类 mtDNA + 黑猩猩 mtDNA + NUMT 的多序列比对结果与派生表格。
- `<SAMPLE_ID>_numts.bed`：由 PSL 生成的 NUMT 区间注释表（4 列 TSV）。

---

## 2. 关键文件解读

### 2.1 `alnHuman/` 目录

#### 2.1.1 `<sample>.humanMTaln.fa`
- 内容：人类 mtDNA 参考 + 多条 NUMT 序列（多序列比对 FASTA）。
- 用途：用于 Human-only 变异判定。

#### 2.1.2 `<sample>.humanMTaln.fa.numt.tsv`
- 全量变异表（未过滤）。
- 每行对应比对中的一个位置。
- 主要字段：
  - `ref`：人类 mtDNA 参考碱基
  - `alt`：NUMT 碱基
  - `pos`：mtDNA 位置（1-based）
  - `numtID`：NUMT/contig ID（来自 FASTA header）
  - `numt_start/numt_end`：来自 `<SAMPLE_ID>_numts.bed` 的对应区间

#### 2.1.3 `<sample>.humanMTaln.fa.numtVar.tsv`
- 仅保留 `ref != alt` 的变异位点（去重后）。
- 适合统计变异数、位点类型分布等。

#### 2.1.4 `<sample>.humanMTaln.fa.numtVarFilterPos.tsv`
- 变异位点进一步限定在 NUMT 区间 `[numt_start, numt_end]` 内。
- 这是最常用的“可信 NUMT 变异集合”。

### 2.2 `alnHumanChimp/` 目录

#### 2.2.1 `<sample>.humanchimpMTaln.fa`
- 内容：人类 mtDNA + 黑猩猩 mtDNA + 多条 NUMT（多序列比对）。
- 用途：判断 NUMT 更接近人类还是黑猩猩，从而推断“年龄”或起源方向。

#### 2.2.2 `<sample>.humanchimpMTaln.fa.full.tsv`
- 全量比对表（Human/Chimp/NUMT 三列）。

#### 2.2.3 `<sample>.humanchimpMTaln.fa.numt.tsv`
- 限定在 NUMT 区间内的比对位点。

#### 2.2.4 `<sample>.humanchimpMTaln.fa.numtDhumanChimp.tsv`
- 只保留 Human 与 Chimp 不同的位点，再看 NUMT 更像哪一边。
- 字段 `category`：
  - `humanM`：NUMT 与人类相同
  - `chimpM`：NUMT 与黑猩猩相同
  - `none`：都不同

#### 2.2.5 `<sample>.humanchimpMTaln.fa.numtDhumanChimp.sum.tsv`
- 汇总表（每个 NUMT 的计数与比例），常用于估计 NUMT 相对“年龄”。
- 关键字段：
  - `humanN`：NUMT 更接近人类的位点数
  - `chimpN`：NUMT 更接近黑猩猩的位点数
  - `pct_human`：humanN / (humanN + chimpN)
  - `age`：脚本中给出的经验估计值（基于人/猩差异比例）

---

## 3. 典型分析思路

### 3.1 NUMT 变异负荷与分布
- 统计每个 NUMT 的变异数（`numtVar.tsv` 或 `numtVarFilterPos.tsv`）。
- 比较不同样本之间的变异负荷差异。
- 可进一步按 mtDNA 区段（D-loop、编码区等）分层统计。

### 3.2 NUMT 年龄推断
- 使用 `numtDhumanChimp.sum.tsv` 中的 `pct_human/age` 指标。
- 结合人群频率或谱系信息（如样本家系），可辅助判断插入发生的时间窗。

### 3.3 跨样本复现与共享事件
- 按 `numtID` 或基因组坐标对齐不同样本的 NUMT。
- 识别共有 NUMT（可能是古老插入）与样本特异 NUMT（可能是新近事件）。

### 3.4 NUMT 与核基因组注释关联
- 将 NUMT 坐标与基因注释、结构变异、重复序列区域叠加分析。
- 判断是否落在功能区或与结构变异相关。

### 3.5 质量控制与异常检查
- 若某些 NUMT 仅含 N 或极短序列，可能是上游组装质量问题。
- 可统计过滤前后 NUMT 数量，评估数据质量。

---

## 4. 实际使用建议

- 常用主结果文件：
  - Human：`*.humanMTaln.fa.numtVarFilterPos.tsv`
  - HumanChimp：`*.humanchimpMTaln.fa.numtDhumanChimp.sum.tsv`
- 若要做统计汇总，建议先把多个样本同名结果文件合并再分析。

---

## 5. 下一步可选扩展

- 以 NUMT 为单位，构建变异特征矩阵（样本 × 变异数）。
- 与人群样本或家系信息联合做聚类/树状图分析。
- 结合外部 mtDNA 变异数据库做功能或频率注释。
- 若有更高质量参考序列，可替换 `REF_HUMAN/REF_CHIMP` 重新跑。



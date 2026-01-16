# GRCh37 -> GRCh38 liftover for NUMTs

## 目录
- `conf/config.yaml`: 配置（liftOver 路径、chain、列名、并发数等）
- `script/1-liftover.sh`: 调度脚本（GNU parallel / xargs）
- `python/liftover_numts.py`: 单文件转换模块

## 使用方法
1. 安装 UCSC liftOver，并准备 `hg19ToHg38.over.chain.gz`。
2. 修改 `conf/config.yaml`：
   - `liftover.bin` 指向 liftOver 可执行文件
   - `liftover.chain` 指向 chain 文件
3. 在 `conf/config.yaml`、`conf/config_confident.yaml`、`conf/config_breakpointinput.yaml`、`conf/config_cluster.yaml` 或 `conf/config_cluster_summary.yaml` 中设置：
   - `input.path` 输入文件
   - `output.path` 输出文件
4. 运行：
   ```bash
   /bin/bash /mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/0-GRCH37-liftover/script/1-liftover.sh --config /mnt/f/Onedrive/文档（科研）/脚本/Download/15-NUMTs-detector-V2/0-GRCH37-liftover/conf/config.yaml
   ```

## 说明
- 仅对常染色体端做 liftover；`MT/chrM` 端保持不变（仅统一命名）。
- `output.on_fail` 可控制失败行处理：`keep` 保留原坐标、`drop` 删除行、`na` 置空坐标。
- 日志记录在 `log/1-liftover.sh.log`，断点续跑只依赖日志判断。

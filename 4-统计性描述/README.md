# 介绍

本文件夹负责对NUMTs的分布进行统计性描述。

其中 `distinct NUMTs` 的统计不再基于 `3-numt-events.tsv` 的事件级 1 kb 聚类，而是基于
`input.distinct_numt_cluster_input_tsv_gz` 指向的多样本原始输入，调用
`python/groupNumtCluster_fromMultipleSamples.py` 生成 `GroupID` 聚类结果。

`conf/Config.yaml` 中的 `analysis.distinct_numt_min_sample_supports` 用于控制
GroupID 至少需要在多少个样本中出现才计入 distinct NUMTs，当前支持一次性并行计算多个阈值。

# 一、查看完整输出日志（重点看是否有SUCCESS提示）
cat /project/liangzhanhao/m6a_analysis/log/star_batch_align.out

# 查看错误日志（无内容=无报错）
cat /project/liangzhanhao/m6a_analysis/log/star_batch_align.err

#二、验证输出目录的文件（核心！确认 BAM 是否生成）
# 1. 查看align_results下是否有6个样本的子目录
ls -l /project/liangzhanhao/m6a_analysis/align_results/

# 2. 查看单个样本的输出文件（以SRR29917562为例）
ls -l /project/liangzhanhao/m6a_analysis/align_results/SRR29917562/

#正常输出:
#SRR29917562_Aligned.sortedByCoord.out.bam  # 核心比对结果（大小至少几百MB）
#SRR29917562_Aligned.sortedByCoord.out.bam.bai  # BAM索引（可选）
#SRR29917562_Log.final.out  # STAR比对统计日志（关键！）
#SRR29917562_Log.out  # STAR运行日志

#三、验证 BAM 文件有效性（确保文件未损坏)
# 用samtools查看BAM的比对统计（以SRR29917562为例）
/project/liangzhanhao/soft/bin/samtools flagstat /project/liangzhanhao/m6a_analysis/align_results/SRR29917562/SRR29917562_Aligned.sortedByCoord.out.bam

#预期输出：
'''
10000000 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
9500000 + 0 mapped (95.00% : N/A)  # 比对率正常≥80%
...
'''

WORK_DIR="/project/liangzhanhao/m6a_analysis"

# 创建主目录及所有子目录（-p 确保父目录不存在时自动创建）
mkdir -p $WORK_DIR/{
  raw_srr,          # 存放原始SRR文件（从NCBI下载的.sra文件）
  raw_fastq,        # 从SRR转换后的原始fastq文件（未过滤），一般下载就是fastq
  clean_fastq,      # 质量过滤后的clean fastq文件（用于比对）
  qc_report/{raw,clean,summary},  # 质量控制报告：原始fastq的QC、过滤后QC、汇总报告
  reference/{genome,gtf,index},   # 参考基因组相关：基因组序列、注释文件、STAR/bowtie2索引。但是gtf和索引一般已经建立在公共库
  bam,              # 比对后的BAM文件及索引（.bam和.bai），是exomePeak2 输入
  result/{exomepeak2,enrichment,plots}  # 分析结果：差异位点、富集分析、可视化图表
}

#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --output=/project/liangzhanhao/m6a_analysis/log/star_align.out
#SBATCH --error=/project/liangzhanhao/m6a_analysis/log/star_align.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12  # 比对用12核，兼顾速度和资源
#SBATCH --mem=60G           # 比对内存需求低于索引构建，60G足够
#SBATCH --time=8:00:00      # 比对耗时根据数据量调整，8小时适配常规样本
#SBATCH --partition=h240    # 沿用能成功调度的分区

# ===================== 需替换的路径（核心！）=====================
# STAR绝对路径（和索引构建一致）
STAR_PATH="/project/liangzhanhao/soft/bin/STAR"
# 已构建好的STAR索引路径
GENOME_DIR="/project/liangzhanhao/m6a_analysis/reference/index/star_hg38/"
# 输入fastq文件（单端/双端二选一，注释掉不用的）
## 单端测序（示例）：
# FASTQ_FILE="/project/liangzhanhao/m6a_analysis/data/sample1.fastq.gz"
## 双端测序（示例）：
FASTQ_FILE1="/project/liangzhanhao/m6a_analysis/data/sample1_R1.fastq.gz"
FASTQ_FILE2="/project/liangzhanhao/m6a_analysis/data/sample1_R2.fastq.gz"
# 输出目录（自动创建，无需手动建）
OUT_DIR="/project/liangzhanhao/m6a_analysis/align_results/sample1"

# ===================== 固定流程（无需修改）=====================
# 创建输出目录
mkdir -p ${OUT_DIR}
echo "===== 开始比对：$(date) =====" >> ${OUT_DIR}/align.log

# 执行STAR比对（单端/双端二选一）
## 双端测序（默认启用，单端则注释这行，启用下面单端行）
${STAR_PATH} \
  --genomeDir ${GENOME_DIR} \
  --readFilesIn ${FASTQ_FILE1} ${FASTQ_FILE2} \
  --readFilesCommand zcat  # 若fastq未压缩，删除此行
  --outSAMtype BAM SortedByCoordinate  # 直接输出排序后的BAM（无需后续samtools排序）
  --outFilterMultimapNmax 1  # 只保留唯一比对的reads（MeRIP-seq常用）
  --outFilterType BySJout \
  --outFilterScoreMinOverLread 0.3 \
  --outFilterMatchNminOverLread 0.3 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \
  --quantMode TranscriptomeSAM  # 可选：输出转录组比对结果（用于定量）
  --outFileNamePrefix ${OUT_DIR}/sample1_ \
  --runThreadN ${SLURM_CPUS_PER_TASK}  # 调用申请的CPU核心

## 单端测序（启用需注释上面双端代码，取消下面注释）
# ${STAR_PATH} \
#   --genomeDir ${GENOME_DIR} \
#   --readFilesIn ${FASTQ_FILE} \
#   --readFilesCommand zcat \
#   --outSAMtype BAM SortedByCoordinate \
#   --outFilterMultimapNmax 1 \
#   --outFilterType BySJout \
#   --outFilterScoreMinOverLread 0.3 \
#   --outFilterMatchNminOverLread 0.3 \
#   --alignIntronMin 20 \
#   --alignIntronMax 1000000 \
#   --quantMode TranscriptomeSAM \
#   --outFileNamePrefix ${OUT_DIR}/sample1_ \
#   --runThreadN ${SLURM_CPUS_PER_TASK}

# 比对完成验证
if [ -f "${OUT_DIR}/sample1_Aligned.sortedByCoord.out.bam" ]; then
  echo "===== 比对成功！BAM文件生成：$(date) =====" >> ${OUT_DIR}/align.log
  # 生成BAM索引（后续IGV可视化/定量需要）
  samtools index ${OUT_DIR}/sample1_Aligned.sortedByCoord.out.bam
else
  echo "===== 比对失败！$(date) =====" >> ${OUT_DIR}/align.log
  exit 1
fi

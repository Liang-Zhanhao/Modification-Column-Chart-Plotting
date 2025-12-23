#!/bin/bash
#SBATCH --job-name=star_batch_align
#SBATCH --output=/project/liangzhanhao/m6a_analysis/log/star_batch_align.out
#SBATCH --error=/project/liangzhanhao/m6a_analysis/log/star_batch_align.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=8:00:00
#SBATCH --partition=hpxg

# ===================== 核心路径配置（和你成功的单个脚本一致） =====================
STAR_PATH="/project/liangzhanhao/soft/bin/STAR"
GENOME_DIR="/project/liangzhanhao/m6a_analysis/reference/index/star_hg38"
FASTQ_DIR="/project/liangzhanhao/m6a_analysis/raw_fasta"
OUT_ROOT="/project/liangzhanhao/m6a_analysis/align_results"
SAMPLE_LIST="/project/liangzhanhao/m6a_analysis/samples.list"  # 样本列表文件
SAMTOOLS="/project/liangzhanhao/soft/bin/samtools"  # 假设samtools路径和STAR同目录，若不同请替换

# ===================== 前置检查 =====================
# 创建日志目录
mkdir -p /project/liangzhanhao/m6a_analysis/log

# 检查STAR是否可执行
if [ ! -x "${STAR_PATH}" ]; then
    echo "ERROR: STAR路径不可执行 - ${STAR_PATH}" >&2
    exit 1
fi

# 检查基因组索引目录
if [ ! -d "${GENOME_DIR}" ] || [ ! -f "${GENOME_DIR}/Genome" ]; then
    echo "ERROR: STAR索引目录无效 - ${GENOME_DIR}" >&2
    exit 1
fi

# 检查样本列表文件
if [ ! -f "${SAMPLE_LIST}" ]; then
    echo "ERROR: 样本列表文件不存在 - ${SAMPLE_LIST}" >&2
    exit 1
fi

# 检查samtools
if [ ! -x "${SAMTOOLS}" ]; then
    echo "WARNING: samtools路径不可执行，跳过BAM索引 - ${SAMTOOLS}" >&2
    SAMTOOLS=""
fi

# ===================== 批量比对样本 =====================
echo "===== 开始批量STAR比对 | 时间：$(date) ====="
echo "样本列表：${SAMPLE_LIST}"
echo "输出目录：${OUT_ROOT}"

# 读取样本列表并循环处理
cat "${SAMPLE_LIST}" | while read -r sample; do
    # 跳过空行和注释行
    if [[ -z "${sample}" || "${sample}" =~ ^# ]]; then
        continue
    fi

    echo -e "\n===== 处理样本：${sample} | 时间：$(date) ====="
    
    # 定义当前样本路径
    FASTQ_FILE="${FASTQ_DIR}/${sample}.fastq.gz"
    OUT_DIR="${OUT_ROOT}/${sample}"
    BAM_FILE="${OUT_DIR}/${sample}_Aligned.sortedByCoord.out.bam"
    
    # 检查fastq文件是否存在
    if [ ! -f "${FASTQ_FILE}" ]; then
        echo "ERROR: 样本${sample}的fastq文件不存在 - ${FASTQ_FILE}" >&2
        continue
    fi

    # 创建输出目录
    mkdir -p "${OUT_DIR}"

    # STAR比对（和你成功的单个脚本参数完全一致）
    ${STAR_PATH} \
        --runMode alignReads \
        --genomeDir "${GENOME_DIR}" \
        --readFilesIn "${FASTQ_FILE}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${OUT_DIR}/${sample}_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within KeepPairs \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --runThreadN 16

    # 检查BAM是否生成
    if [ ! -f "${BAM_FILE}" ]; then
        echo "ERROR: 样本${sample}比对失败，未生成BAM文件 - ${BAM_FILE}" >&2
        continue
    fi

    echo "SUCCESS: 样本${sample}比对完成，BAM文件：${BAM_FILE}"

    # samtools构建BAM索引（可选，和你成功的脚本一致）
    if [ -n "${SAMTOOLS}" ]; then
        echo "开始构建${sample}的BAM索引..."
        ${SAMTOOLS} index "${BAM_FILE}"
        
        if [ -f "${BAM_FILE}.bai" ]; then
            echo "SUCCESS: ${sample}的BAM索引构建完成"
        else
            echo "ERROR: ${sample}的BAM索引构建失败" >&2
        fi
    fi

done

echo -e "\n===== 批量比对结束 | 时间：$(date) ====="
echo "所有结果已输出至：${OUT_ROOT}"

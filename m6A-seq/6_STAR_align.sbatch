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

# ===================== 核心修正：去掉子文件夹路径 =====================
STAR_PATH="/project/liangzhanhao/soft/bin/STAR"
GENOME_DIR="/project/liangzhanhao/m6a_analysis/reference/index/star_hg38"
FASTQ_DIR="/project/liangzhanhao/m6a_analysis/raw_fastq"  
OUT_ROOT="/project/liangzhanhao/m6a_analysis/align_results"
SAMPLE_LIST="/project/liangzhanhao/m6a_analysis/samples.list"
SAMTOOLS="/project/liangzhanhao/soft/bin/samtools"

# ===================== 前置检查 =====================
mkdir -p /project/liangzhanhao/m6a_analysis/log

if [ ! -x "${STAR_PATH}" ]; then
    echo "ERROR: STAR路径不可执行 - ${STAR_PATH}" >&2
    exit 1
fi
if [ ! -d "${GENOME_DIR}" ] || [ ! -f "${GENOME_DIR}/Genome" ]; then
    echo "ERROR: STAR索引目录无效 - ${GENOME_DIR}" >&2
    exit 1
fi
if [ ! -f "${SAMPLE_LIST}" ]; then
    echo "ERROR: 样本列表文件不存在 - ${SAMPLE_LIST}" >&2
    exit 1
fi
if [ ! -d "${FASTQ_DIR}" ]; then
    echo "ERROR: raw_fastq目录不存在 - ${FASTQ_DIR}" >&2
    exit 1
fi

# ===================== 批量比对（路径直接拼接：raw_fastq/样本名.fastq.gz） =====================
echo "===== 开始批量STAR比对 | 时间：$(date) ====="
echo "fastq根目录：${FASTQ_DIR}"
echo "样本列表：${SAMPLE_LIST}"
echo "输出目录：${OUT_ROOT}"

cat "${SAMPLE_LIST}" | while read -r sample; do
    if [[ -z "${sample}" || "${sample}" =~ ^# ]]; then
        continue
    fi

    echo -e "\n===== 处理样本：${sample} | 时间：$(date) ====="
    

    FASTQ_FILE="${FASTQ_DIR}/${sample}.fastq.gz"
    OUT_DIR="${OUT_ROOT}/${sample}"
    BAM_FILE="${OUT_DIR}/${sample}_Aligned.sortedByCoord.out.bam"
    
    # 打印实际查找的路径，方便核对
    echo "正在查找fastq文件：${FASTQ_FILE}"
    
    if [ ! -f "${FASTQ_FILE}" ]; then
        echo "ERROR: 样本${sample}的fastq文件不存在！" >&2
        # 自动检测后缀是否是.fq.gz（备选）
        FASTQ_FILE_FQ="${FASTQ_DIR}/${sample}.fq.gz"
        if [ -f "${FASTQ_FILE_FQ}" ]; then
            echo "提示：找到.fq.gz后缀的文件，自动切换路径：${FASTQ_FILE_FQ}" >&2
            FASTQ_FILE="${FASTQ_FILE_FQ}"
        else
            continue
        fi
    fi

    mkdir -p "${OUT_DIR}"

    # STAR比对（参数和你成功的脚本一致）
    echo "开始STAR比对..."
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

    if [ ! -f "${BAM_FILE}" ]; then
        echo "ERROR: 样本${sample}比对失败，未生成BAM文件" >&2
        continue
    fi

    echo "SUCCESS: 样本${sample}比对完成！BAM文件：${BAM_FILE}"

    # samtools索引
    if [ -x "${SAMTOOLS}" ]; then
        echo "开始构建BAM索引..."
        ${SAMTOOLS} index "${BAM_FILE}"
        if [ -f "${BAM_FILE}.bai" ]; then
            echo "SUCCESS: BAM索引构建完成！"
        fi
    fi

done

echo -e "\n===== 批量比对结束 | 时间：$(date) ====="
echo "所有结果输出至：${OUT_ROOT}"

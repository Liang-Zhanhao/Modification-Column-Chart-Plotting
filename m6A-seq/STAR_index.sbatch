#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --output=/project/liangzhanhao/m6a_analysis/log/star_index.out
#SBATCH --error=/project/liangzhanhao/m6a_analysis/log/star_index.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH --partition=hpxg

# 替换为正确的STAR绝对路径
STAR_PATH="/project/liangzhanhao/soft/bin/STAR"

# 1. 验证STAR路径是否存在（写入日志，方便排查）
echo "===== 脚本启动：$(date) =====" >> /project/liangzhanhao/m6a_analysis/log/star_index.out
echo "STAR绝对路径：$STAR_PATH" >> /project/liangzhanhao/m6a_analysis/log/star_index.out
if [ -f "$STAR_PATH" ]; then
  echo "STAR文件存在！版本信息：$($STAR_PATH --version)" >> /project/liangzhanhao/m6a_analysis/log/star_index.out
else
  echo "ERROR：STAR文件不存在！路径错误！" >> /project/liangzhanhao/m6a_analysis/log/star_index.out
  exit 1  # 路径错误直接终止脚本
fi

# 2. 创建索引输出目录（确保可写）
mkdir -p /project/liangzhanhao/m6a_analysis/reference/index/star_hg38/
echo "索引输出目录：/project/liangzhanhao/m6a_analysis/reference/index/star_hg38/" >> /project/liangzhanhao/m6a_analysis/log/star_index.out

# 3. 执行STAR索引构建（输出所有日志到文件）
echo "开始构建STAR索引：$(date)" >> /project/liangzhanhao/m6a_analysis/log/star_index.out
$STAR_PATH --runMode genomeGenerate \
  --genomeDir /project/liangzhanhao/m6a_analysis/reference/index/star_hg38/ \
  --genomeFastaFiles /project/liangzhanhao/00.DATABASE/hg38/hg38.fa \
  --sjdbGTFfile /project/liangzhanhao/00.DATABASE/hg38/gencode.v49.annotation.gtf \
  --sjdbOverhang 100 \    #测序reads长度-1
  --runThreadN 16 >> /project/liangzhanhao/m6a_analysis/log/star_index.out 2>&1

# 4. 执行完成后输出状态
echo "索引构建完成/终止：$(date)，退出码：$?" >> /project/liangzhanhao/m6a_analysis/log/star_index.out
echo "===== 脚本结束 =====" >> /project/liangzhanhao/m6a_analysis/log/star_index.out

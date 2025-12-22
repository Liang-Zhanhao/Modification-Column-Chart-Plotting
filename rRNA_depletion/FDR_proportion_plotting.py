import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# ----------------------
# 全局参数配置（非路径参数）
# ----------------------
file_suffix = ".csv"  # csv格式
sep = "\t"  # csv分隔符\t或者，
INCLUDE_UNKNOWN_TYPES = False  # 是否包含未定义转换类型的基因

# 关键：根据你的表头映射列名
column_mapping = {
    "Sites": "site",
    "Gene_Type": "rna_type"
}

# 基因类型→RNA类型的转换规则
gene_type_to_rna = {
    "protein_coding": "mRNA",
    "intergenic": "ncRNA",
    "pseudogene": "unaligned",
    "tRNA": "ncRNA",
    "rRNA_gene": "rRNA"
}

# 颜色映射
color_map = {
    "mRNA": "#E41A1C",
    "rRNA": "#377EB8",
    "ncRNA": "#FFD92F",
    "unaligned": "#999999"
}


# ----------------------
# 数据加载函数
# ----------------------
def load_samples(folder, group_name):
    if not folder.exists():
        raise FileNotFoundError(f"文件夹不存在：{folder}")
    if not folder.is_dir():
        raise NotADirectoryError(f"不是文件夹：{folder}")

    files = list(folder.glob(f"*{file_suffix}"))
    if not files:
        all_files = [f.name for f in folder.glob("*")]
        raise ValueError(f"在 {folder} 中未找到 {file_suffix} 文件！文件夹内文件：{all_files}")

    print(f"\n在 {group_name} 文件夹中找到 {len(files)} 个文件：")
    for f in files:
        print(f"  - {f.name}")

    all_samples = []
    for file in files:
        sample_name = file.stem
        try:
            df = pd.read_csv(file, sep=sep)
        except Exception as e:
            raise RuntimeError(f"读取文件 {file} 失败：{str(e)}")

        # 检查列名
        missing_cols = [col for col in column_mapping.keys() if col not in df.columns]
        if missing_cols:
            raise ValueError(f"文件 {file} 缺少表头中的列：{missing_cols}")

        df_renamed = df.rename(columns=column_mapping)
        original_count = len(df_renamed)

        # 过滤未定义类型
        unknown_types = set(df_renamed['rna_type']) - set(gene_type_to_rna.keys())
        if unknown_types and not INCLUDE_UNKNOWN_TYPES:
            df_filtered = df_renamed[df_renamed['rna_type'].isin(gene_type_to_rna.keys())]
            filtered_count = original_count - len(df_filtered)
            print(f"  文件 {file.name}：过滤掉 {filtered_count} 条未定义类型（{unknown_types}）的数据")
            df_renamed = df_filtered
        elif unknown_types:
            print(f"  警告：文件 {file.name} 中存在未定义转换规则的基因类型：{unknown_types}")

        # 转换RNA类型
        df_renamed['rna_type'] = df_renamed['rna_type'].map(
            lambda x: gene_type_to_rna.get(x, x)
        )

        df_renamed["sample"] = f"{group_name}_{sample_name}"
        df_renamed["group"] = group_name
        all_samples.append(df_renamed)

    combined_df = pd.concat(all_samples, ignore_index=True)
    print(f"  {group_name} 组最终保留数据量：{len(combined_df)} 条")
    return combined_df


# ----------------------
# 绘图函数
# ----------------------
def plot_rna_distribution(all_data, output_file="rRNA_depletion_effect_final.png"):
    # 统计占比
    counts = all_data.groupby(["sample", "rna_type"]).size().reset_index(name="count")
    total_counts = counts.groupby("sample")["count"].transform("sum")
    counts["percentage"] = (counts["count"] / total_counts) * 100
    pivot_data = counts.pivot(index="sample", columns="rna_type", values="percentage").fillna(0)

    # 补充未定义颜色
    for rna_type in pivot_data.columns:
        if rna_type not in color_map:
            color_map[rna_type] = "#AAAAAA"
            print(f"提示：RNA类型 {rna_type} 使用默认颜色")

    # 绘图设置
    plt.rcParams["font.family"] = ["Arial"]
    plt.rcParams['axes.unicode_minus'] = False

    fig, ax = plt.subplots(figsize=(10, 4))  # 画布大小
    bar_width = 0.6  # 柱子宽度
    bottom = np.zeros(len(pivot_data))

    # 堆叠柱状图
    for rna in pivot_data.columns:
        ax.bar(pivot_data.index, pivot_data[rna], width=bar_width, bottom=bottom,
               color=color_map[rna], label=rna)
        bottom += pivot_data[rna].values

    # 图形配置
    #隐藏XY轴线=False
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)

    # 调整Y轴刻度（仅保留标签，隐藏轴线）
    ax.set_yticks(np.arange(0, 101, 25))
    ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'], fontsize=7)
    #ax.tick_params(axis='y', length=0)  # 隐藏Y轴刻度线

    # 调整X轴标签（旋转且无轴线）
    ax.set_xticks(np.arange(len(pivot_data)))
    ax.set_xticklabels(pivot_data.index, fontsize=7, rotation=45, ha="right")
    #ax.tick_params(axis='x', length=0)  # 隐藏X轴刻度线

    # xy轴标题
    ax.set_title("rRNA depletion methods", fontsize=7)
    ax.set_xlabel('', fontsize=7)
    ax.set_ylabel('Relative Proportion (%)', fontsize=7)

    # 图例
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)

    # 坐标轴范围
    ax.set_xlim(-0.5, len(pivot_data) - 0.5)
    ax.set_ylim(0, 100)  # Y轴0-100

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n图形已保存至：{output_file}")
    plt.show()


# ----------------------
# 主函数
# ----------------------
def main(total_rna_dir, remove_rna_dir):
    try:
        print("开始读取total RNA样本...")
        total_rna_data = load_samples(total_rna_dir, "total")

        print("\n开始读取去rRNA样本...")
        remove_rna_data = load_samples(remove_rna_dir, "removed")

        all_data = pd.concat([total_rna_data, remove_rna_data], ignore_index=True)
        print(f"\n数据读取完成，共包含 {len(all_data)} 条记录")

        # 查看转换后的RNA类型分布
        print("\n转换后的RNA类型分布：")
        print(all_data['rna_type'].value_counts())

        # 绘图
        plot_rna_distribution(all_data)

    except Exception as e:
        print(f"\n运行出错：{str(e)}")


# ----------------------
# 执行入口（在这里修改文件路径）
# ----------------------
if __name__ == "__main__":
    # ========== 在这里修改你的文件路径 ==========
    TOTAL_RNA_DIR = Path("E:/private work/20251123_project2/rRNA_analysis/totalRNA")
    REMOVE_RNA_DIR = Path("E:/private work/20251123_project2/rRNA_analysis/remove_rRNA")

    # 调用主函数
    main(TOTAL_RNA_DIR, REMOVE_RNA_DIR)

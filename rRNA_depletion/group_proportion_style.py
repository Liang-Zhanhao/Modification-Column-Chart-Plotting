import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import warnings

# 忽略非致命警告（如libpng警告）
warnings.filterwarnings('ignore')

# ----------------------
# 全局参数配置（非路径参数）
# ----------------------
file_suffix = ".csv"  # csv格式
sep = "\t"  # csv分隔符\t或者，
INCLUDE_UNKNOWN_TYPES = False  # 是否包含未定义转换类型的基因

# 关键：根据你的表头映射列名（修改：将映射后的列名改为gene_type，避免命名混淆）
column_mapping = {
    "Sites": "site",
    "Gene_Type": "gene_type"  # 原："rna_type"，容易和转换后的RNA类型混淆
}

# 生信分析常用 基因类型→RNA类型 完整转换字典（统一小写键，方便大小写兼容）
# 恢复关键的ncRNA类型注释，确保匹配范围完整
gene_type_to_rna = {
    # 一、编码RNA（Protein-coding RNA）
    "protein_coding": "mRNA",

    # 二、非编码RNA（ncRNA）- 转运RNA相关
    "trna": "ncRNA",
    "trna_pseudogene": "ncRNA",

    # 三、非编码RNA（ncRNA）- 核糖体RNA相关
    "rrna_gene": "rRNA",
    "rrna_pseudogene": "ncRNA",

    # 四、非编码RNA（ncRNA）- 微小RNA相关
    "mirna": "ncRNA",
    "mirna_pseudogene": "ncRNA",

    # 五、非编码RNA（ncRNA）- 小核RNA相关
    "snrna": "ncRNA",
    "snrna_pseudogene": "ncRNA",

    # 六、非编码RNA（ncRNA）- 小核仁RNA相关
    "snorna": "ncRNA",
    "snorna_pseudogene": "ncRNA",

    # 七、非编码RNA（ncRNA）- 长链非编码RNA相关
    "lncrna": "ncRNA",
    "lincrna": "ncRNA",  # 基因间长链非编码RNA（intergenic相关核心类型）
    "lncrna_pseudogene": "ncRNA",

    # 八、非编码RNA（ncRNA）- 其他常见亚型
    "pirna": "ncRNA",
    "scrna": "ncRNA",
    "circrna": "ncRNA",
    "sirna": "ncRNA",
    "vaultrna": "ncRNA",
    "yrna": "ncRNA",
    "scarna": "ncRNA",
    "7skrna": "ncRNA",
    "rnase_p_rna": "ncRNA",
    "rnase_mrp_rna": "ncRNA",


    # 九、假基因（恢复注释，假基因属于ncRNA）
 #   "pseudogene": "ncRNA",
 #   "processed_pseudogene": "ncRNA",
 #   "unprocessed_pseudogene": "ncRNA",
 #   "polymorphic_pseudogene": "ncRNA",

    # 十、基因间（intergenic）相关RNA（恢复注释，核心ncRNA类型）
  #  "intergenic_rna": "ncRNA",
  #  "intergenic_ncrna": "ncRNA",
  #  "long_intergenic_rna": "ncRNA",
  #  "igrna": "ncRNA",

    # 十一、其他特殊RNA类型（恢复注释，补充ncRNA类型）
  #  "antisense_rna": "ncRNA",
  #  "sense_intronic_rna": "ncRNA",
 #   "sense_overlapping_rna": "ncRNA",
 #   "non_coding": "ncRNA"
}


def convert_gene_to_rna(gene_type, default="unknown"):
    """
    基因类型→RNA类型转换函数（大小写不敏感）
    参数：
        gene_type: 输入的基因类型字符串（支持任意大小写，如"LincRNA"、"MIRNA"）
        default: 未匹配到时的默认返回值（默认"unknown"）
    返回：
        对应的RNA类型（"mRNA"、"ncRNA"、"rRNA"或default）
    """
    # 新增：处理空值（如果基因类型为空，直接返回默认值）
    if pd.isna(gene_type) or gene_type.strip() == "":
        return default
    # 核心：将输入转为小写，再与字典匹配（解决大小写差异问题）
    lower_gene_type = gene_type.strip().lower()
    return gene_type_to_rna.get(lower_gene_type, default)


# 颜色映射（保持不变）
color_map = {
    "mRNA": "#E41A1C",
    "rRNA": "#377EB8",
    "ncRNA": "#FFD92F",
    "unknown": "#999999",  # 修改：将"unaligned"改为"unknown"，与函数默认值一致
    "": "#999999"  # 处理空值颜色
}

# ----------------------
# 分组规则配置（保持不变）
# ----------------------
grouping_rules = [
    ("removed", ["m6A_03"], "OD600=0.3(-rRNA)"),  # 第一组
    ("removed", ["m6A_1"], "OD600=1(-rRNA)"),  # 第二组
    ("total", ["m6A_1", "m6A_2", "m6A_3", "m6A_4"], "OD600=0.3"),  # 第三组
    ("total", ["m6A_5", "m6A_6", "m6A_7", "m6A_8"], "OD600=1")  # 第四组
]


# ----------------------
# 辅助函数：根据文件名和文件夹类型获取组名（保持不变）
# ----------------------
def get_group_name(file_name, folder_type):
    """
    根据文件名前缀和文件夹类型获取目标组名
    folder_type: "total" 或 "removed"
    """
    for rule_folder, prefixes, target_group in grouping_rules:
        if rule_folder == folder_type:
            if any(file_name.startswith(prefix) for prefix in prefixes):
                return target_group
    raise ValueError(f"文件 {file_name} 无法匹配任何分组规则（文件夹类型：{folder_type}）")


# ----------------------
# 数据加载函数（关键修改：使用convert_gene_to_rna函数）
# ----------------------
def load_samples(folder, folder_type):
    if not folder.exists():
        raise FileNotFoundError(f"文件夹不存在：{folder}")
    if not folder.is_dir():
        raise NotADirectoryError(f"不是文件夹：{folder}")

    files = list(folder.glob(f"*{file_suffix}"))
    if not files:
        all_files = [f.name for f in folder.glob("*")]
        raise ValueError(f"在 {folder} 中未找到 {file_suffix} 文件！文件夹内文件：{all_files}")

    print(f"\n在 {folder} 文件夹中找到 {len(files)} 个文件：")
    for f in files:
        print(f"  - {f.name}")

    group_data = {}

    for file in files:
        file_name = file.stem
        try:
            target_group = get_group_name(file_name, folder_type)
        except ValueError as e:
            print(f"警告：{e}，该文件将被跳过")
            continue

        try:
            df = pd.read_csv(file, sep=sep)
        except Exception as e:
            raise RuntimeError(f"读取文件 {file} 失败：{str(e)}")

        missing_cols = [col for col in column_mapping.keys() if col not in df.columns]
        if missing_cols:
            raise ValueError(f"文件 {file} 缺少表头中的列：{missing_cols}")

        df_renamed = df.rename(columns=column_mapping)
        original_count = len(df_renamed)

        # 关键修改1：使用转换函数预处理基因类型，获取所有可能的转换结果（含大小写兼容）
        df_renamed['converted_rna_type'] = df_renamed['gene_type'].apply(convert_gene_to_rna)

        # 关键修改2：基于转换结果过滤，而不是原始基因类型
        unknown_types = set(df_renamed['converted_rna_type']) - {"mRNA", "ncRNA", "rRNA"}
        if unknown_types and not INCLUDE_UNKNOWN_TYPES:
            df_filtered = df_renamed[df_renamed['converted_rna_type'].isin({"mRNA", "ncRNA", "rRNA"})]
            filtered_count = original_count - len(df_filtered)
            print(f"  文件 {file.name}：过滤掉 {filtered_count} 条未定义类型（{unknown_types}）的数据")
            df_renamed = df_filtered
        elif unknown_types:
            print(f"  警告：文件 {file.name} 中存在未定义转换规则的基因类型：{unknown_types}")

        # 关键修改3：最终使用转换后的RNA类型（大小写兼容）
        df_renamed["rna_type"] = df_renamed["converted_rna_type"]
        df_renamed = df_renamed.drop(columns=["converted_rna_type"])  # 删除中间列
        df_renamed["group"] = target_group

        if target_group not in group_data:
            group_data[target_group] = []
        group_data[target_group].append(df_renamed)

    all_samples = []
    for group_name, dfs in group_data.items():
        combined_group_df = pd.concat(dfs, ignore_index=True)
        print(f"  组 {group_name}：合并 {len(dfs)} 个文件，保留 {len(combined_group_df)} 条数据")
        all_samples.append(combined_group_df)

    if not all_samples:
        raise ValueError(f"{folder} 文件夹中没有符合分组规则的有效文件")

    final_df = pd.concat(all_samples, ignore_index=True)
    print(f"  {folder_type} 类型文件夹最终总数据量：{len(final_df)} 条")
    return final_df


# ----------------------
# 绘图函数（保持不变，仅适配颜色映射）
# ----------------------
def plot_rna_distribution(all_data, output_file="rRNA_depletion_effect_final.png"):
    # 按组统计各类RNA的数量和占比
    counts = all_data.groupby(["group", "rna_type"]).size().reset_index(name="count")
    total_counts = counts.groupby("group")["count"].transform("sum")
    counts["percentage"] = (counts["count"] / total_counts) * 100

    # 透视表整理数据
    pivot_data = counts.pivot(index="group", columns="rna_type", values="percentage").fillna(0)

    # 补充未定义颜色
    for rna_type in pivot_data.columns:
        if rna_type not in color_map:
            color_map[rna_type] = "#AAAAAA"
            print(f"提示：RNA类型 {rna_type} 使用默认颜色 #AAAAAA")

    # 新增：文字输出每个柱子每一层的比例
    print("\n" + "=" * 60)
    print("每个组（柱子）的RNA类型占比明细（保留2位小数）：")
    print("=" * 60)

    # 按指定顺序输出各组
    group_order = [
        "OD600=0.3(-rRNA)",
        "OD600=1(-rRNA)",
        "OD600=0.3",
        "OD600=1"
    ]

    for group in group_order:
        if group in pivot_data.index:
            print(f"\n【{group}】")
            print(f"  总数据量：{int(total_counts[counts['group'] == group].iloc[0])} 条")
            # 遍历该组的所有RNA类型，输出占比
            for rna_type in pivot_data.columns:
                percentage = pivot_data.loc[group, rna_type]
                print(f"  - {rna_type}：{percentage:.2f}%")
            # 验证总和是否为100%（避免计算误差）
            total_percent = pivot_data.loc[group].sum()
            print(f"  占比总和：{total_percent:.2f}%（理论应为100%）")

    print("=" * 60 + "\n")

    # 绘图设置
    plt.rcParams["font.family"] = ["Arial"]
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['text.usetex'] = False

    # 确保组的显示顺序正确
    pivot_data = pivot_data.reindex(group_order)

    fig, ax = plt.subplots(figsize=(5.5, 6))
    bar_width = 0.6
    bottom = np.zeros(len(pivot_data))

    # 绘制堆叠柱状图
    for rna_type in pivot_data.columns:
        ax.bar(
            pivot_data.index,
            pivot_data[rna_type],
            width=bar_width,
            bottom=bottom,
            color=color_map[rna_type],
            label=rna_type,
            edgecolor='white',
            linewidth=0.5
        )
        bottom += pivot_data[rna_type].values

    # 图形配置
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)

    # 调整Y轴
    ax.set_yticks(np.arange(0, 101, 20))
    ax.set_yticklabels(['0%', '20%', '40%', '60%', '80%', '100%'], fontsize=10)
    ax.set_ylim(0, 100)

    # 调整X轴
    ax.set_xticks(np.arange(len(pivot_data)))
    ax.set_xticklabels(pivot_data.index, fontsize=10, rotation=45, ha="right")
    ax.tick_params(axis='x', pad=8)

    # 轴标题
    ax.set_title("RNA type", fontsize=14, pad=15)
    ax.set_xlabel('Group', fontsize=12)
    ax.set_ylabel('Proportion (%)', fontsize=12, labelpad=10)

    # 图例
    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        fontsize=10,
        frameon=False
    )

    # 调整布局并保存
    plt.tight_layout()
    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches='tight',
        facecolor='white',
        edgecolor='none'
    )
    print(f"\n图形已保存至：{output_file}")
    plt.show()


# ----------------------
# 主函数（保持不变）
# ----------------------
def main(total_rna_dir, remove_rna_dir):
    try:
        print("开始读取total RNA样本（对应分组规则中的 'total' 类型）...")
        total_rna_data = load_samples(total_rna_dir, folder_type="total")

        print("\n开始读取去rRNA样本（对应分组规则中的 'removed' 类型）...")
        remove_rna_data = load_samples(remove_rna_dir, folder_type="removed")

        all_data = pd.concat([total_rna_data, remove_rna_data], ignore_index=True)
        print(f"\n数据读取完成，共包含 {len(all_data)} 条记录")

        print("\n各组数据量分布：")
        group_counts = all_data['group'].value_counts()
        for group, count in group_counts.items():
            print(f"  {group}：{count} 条记录")

        print("\n转换后的RNA类型分布（关键：查看ncRNA是否存在）：")
        rna_distribution = all_data['rna_type'].value_counts()
        for rna_type, count in rna_distribution.items():
            print(f"  {rna_type}：{count} 条记录（{count / len(all_data) * 100:.2f}%）")

        # 绘图（会自动输出比例明细）
        plot_rna_distribution(all_data)

    except Exception as e:
        print(f"\n运行出错：{str(e)}")


# ----------------------
# 执行入口（保持不变）
# ----------------------
if __name__ == "__main__":
    TOTAL_RNA_DIR = Path("E:/private work/20251123_project2/rRNA_analysis/totalRNA")
    REMOVE_RNA_DIR = Path("E:/private work/20251123_project2/rRNA_analysis/remove_rRNA")

    main(TOTAL_RNA_DIR, REMOVE_RNA_DIR)

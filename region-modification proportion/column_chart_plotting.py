import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

# 设置中文字体
plt.rcParams["font.family"] = ["Arial"]
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题


def parse_region(gene_type, gene_name):
    """解析位点区域（根据Gene_Type和Gene_Name）"""
    gene_type = str(gene_type).lower()
    gene_name = str(gene_name).lower()

    if 'upstream' in gene_name:
        return "5' UTR"
    elif 'downstream' in gene_name:
        return "3' UTR"

    if gene_type in ['trna', 'rrna', 'ncrna', 'non-coding', 'pseudogene']:
        return "ncRNA"
    elif gene_type == 'protein_coding':
        return "CDS"
    elif gene_type == 'intergenic':
        return "Intergenic"
    else:
        return "Other"


def process_sites_file(file_path):
    """处理位点文件，返回：
    - site_regions: {位点ID: 区域}
    - site_to_rna: {位点ID: 所属RNAID}（关键：关联位点到RNA分子）
    """
    df = pd.read_csv(file_path)
    site_regions = {}
    site_to_rna = {}  # 新增：位点→RNA的映射（实际需从文件中提取RNAID）

    gene_type_cols = [col for col in df.columns if col.endswith('_Gene_Type')]
    gene_name_cols = [col for col in df.columns if col.endswith('_Gene_Name')]
    rna_id_col = 'RNA_ID' if 'RNA_ID' in df.columns else df.columns[2]  # 假设RNAID列存在

    for _, row in df.iterrows():
        site_id = f"{row['Chr']}:{row['Sites']}"
        # 提取RNAID（实际数据中需替换为真实RNAID列）
        rna_id = row[rna_id_col] if pd.notna(row[rna_id_col]) else site_id.split(':')[0] + '_RNA' + str(row.name)

        gene_types = [row[col] for col in gene_type_cols if pd.notna(row[col])]
        gene_names = [row[col] for col in gene_name_cols if pd.notna(row[col])]

        if gene_types and gene_names:
            region = parse_region(gene_types[0], gene_names[0])
            site_regions[site_id] = region
            site_to_rna[site_id] = rna_id  # 关联位点到RNA

    return site_regions, site_to_rna


def main(file1, file2, output_file='modification_region_distribution_correct_biology.png'):
    # 处理文件，获取位点-区域、位点-RNA映射
    sites_37_regions, sites_37_rna = process_sites_file(file1)
    sites_45_regions, sites_45_rna = process_sites_file(file2)

    # 定义区域顺序
    region_order = ["ncRNA", "3' UTR", "Intergenic", "5' UTR", "CDS"]

    # -------------------------- 第一步：鉴定RNA的共有/特有 --------------------------
    # 提取两组的RNA集合
    rna_37 = set(sites_37_rna.values())  # 37°C组表达的RNA
    rna_45 = set(sites_45_rna.values())  # 45°C组表达的RNA

    # 分类RNA：共有RNA、37°C特有RNA、45°C特有RNA
    shared_rna = rna_37 & rna_45  # 两组均表达的RNA
    unique_rna_37 = rna_37 - rna_45  # 仅37°C表达的RNA
    unique_rna_45 = rna_45 - rna_37  # 仅45°C表达的RNA

    # -------------------------- 第二步：按生物逻辑统计各层位点 --------------------------
    # 初始化各层计数
    layer_counts = {
        'mod_unique37': [0] * len(region_order),  # 37°C特有RNA上的修饰位点
        'unique_37': [0] * len(region_order),  # 共有RNA上仅37°C有的修饰位点
        'shared': [0] * len(region_order),  # 共有RNA上的共有修饰位点
        'unique_45': [0] * len(region_order),  # 共有RNA上仅45°C有的修饰位点
        'mod_unique45': [0] * len(region_order)  # 45°C特有RNA上的修饰位点
    }

    # 1. 统计共有RNA上的修饰位点（区分共有/特有修饰）
    shared_site_37 = {s for s in sites_37_regions if sites_37_rna[s] in shared_rna}
    shared_site_45 = {s for s in sites_45_regions if sites_45_rna[s] in shared_rna}

    # 共有修饰位点（同一RNA+同位点在两组均出现）
    shared_mod_sites = set()
    for s37 in shared_site_37:
        chr_pos = s37.split(':')[-1]
        for s45 in shared_site_45:
            if s45.split(':')[-1] == chr_pos and sites_37_regions[s37] == sites_45_regions[s45]:
                shared_mod_sites.add(s37)
                shared_mod_sites.add(s45)

    # 2. 逐区域统计各层
    for idx, region in enumerate(region_order):
        # a. Mod_unique45：45°C特有RNA上的修饰位点
        mod_unique45_sites = [s for s in sites_45_regions
                              if sites_45_regions[s] == region and sites_45_rna[s] in unique_rna_45]
        layer_counts['mod_unique45'][idx] = len(mod_unique45_sites)

        # b. Unique_45：共有RNA上仅45°C有的修饰位点
        unique_45_sites = [s for s in shared_site_45
                           if sites_45_regions[s] == region and s not in shared_mod_sites]
        layer_counts['unique_45'][idx] = len(unique_45_sites)

        # c. Shared：共有RNA上的共有修饰位点
        shared_sites = [s for s in shared_mod_sites if sites_37_regions.get(s, '') == region]
        layer_counts['shared'][idx] = len(shared_sites)

        # d. Unique_37：共有RNA上仅37°C有的修饰位点
        unique_37_sites = [s for s in shared_site_37
                           if sites_37_regions[s] == region and s not in shared_mod_sites]
        layer_counts['unique_37'][idx] = len(unique_37_sites)

        # e. Mod_unique37：37°C特有RNA上的修饰位点
        mod_unique37_sites = [s for s in sites_37_regions
                              if sites_37_regions[s] == region and sites_37_rna[s] in unique_rna_37]
        layer_counts['mod_unique37'][idx] = len(mod_unique37_sites)

    # -------------------------- 第三步：归一化到Y轴100（保持比例） --------------------------
    max_height = 100
    normalized = {}
    for layer in layer_counts:
        normalized[layer] = []
        for idx in range(len(region_order)):
            total = sum([layer_counts[l][idx] for l in layer_counts])
            if total == 0:
                normalized[layer].append(0)
            else:
                normalized[layer].append((layer_counts[layer][idx] / total) * max_height)

    # -------------------------- 第四步：正确堆叠绘制（生物逻辑顺序） --------------------------
    colors = {
        'Unique modification (OD1)': '#aec7e8',  # 共有RNA上仅37°C的修饰
        'Shared modification': '#c7c7c7',  # 共有RNA上的共有修饰
        'Unique modification (OD0.3)': '#ffbb78',  # 共有RNA上仅45°C的修饰
    }

    fig, ax = plt.subplots(figsize=(6.5, 5.4))
    bar_width = 0.9
    index = np.arange(len(region_order))
    bottom = np.zeros(len(region_order))

    # 严格按生物逻辑堆叠：从下到上→45°C特有RNA→共有RNA特有修饰→共有修饰→37°C特有修饰→37°C特有RNA
    layers_draw_order = ['unique_45', 'shared', 'unique_37']
    labels = {
        'unique_45': 'Unique modification (OD0.3)',
        'shared': 'Shared modification',
        'unique_37': 'Unique modification (OD1)',
    }

    for layer in layers_draw_order:
        # 添加黑色边框：edgecolor='black', linewidth=0.8
        ax.bar(index, normalized[layer], bar_width, bottom=bottom,
               color=colors[labels[layer]], label=labels[layer],
               edgecolor='black', linewidth=0.8)  # 关键修改：添加黑色边框
        bottom += np.array(normalized[layer])  # 正确更新bottom

    # -------------------------- 图形配置 --------------------------
    # 坐标轴标签：设为7号
    ax.set_xlabel('Region', fontsize=7)
    ax.set_ylabel('Relative Proportion (%)', fontsize=7)

    # 坐标轴刻度标签：设为7号
    ax.set_xticks(index)
    ax.set_xticklabels(region_order, fontsize=7)
    ax.tick_params(axis='y', labelsize=7)  # Y轴刻度标签设为7号

    # 图例字体：设为7号
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)

    # 调整x轴范围
    ax.set_xlim(-0.5, len(region_order)-0.5)

    # Y轴固定100
    ax.set_ylim(0, max_height)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"图形已保存至：{output_file}")
    plt.show()


if __name__ == "__main__":
    file_37 = "OD_1_0.75_sites.txt"
    file_45 = "OD_03_0.75_sites.txt"
    main(file_37, file_45)

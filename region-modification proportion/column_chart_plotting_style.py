import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict


# -------------------------- 独立的风格配置函数 --------------------------
def set_nature_methods_style():
    """设置Nature Methods期刊风格的matplotlib参数"""
    plt.style.use('default')
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.sans-serif': 'Arial',
        'axes.unicode_minus': False,
        'font.size': 7,
        'axes.titlesize': 7,
        'axes.labelsize': 7,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 7,
        'axes.linewidth': 0.5,
        'grid.linewidth': 0.3,
        'lines.linewidth': 1.5,
        'legend.frameon': False,
        'pdf.fonttype': 42,  # 关键：PDF字体嵌入，避免乱码
        'figure.dpi': 600,
        'savefig.dpi': 600
    })


# -------------------------- 功能函数 --------------------------
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
    """处理位点文件，返回位点-区域、位点-RNA映射"""
    df = pd.read_csv(file_path)
    site_regions = {}
    site_to_rna = {}

    gene_type_cols = [col for col in df.columns if col.endswith('_Gene_Type')]
    gene_name_cols = [col for col in df.columns if col.endswith('_Gene_Name')]
    rna_id_col = 'RNA_ID' if 'RNA_ID' in df.columns else df.columns[2]

    for _, row in df.iterrows():
        site_id = f"{row['Chr']}:{row['Sites']}"
        rna_id = row[rna_id_col] if pd.notna(row[rna_id_col]) else site_id.split(':')[0] + '_RNA' + str(row.name)

        gene_types = [row[col] for col in gene_type_cols if pd.notna(row[col])]
        gene_names = [row[col] for col in gene_name_cols if pd.notna(row[col])]

        if gene_types and gene_names:
            region = parse_region(gene_types[0], gene_names[0])
            site_regions[site_id] = region
            site_to_rna[site_id] = rna_id

    return site_regions, site_to_rna


def main(file1, file2, output_file='modification_region_distribution.pdf'):
    # 应用绘图风格
    set_nature_methods_style()

    # 处理文件
    sites_37_regions, sites_37_rna = process_sites_file(file1)
    sites_45_regions, sites_45_rna = process_sites_file(file2)

    # 定义区域顺序
    region_order = ["ncRNA", "3' UTR", "Intergenic", "5' UTR", "CDS"]

    # 鉴定RNA的共有/特有
    rna_37 = set(sites_37_rna.values())
    rna_45 = set(sites_45_rna.values())
    shared_rna = rna_37 & rna_45
    unique_rna_37 = rna_37 - rna_45
    unique_rna_45 = rna_45 - rna_37

    # 初始化各层计数
    layer_counts = {
        'mod_unique37': [0] * len(region_order),
        'unique_37': [0] * len(region_order),
        'shared': [0] * len(region_order),
        'unique_45': [0] * len(region_order),
        'mod_unique45': [0] * len(region_order)
    }

    # 统计共有RNA上的修饰位点
    shared_site_37 = {s for s in sites_37_regions if sites_37_rna[s] in shared_rna}
    shared_site_45 = {s for s in sites_45_regions if sites_45_rna[s] in shared_rna}

    # 筛选共有修饰位点
    shared_mod_sites = set()
    for s37 in shared_site_37:
        chr_pos = s37.split(':')[-1]
        for s45 in shared_site_45:
            if s45.split(':')[-1] == chr_pos and sites_37_regions[s37] == sites_45_regions[s45]:
                shared_mod_sites.add(s37)
                shared_mod_sites.add(s45)

    # 逐区域统计各层
    for idx, region in enumerate(region_order):
        mod_unique45_sites = [s for s in sites_45_regions if
                              sites_45_regions[s] == region and sites_45_rna[s] in unique_rna_45]
        layer_counts['mod_unique45'][idx] = len(mod_unique45_sites)

        unique_45_sites = [s for s in shared_site_45 if sites_45_regions[s] == region and s not in shared_mod_sites]
        layer_counts['unique_45'][idx] = len(unique_45_sites)

        shared_sites = [s for s in shared_mod_sites if sites_37_regions.get(s, '') == region]
        layer_counts['shared'][idx] = len(shared_sites)

        unique_37_sites = [s for s in shared_site_37 if sites_37_regions[s] == region and s not in shared_mod_sites]
        layer_counts['unique_37'][idx] = len(unique_37_sites)

        mod_unique37_sites = [s for s in sites_37_regions if
                              sites_37_regions[s] == region and sites_37_rna[s] in unique_rna_37]
        layer_counts['mod_unique37'][idx] = len(mod_unique37_sites)

    # 归一化到Y轴100
    max_height = 100
    normalized = {}
    for layer in layer_counts:
        normalized[layer] = []
        for idx in range(len(region_order)):
            total = sum([layer_counts[l][idx] for l in layer_counts])
            normalized[layer].append((layer_counts[layer][idx] / total) * max_height if total != 0 else 0)

    # 绘制堆叠柱状图
    colors = {
        'Unique modification (OD1)': '#aec7e8',
        'Shared modification': '#c7c7c7',
        'Unique modification (OD0.3)': '#ffbb78',
    }

    fig, ax = plt.subplots(figsize=(6.5, 5.4))
    bar_width = 0.9
    index = np.arange(len(region_order))
    bottom = np.zeros(len(region_order))

    layers_draw_order = ['unique_45', 'shared', 'unique_37']
    labels = {
        'unique_45': 'Unique modification (OD0.3)',
        'shared': 'Shared modification',
        'unique_37': 'Unique modification (OD1)',
    }

    for layer in layers_draw_order:
        ax.bar(index, normalized[layer], bar_width, bottom=bottom,
               color=colors[labels[layer]], label=labels[layer],
               edgecolor='black', linewidth=0.8)
        bottom += np.array(normalized[layer])

    # 基础图形配置（保留原有逻辑）
    ax.set_xlabel('Region')
    ax.set_ylabel('Relative Proportion (%)')
    ax.set_xticks(index)
    ax.set_xticklabels(region_order)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xlim(-0.5, len(region_order) - 0.5)
    ax.set_ylim(0, max_height)

    plt.tight_layout()

    # 保存文件：先保存PNG，再保存PDF（默认输出为PDF）
    plt.savefig(output_file.replace('.pdf', '.png'), dpi=600, bbox_inches='tight', pad_inches=0.2)
    plt.savefig(output_file, dpi=600, bbox_inches='tight', pad_inches=0.2, format='pdf')
    print(f"PDF文件已保存至：{output_file}")
    print(f"PNG文件已保存至：{output_file.replace('.pdf', '.png')}")
    plt.show()


if __name__ == "__main__":
    file_37 = "OD_1_0.75_sites.txt"
    file_45 = "OD_03_0.75_sites.txt"
    main(file_37, file_45)

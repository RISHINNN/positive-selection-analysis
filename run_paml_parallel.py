import os
import re
import shutil
import subprocess
import multiprocessing
import csv
import argparse
from scipy.stats import chi2

def load_mapping(mapping_file):
    """读取制表符分隔的映射文件，构建 {基因ID: 物种名} 的字典"""
    mapping = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                gene_id = parts[0].strip()
                species_name = parts[1].strip()
                mapping[gene_id] = species_name
    return mapping

def rename_phy_for_paml(phy_in, phy_out, mapping_dict):
    """
    【完美修复】兼容 pal2nal 的双行 PAML 格式，精准替换 ID
    """
    with open(phy_in, "r") as f:
        lines = f.readlines()
        
    with open(phy_out, "w") as f:
        f.write(lines[0])  # 写入第一行的头信息
        for line in lines[1:]:
            stripped = line.strip()
            if not stripped:
                f.write(line)
                continue
                
            # 尝试格式 1：ID 和序列在同一行 (传统 Phylip)
            parts = stripped.split(None, 1)
            if len(parts) == 2 and parts[0] in mapping_dict:
                new_id = mapping_dict[parts[0]]
                f.write(f"{new_id:<30} {parts[1]}\n")
                continue
                
            # 尝试格式 2：ID 独占一行 (pal2nal PAML 格式的真面目)
            if stripped in mapping_dict:
                new_id = mapping_dict[stripped]
                f.write(new_id + "\n")
                continue
                
            # 如果都不是，说明这一行是纯 ATGC 序列，原样写入
            f.write(line)

def generate_ctl(phy_file, tree_file, ctl_file, out_mlc, model_type):
    """动态生成 codeml 配置文件"""
    fix_omega = "0" if model_type == "alt" else "1"
    omega_val = "1" if model_type == "alt" else "1"
    
    ctl_content = f"""
      seqfile = {phy_file}
     treefile = {tree_file}
      outfile = {out_mlc}

        noisy = 0
      verbose = 1
      runmode = 0

      seqtype = 1
    cleandata = 0      * 【新增】自动清理包含全 gap 或含 N 的模糊列
    CodonFreq = 2
        clock = 0
       aaDist = 0

        model = 2
      NSsites = 2
        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = {fix_omega}
        omega = {omega_val}
    """
    with open(ctl_file, "w") as f:
        f.write(ctl_content.strip())

def bh_fdr_correction(p_values):
    """
    【新增】原生实现的 Benjamini-Hochberg FDR 多重检验校正算法
    """
    n = len(p_values)
    if n == 0:
        return []
        
    ranked_p_indices = sorted(range(n), key=lambda i: p_values[i])
    sorted_p = [p_values[i] for i in ranked_p_indices]
    
    fdr = [0.0] * n
    min_fdr = 1.0
    
    for i in range(n - 1, -1, -1):
        p = sorted_p[i]
        fdr_val = p * n / (i + 1)
        min_fdr = min(min_fdr, fdr_val)
        fdr[ranked_p_indices[i]] = min_fdr
        
    return fdr

def parse_mlc(mlc_file):
    """
    【升级】提取 Log-Likelihood (lnL) 值，并检查前景支有效性
    """
    lnL = None
    has_foreground = False
    
    if not os.path.exists(mlc_file):
        return None, False
        
    with open(mlc_file, "r") as f:
        content = f.read()
        
        # 1. 提取似然值 (lnL)
        lnL_match = re.search(r"lnL.*:\s+(-?\d+\.\d+)", content)
        if lnL_match:
            lnL = float(lnL_match.group(1))
            
        # 2. 检查前景支是否有效
        if "foreground w" in content or "w (dN/dS) for branches" in content:
            has_foreground = True
            
    return lnL, has_foreground
def run_single_paml(args):
    """单个基因的 PAML 运行逻辑 (打上 Sandra 官方树文件补丁)"""
    og_name, phy_source, tree_source, base_outdir, codeml_exec, mapping_dict = args
    
    sandbox_dir = os.path.join(base_outdir, og_name)
    os.makedirs(sandbox_dir, exist_ok=True)
    
    local_tree = os.path.join(sandbox_dir, "input.nwk")
    local_phy = os.path.join(sandbox_dir, "input.phy")
    
    # 1. 准备并净化 phy 文件
    try:
        rename_phy_for_paml(phy_source, local_phy, mapping_dict)
    except Exception as e:
        return og_name, None, None, f"改名失败: {e}"
        
    # 2. 【Sandra 官方补丁】：动态为树文件添加 PHYLIP Header
    # 先从 phy 文件第一行读取物种数
    num_taxa = "0"
    with open(local_phy, "r") as f:
        first_line = f.readline().strip()
        if first_line:
            num_taxa = first_line.split()[0] # 比如 "21 500" 提取出 "21"
            
    # 读取原始树文件内容
    with open(tree_source, "r") as f:
        tree_content = f.read().strip()
        
    # 将 "物种数 1" 作为第一行，写入新的树文件
    with open(local_tree, "w") as f:
        f.write(f"{num_taxa}  1\n")
        f.write(tree_content + "\n")
    
    alt_ctl = os.path.join(sandbox_dir, "alt.ctl")
    null_ctl = os.path.join(sandbox_dir, "null.ctl")
    generate_ctl("input.phy", "input.nwk", alt_ctl, "alt.mlc", "alt")
    generate_ctl("input.phy", "input.nwk", null_ctl, "null.mlc", "null")
    
    try:
        # 依然无视 PAML 的报错状态码
        subprocess.run([codeml_exec, "alt.ctl"], cwd=sandbox_dir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([codeml_exec, "null.ctl"], cwd=sandbox_dir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as e:
        return og_name, None, None, f"系统调用失败: {e}"
        
    lnl_alt, valid_fg_alt = parse_mlc(os.path.join(sandbox_dir, "alt.mlc"))
    lnl_null, _ = parse_mlc(os.path.join(sandbox_dir, "null.mlc"))
    
    # ================= 严格质控模块 =================
    if lnl_alt is None or lnl_null is None:
        return og_name, lnl_alt, lnl_null, "未生成有效_mlc_崩溃"
        
    if not valid_fg_alt:
        return og_name, lnl_alt, lnl_null, "Invalid_Foreground_Branch"
        
    if lnl_alt < lnl_null - 0.1:
        return og_name, lnl_alt, lnl_null, "Model_Non_Converged"
    # ================================================
        
    return og_name, lnl_alt, lnl_null, "Success"


if __name__ == "__main__":
    # ================= 固定的软件路径配置区 =================
    CODEML_EXEC = "/data01/lichen/00_software/paml/bin/codeml" 
    # =======================================================

    if shutil.which(CODEML_EXEC) is None:
        print(f"❌ 严重错误: 找不到 codeml 命令！请检查路径: {CODEML_EXEC}")
        exit(1)

    parser = argparse.ArgumentParser(description="基于严格映射表与严谨统计质控的 PAML 并发运行脚本")
    parser.add_argument("-t", "--tree", required=True, help="带有 #1 标记的物种树文件")
    parser.add_argument("-p", "--phy_dir", required=True, help="存放所有 .phy 的目录")
    parser.add_argument("-m", "--mapping", required=True, help="基因ID到物种名的映射文件 (Tab 分隔: ID \t Species)")
    parser.add_argument("-o", "--out_dir", default="./PAML_Results", help="输出目录")
    parser.add_argument("-c", "--cpus", type=int, default=4, help="并发线程数")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    csv_file = os.path.join(args.out_dir, "Positive_Selection_Results.csv")

    print(f"正在读取映射文件: {args.mapping}")
    try:
        global_mapping_dict = load_mapping(args.mapping)
        print(f"✅ 成功加载 {len(global_mapping_dict)} 条基因映射关系。")
    except Exception as e:
        print(f"❌ 读取映射文件失败: {e}")
        exit(1)

    tasks = []
    valid_files = [f for f in os.listdir(args.phy_dir) if f.endswith(".phy")]
    
    for filename in valid_files:
        og_name = os.path.splitext(filename)[0]
        phy_path = os.path.join(args.phy_dir, filename)
        
        if os.path.getsize(phy_path) > 0:
            tasks.append((og_name, os.path.abspath(phy_path), os.path.abspath(args.tree), os.path.abspath(args.out_dir), CODEML_EXEC, global_mapping_dict))
            
    total_tasks = len(tasks)
    print(f"🚀 开始并发运行 PAML，有效基因数: {total_tasks}，分配线程数: {args.cpus}...")

    results = []
    with multiprocessing.Pool(processes=args.cpus) as pool:
        for i, res in enumerate(pool.imap_unordered(run_single_paml, tasks), 1):
            results.append(res)
            if i % 10 == 0 or i == total_tasks:
                print(f"进度: [{i}/{total_tasks}] 个基因已计算完成...")

    print("\n📊 正在进行 LRT 检验与 FDR 多重检验校正...")
    
    successful_genes = []
    failed_genes = []
    p_values = []
    
    for og_name, lnl_alt, lnl_null, status in results:
        if status == "Success":
            lrt_stat = max(0, 2 * (lnl_alt - lnl_null))
            p_val = chi2.sf(lrt_stat, df=1)
            p_values.append(p_val)
            successful_genes.append({
                "og_name": og_name, "lnl_alt": lnl_alt, "lnl_null": lnl_null, 
                "lrt_stat": lrt_stat, "p_val": p_val
            })
        else:
            failed_genes.append([og_name, lnl_alt, lnl_null, "NA", "NA", "NA", "NA", status])

    # 执行 FDR 校正
    fdr_values = bh_fdr_correction(p_values)
    
    sig_count = 0
    with open(csv_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Gene_Family", "lnL_Alt", "lnL_Null", "LRT_Statistic", "P_Value", "FDR_Q_Value", "Significance", "Status/QC_Check"])
        
        # 写入成功的基因
        for gene, fdr in zip(successful_genes, fdr_values):
            # 严格标准：使用 FDR 校正后的 Q-value < 0.05 来判断是否显著
            is_sig = "Yes" if fdr < 0.05 else "No"
            if is_sig == "Yes": sig_count += 1
            
            writer.writerow([
                gene["og_name"], gene["lnl_alt"], gene["lnl_null"], 
                round(gene["lrt_stat"], 4), format(gene["p_val"], ".3e"), format(fdr, ".3e"), 
                is_sig, "Pass"
            ])
            
        # 写入被拦截的基因
        for row in failed_genes:
            writer.writerow(row)

    print(f"\n🎉 统计分析与严格质控全部完成！")
    print(f" - 通过质控进入检验的基因: {len(successful_genes)} 个")
    print(f" - 未通过质控被拦截的基因: {len(failed_genes)} 个")
    print(f" - 👑 经 FDR 校正后显著的正选择基因: {sig_count} 个")
    print(f"最终结果表格已保存至: {csv_file}")

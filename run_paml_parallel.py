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

def parse_lnl(mlc_file):
    """提取 Log-Likelihood (lnL) 值"""
    if not os.path.exists(mlc_file):
        return None
    with open(mlc_file, "r") as f:
        for line in f:
            if line.startswith("lnL"):
                match = re.search(r"lnL.*:\s+(-?\d+\.\d+)", line)
                if match:
                    return float(match.group(1))
    return None

def run_single_paml(args):
    """单个基因的 PAML 运行逻辑 (沙盒模式)"""
    # 显式接收 mapping_dict，防止多进程变量丢失
    og_name, phy_source, tree_source, base_outdir, codeml_exec, mapping_dict = args
    
    sandbox_dir = os.path.join(base_outdir, og_name)
    os.makedirs(sandbox_dir, exist_ok=True)
    
    local_tree = os.path.join(sandbox_dir, "input.nwk")
    local_phy = os.path.join(sandbox_dir, "input.phy")
    
    shutil.copy2(tree_source, local_tree)
    
    try:
        rename_phy_for_paml(phy_source, local_phy, mapping_dict)
    except Exception as e:
        return og_name, None, None, f"改名失败: {e}"
    
    alt_ctl = os.path.join(sandbox_dir, "alt.ctl")
    null_ctl = os.path.join(sandbox_dir, "null.ctl")
    generate_ctl("input.phy", "input.nwk", alt_ctl, "alt.mlc", "alt")
    generate_ctl("input.phy", "input.nwk", null_ctl, "null.mlc", "null")
    
    try:
        subprocess.run([codeml_exec, "alt.ctl"], cwd=sandbox_dir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        subprocess.run([codeml_exec, "null.ctl"], cwd=sandbox_dir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except Exception as e:
        return og_name, None, None, f"codeml 运行失败: {e}"
        
    lnl_alt = parse_lnl(os.path.join(sandbox_dir, "alt.mlc"))
    lnl_null = parse_lnl(os.path.join(sandbox_dir, "null.mlc"))
    
    if lnl_alt is None or lnl_null is None:
        return og_name, None, None, "未能提取 lnL (序列可能差异过大导致 PAML 崩溃)"
        
    return og_name, lnl_alt, lnl_null, "Success"

if __name__ == "__main__":
    # ================= 固定的软件路径配置区 =================
    # 请确保这里是您测试通过的 codeml 绝对路径！
    CODEML_EXEC = "/data01/lichen/00_software/paml/bin/codeml" 
    # =======================================================

    if shutil.which(CODEML_EXEC) is None:
        print(f"❌ 严重错误: 找不到 codeml 命令！请检查路径: {CODEML_EXEC}")
        exit(1)

    parser = argparse.ArgumentParser(description="基于严格映射表的 PAML 并发运行脚本")
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
            # 显式把字典传进去
            tasks.append((og_name, os.path.abspath(phy_path), os.path.abspath(args.tree), os.path.abspath(args.out_dir), CODEML_EXEC, global_mapping_dict))
            
    total_tasks = len(tasks)
    print(f"🚀 开始并发运行 PAML，有效基因数: {total_tasks}，分配线程数: {args.cpus}...")

    results = []
    with multiprocessing.Pool(processes=args.cpus) as pool:
        for i, res in enumerate(pool.imap_unordered(run_single_paml, tasks), 1):
            results.append(res)
            if i % 10 == 0 or i == total_tasks:
                print(f"进度: [{i}/{total_tasks}] 个基因已计算完成...")

    print("\n📊 运行结束，正在计算似然比检验 (LRT)...")
    
    sig_count = 0
    with open(csv_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Gene_Family", "lnL_Alt", "lnL_Null", "LRT_Statistic", "P_Value", "Significance", "Status"])
        
        for og_name, lnl_alt, lnl_null, status in results:
            if lnl_alt is not None and lnl_null is not None:
                lrt_stat = max(0, 2 * (lnl_alt - lnl_null))
                p_value = chi2.sf(lrt_stat, df=1)
                
                is_sig = "Yes" if p_value < 0.05 else "No"
                if is_sig == "Yes": sig_count += 1
                
                writer.writerow([og_name, lnl_alt, lnl_null, round(lrt_stat, 4), format(p_value, ".2e"), is_sig, status])
            else:
                writer.writerow([og_name, "NA", "NA", "NA", "NA", "NA", status])

    print(f"\n🎉 在 {total_tasks} 个基因中，初步发现 {sig_count} 个基因受到显著正选择 (P < 0.05)。")
    print(f"最终结果表格已保存至: {csv_file}")

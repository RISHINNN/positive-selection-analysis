import os
import subprocess
import argparse
import shutil
import multiprocessing
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

global_pep_dict = {}
global_cds_dict = {}

def fasta_to_paml_phylip(fasta_in, phy_out):
    """将 FASTA 格式安全转换为 PAML 兼容的双行 PHYLIP 格式，防止 ID 被截断"""
    records = list(SeqIO.parse(fasta_in, "fasta"))
    if not records:
        return False
        
    seq_len = len(records[0].seq)
    # 如果修剪后序列全被删光了 (长度为0)，返回 False
    if seq_len == 0:
        return False
        
    with open(phy_out, "w") as f:
        f.write(f"  {len(records)}  {seq_len}\n")
        for rec in records:
            # PAML 格式：第一行 ID，第二行序列
            f.write(f"{rec.id}\n{str(rec.seq)}\n")
    return True

def process_single_og(og_name, unaligned_fa_path, out_dir, mafft_exec, pal2nal_script, trimal_exec, gap_threshold, error_log):
    """包含 trimAl 智能修剪的单基因处理核心"""
    temp_dir = os.path.join(out_dir, "temp_files")
    raw_pep = os.path.join(temp_dir, f"{og_name}_raw.pep")
    raw_cds = os.path.join(temp_dir, f"{og_name}_raw.cds")
    aln_pep = os.path.join(temp_dir, f"{og_name}_aligned.pep")
    p2n_fa = os.path.join(temp_dir, f"{og_name}_p2n.fasta")
    trimmed_fa = os.path.join(temp_dir, f"{og_name}_trimmed.fasta")
    out_phy = os.path.join(out_dir, f"{og_name}.phy")
    
    try:
        og_ids = [record.id for record in SeqIO.parse(unaligned_fa_path, "fasta")]
        pep_recs, cds_recs = [], []
        
        for pid in og_ids:
            if pid not in global_pep_dict or pid not in global_cds_dict:
                return "error", f"在原始库中找不到 ID: {pid}"
            pep_recs.append(SeqRecord(global_pep_dict[pid].seq, id=pid, description=""))
            cds_recs.append(SeqRecord(global_cds_dict[pid].seq, id=pid, description=""))
            
        SeqIO.write(pep_recs, raw_pep, "fasta")
        SeqIO.write(cds_recs, raw_cds, "fasta")
        
        # 1. MAFFT 比对
        subprocess.run([mafft_exec, "--quiet", "--auto", raw_pep], stdout=open(aln_pep, "w"), stderr=subprocess.PIPE, check=True)
                
        # 2. Pal2Nal 映射 (取消 -nogap，输出 fasta 格式以供修剪)
        pal2nal_cmd = ["perl", pal2nal_script, aln_pep, raw_cds, "-output", "fasta"]
        run_p2n = subprocess.run(pal2nal_cmd, stdout=open(p2n_fa, "w"), stderr=subprocess.PIPE, text=True)
        
        if os.path.getsize(p2n_fa) == 0:
            with open(error_log, "a") as err_f: err_f.write(f"\n[{og_name}] PAL2NAL 失败: {run_p2n.stderr}\n")
            return "failed", None
            
        # 3. trimAl 智能修剪 (-gt 参数决定保留列的阈值)
        trimal_cmd = [trimal_exec, "-in", p2n_fa, "-out", trimmed_fa, "-gt", str(gap_threshold)]
        subprocess.run(trimal_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=True)
        
        # 4. 转换为 PAML 格式
        success = fasta_to_paml_phylip(trimmed_fa, out_phy)
        
        if not success:
            return "trimmed_to_zero", None
            
        return "success", None
    
    except Exception as e:
        return "error", str(e)

def unpack_args(args):
    return process_single_og(*args)

if __name__ == "__main__":
    # ================= 固定的软件路径配置区 =================
    PAL2NAL_SCRIPT = "/data01/lichen/03_tools/pal2nal.pl" 
    MAFFT_EXEC = "/public/home/wangwen_lab/duanlei/miniconda3/envs/ortholog/bin/mafft" 
    TRIMAL_EXEC = "/data01/lichen/03_tools/trimal"
    # =======================================================

    if not all(shutil.which(cmd) for cmd in [MAFFT_EXEC, TRIMAL_EXEC]):
        print(f"❌ 找不到 mafft 或 trimal，请检查环境变量！")
        exit(1)

    parser = argparse.ArgumentParser(description="带 trimAl 修剪的正选择前端数据准备管线")
    parser.add_argument("-p", "--pep", required=True)
    parser.add_argument("-c", "--cds", required=True)
    parser.add_argument("-i", "--input_dir", required=True)
    parser.add_argument("-o", "--out_dir", default="./Codon_Alignments_Trimmed")
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("-gt", "--gap_threshold", type=float, default=0.7, 
                        help="trimAl 的 -gt 参数 (默认 0.7: 保留至少 70% 物种有真实碱基的列)")
    args = parser.parse_args()

    os.makedirs(os.path.join(args.out_dir, "temp_files"), exist_ok=True)
    error_log = os.path.join(args.out_dir, "pipeline_error.log")
    open(error_log, "w").close()

    print("正在加载序列数据库...")
    global_pep_dict = SeqIO.to_dict(SeqIO.parse(args.pep, "fasta"))
    global_cds_dict = SeqIO.to_dict(SeqIO.parse(args.cds, "fasta"))

    valid_files = [f for f in os.listdir(args.input_dir) if f.endswith((".fa", ".fasta"))]
    total = len(valid_files)
    
    tasks = [(os.path.splitext(f)[0], os.path.join(args.input_dir, f), args.out_dir, MAFFT_EXEC, PAL2NAL_SCRIPT, TRIMAL_EXEC, args.gap_threshold, error_log) for f in valid_files]

    print(f"\n🚀 开始 {args.threads} 线程并发处理 (修剪阈值: -gt {args.gap_threshold})...")
    stats = {"success": 0, "failed": 0, "trimmed_to_zero": 0, "error": 0}
    
    with multiprocessing.Pool(processes=args.threads) as pool:
        for i, (status, err_msg) in enumerate(pool.imap_unordered(unpack_args, tasks), 1):
            if status == "success": stats["success"] += 1
            elif status == "failed": stats["failed"] += 1
            elif status == "trimmed_to_zero": stats["trimmed_to_zero"] += 1
            else:
                stats["error"] += 1
                with open(error_log, "a") as err_f: err_f.write(f"\n[Error] {err_msg}\n")
                
            if i % 100 == 0 or i == total:
                print(f"进度: [{i}/{total}] 已完成...")

    print("\n🎉 修剪管线执行完毕！")
    print(f" - 👑 成功生成高质量比对: {stats['success']} 个")
    print(f" - ⚠️ 序列有硬伤 (pal2nal 报错): {stats['failed']} 个")
    print(f" - 🗑️ 差异过大被全量修剪 (全垃圾序列): {stats['trimmed_to_zero']} 个")
    print(f" - ❌ 提取异常: {stats['error']} 个")

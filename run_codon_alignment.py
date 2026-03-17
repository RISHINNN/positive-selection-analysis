import os
import subprocess
import argparse
import shutil
import multiprocessing
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# 使用全局变量让子进程共享庞大的字典内存 (利用 Linux 的 Copy-On-Write 机制，极大地节省内存)
global_pep_dict = {}
global_cds_dict = {}

def process_single_og(og_name, unaligned_fa_path, out_dir, mafft_exec, pal2nal_script, error_log):
    """单基因处理核心函数"""
    temp_dir = os.path.join(out_dir, "temp_files")
    raw_pep_path = os.path.join(temp_dir, f"{og_name}_raw.pep")
    raw_cds_path = os.path.join(temp_dir, f"{og_name}_raw.cds")
    aln_pep_path = os.path.join(temp_dir, f"{og_name}_aligned.pep")
    out_phy_path = os.path.join(out_dir, f"{og_name}.phy")
    
    try:
        og_ids = [record.id for record in SeqIO.parse(unaligned_fa_path, "fasta")]
        pep_records_to_write = []
        cds_records_to_write = []
        
        for pid in og_ids:
            if pid not in global_pep_dict:
                return "error", f"在原始 all.pep 中找不到 ID: {pid}"
            if pid not in global_cds_dict:
                return "error", f"在原始 all.cds 中找不到 ID: {pid}"
                
            pep_records_to_write.append(SeqRecord(global_pep_dict[pid].seq, id=pid, description=""))
            cds_records_to_write.append(SeqRecord(global_cds_dict[pid].seq, id=pid, description=""))
            
        SeqIO.write(pep_records_to_write, raw_pep_path, "fasta")
        SeqIO.write(cds_records_to_write, raw_cds_path, "fasta")
        
        # MAFFT 保持单线程高效运行，把多核资源留给外层的 Python 并发
        mafft_cmd = [mafft_exec, "--quiet", "--auto", raw_pep_path]
        with open(aln_pep_path, "w") as f_out:
            result = subprocess.run(mafft_cmd, stdout=f_out, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                return "error", f"MAFFT 失败: {result.stderr}"
                
        # Pal2Nal 运行
        pal2nal_cmd = ["perl", pal2nal_script, aln_pep_path, raw_cds_path, "-output", "paml", "-nogap"]
        run_p2n = subprocess.run(pal2nal_cmd, stdout=open(out_phy_path, "w"), stderr=subprocess.PIPE, text=True)
            
        if os.path.getsize(out_phy_path) == 0:
            pal2nal_loose = ["perl", pal2nal_script, aln_pep_path, raw_cds_path, "-output", "paml"]
            subprocess.run(pal2nal_loose, stdout=open(out_phy_path, "w"), stderr=subprocess.PIPE)
            if os.path.getsize(out_phy_path) == 0:
                with open(error_log, "a") as err_f:
                    err_f.write(f"\n[{og_name}] PAL2NAL 失败。\n{run_p2n.stderr}\n")
                return "failed", None
            return "loose_success", None
            
        return "strict_success", None
    
    except Exception as e:
        return "error", str(e)


# 用于解包多进程参数的包装器
def unpack_args(args):
    return process_single_og(*args)


if __name__ == "__main__":
    PAL2NAL_SCRIPT = "/data01/lichen/03_tools/pal2nal.pl" 
    MAFFT_EXEC = "/public/home/wangwen_lab/duanlei/miniconda3/envs/ortholog/bin/mafft" 

    if shutil.which(MAFFT_EXEC) is None:
        print(f"❌ 找不到 mafft！")
        exit(1)

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pep", required=True)
    parser.add_argument("-c", "--cds", required=True)
    parser.add_argument("-i", "--input_dir", required=True)
    parser.add_argument("-o", "--out_dir", default="./Codon_Alignments_Rebuild")
    parser.add_argument("-t", "--threads", type=int, default=4, help="并发运行的线程/进程数 (默认: 4)")
    args = parser.parse_args()

    os.makedirs(os.path.join(args.out_dir, "temp_files"), exist_ok=True)
    error_log = os.path.join(args.out_dir, "pipeline_error.log")
    open(error_log, "w").close()

    print("正在加载全基因组序列到共享内存...")
    global_pep_dict = SeqIO.to_dict(SeqIO.parse(args.pep, "fasta"))
    global_cds_dict = SeqIO.to_dict(SeqIO.parse(args.cds, "fasta"))

    valid_files = [f for f in os.listdir(args.input_dir) if f.endswith((".fa", ".fasta"))]
    total = len(valid_files)
    
    print(f"\n🚀 开始 {args.threads} 线程并发处理，共 {total} 个基因家族...")

    # 准备任务列表
    tasks = []
    for filename in valid_files:
        og_name = os.path.splitext(filename)[0]
        unaligned_fa_path = os.path.join(args.input_dir, filename)
        tasks.append((og_name, unaligned_fa_path, args.out_dir, MAFFT_EXEC, PAL2NAL_SCRIPT, error_log))

    stats = {"strict": 0, "loose": 0, "failed": 0, "error": 0}
    
    # 核心：启动多进程池
    with multiprocessing.Pool(processes=args.threads) as pool:
        for i, (status, err_msg) in enumerate(pool.imap_unordered(unpack_args, tasks), 1):
            if status == "strict_success": stats["strict"] += 1
            elif status == "loose_success": stats["loose"] += 1
            elif status == "failed": stats["failed"] += 1
            else:
                stats["error"] += 1
                with open(error_log, "a") as err_f: err_f.write(f"\n[Error] {err_msg}\n")
                
            if i % 100 == 0 or i == total:
                print(f"进度: [{i}/{total}] 已完成...")

    print(f"\n🎉 极速并发执行完毕！(优质 .phy 保存在: {args.out_dir})")

import os
import argparse
from ete3 import Tree

def mark_foreground_branch(tree_path, target_species, output_tree_path):
    """
    读取物种树，寻找目标物种的最近共同祖先，并加上 #1 标记。
    使用物理切除法，彻底解决 ete3 顽固的 NoName 问题。
    """
    print(f"正在读取物种树: {tree_path}")
    t = Tree(tree_path)

    # 检查指定的物种是否都在树中
    tree_leaves = [leaf.name for leaf in t.get_leaves()]
    for sp in target_species:
        if sp not in tree_leaves:
            raise ValueError(f"错误: 物种 '{sp}' 不在物种树中！请检查拼写。")

    # 获取目标物种的最近共同祖先 (MRCA)
    if len(target_species) == 1:
        ancestor = t.search_nodes(name=target_species[0])[0]
    else:
        ancestor = t.get_common_ancestor(*target_species)

    # 遍历所有节点，进行精准标记和初步清理
    for node in t.traverse():
        if node.is_leaf():
            if node == ancestor and "#1" not in node.name:
                node.name = node.name + "#1"
        else:
            if node == ancestor:
                node.name = "#1"
            else:
                node.name = ""  # 清空其他内部节点名字

    # 【终极物理切除】
    # 生成字符串
    raw_newick = t.write(format=8)
    
    # 无情地把所有 ")NoName" 替换成单纯的 ")"
    clean_newick = raw_newick.replace(")NoName", ")")

    # 将真正纯净的树写入文件
    with open(output_tree_path, "w") as f:
        f.write(clean_newick)
        
    print(f"✅ 前景支已成功标记！纯净的树已保存至: {output_tree_path}")



def generate_codeml_ctl(alignment_file, tree_file, output_ctl, out_file="mlc"):
    """
    生成 PAML 分支-位点模型 (Branch-Site Model A) 的控制文件。
    """
    ctl_content = f"""
      seqfile = {alignment_file}   * 密码子比对序列文件 (PHYLIP格式)
     treefile = {tree_file}        * 带有 #1 标记的树文件
      outfile = {out_file}         * 结果输出文件

        noisy = 9      * 控制台输出级别 (0,1,2,3,9)
      verbose = 1      * 详细输出 (0: 简洁, 1: 详细)
      runmode = 0      * 用户自定义树

      seqtype = 1      * 1: 密码子序列 (必须)
    CodonFreq = 2      * 0: 1/61, 1: F1X4, 2: F3X4, 3: F61
        clock = 0      * 0: 无分子钟限制
       aaDist = 0      * 氨基酸距离

        model = 2      * 2: 分支-位点模型 (Branch-site models)
      NSsites = 2      * 2: selection (与 model=2 结合使用即为 Branch-site Model A)
        icode = 0      * 0: 标准遗传密码
    fix_kappa = 0      * 0: 估算 kappa, 1: 固定 kappa
        kappa = 2      * 初始或固定的 kappa 值
    fix_omega = 0      * 0: 估算 omega (dN/dS), 1: 固定 omega (用于 null model)
        omega = 1      * 初始 omega 值
    """
    with open(output_ctl, "w") as f:
        f.write(ctl_content.strip())
    print(f"✅ PAML 配置文件已生成: {output_ctl}")


if __name__ == "__main__":
    # ================= 命令行参数解析 =================
    parser = argparse.ArgumentParser(description="PAML 分支位点模型前期准备脚本：标记系统发育树并生成 ctl 文件")
    
    parser.add_argument("-t", "--tree", required=True, 
                        help="输入的物种树文件路径 (例如: SpeciesTree_rooted.txt)")
                        
    parser.add_argument("-c", "--clade", required=True, nargs="+", 
                        help="指定前景支的物种名称。可以是一个物种，也可以是多个物种（用空格分隔，脚本会自动寻找最近共同祖先）")
                        
    parser.add_argument("-o", "--out_tree", default="labeled_tree.nwk", 
                        help="输出的带标记的树文件路径 (默认: labeled_tree.nwk)")
                        
    parser.add_argument("-a", "--aln", default="dummy.phy", 
                        help="用于生成配置文件的比对文件占位符 (默认: dummy.phy，如果你使用批处理脚本，此项可忽略)")

    args = parser.parse_args()
    # =================================================

    print(f"输入的物种树: {args.tree}")
    print(f"设定的前景支物种: {args.clade}")

    # 步骤 1: 标记系统发育树
    try:
        mark_foreground_branch(args.tree, args.clade, args.out_tree)
    except Exception as e:
        print(e)
        exit(1)

    # 步骤 2: 生成 codeml 配置文件
    # 备择假设 (Alternative Model)
    generate_codeml_ctl(args.aln, args.out_tree, "alt_model.ctl", "alt_model.mlc")
    
    # 零假设 (Null Model)
    generate_codeml_ctl(args.aln, args.out_tree, "null_model.ctl", "null_model.mlc")
    # 修改 Null 文件的最后两行
    with open("null_model.ctl", "r") as f:
        lines = f.readlines()
    with open("null_model.ctl", "w") as f:
        for line in lines:
            if line.startswith("    fix_omega"):
                f.write("    fix_omega = 1      * 1: 固定 omega\n")
            elif line.startswith("        omega"):
                f.write("        omega = 1      * 将 omega 固定为 1\n")
            else:
                f.write(line)
    print("✅ 零假设 (Null Model) 的 PAML 配置文件已生成: null_model.ctl")
    print("\n🎉 准备完成！模板文件已生成。")

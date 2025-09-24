#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本用于只从moltype = "lipids"部分移除缺失的脂质条目
"""

import re

def remove_missing_lipids_from_lipids_section(py_file_path, missing_lipids):
    """只从moltype = "lipids"部分移除缺失的脂质条目"""
    
    with open(py_file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 找到moltype = "lipids"部分的开始和结束
    lipids_start_pattern = r'moltype = "lipids"'
    lipids_end_pattern = r'moltype = "PI45P2"'
    
    # 找到lipids部分的开始位置
    lipids_start_match = re.search(lipids_start_pattern, content)
    if not lipids_start_match:
        print("未找到 moltype = \"lipids\" 部分")
        return content, 0
    
    # 找到lipids部分的结束位置
    lipids_end_match = re.search(lipids_end_pattern, content)
    if not lipids_end_match:
        print("未找到 moltype = \"PI45P2\" 部分")
        return content, 0
    
    start_pos = lipids_start_match.start()
    end_pos = lipids_end_match.start()
    
    # 提取lipids部分
    lipids_section = content[start_pos:end_pos]
    
    # 统计移除的条目数量
    removed_count = 0
    
    # 对每个缺失的脂质，移除所有相关的条目
    for lipid in missing_lipids:
        # 创建匹配模式，匹配整个脂质条目行
        pattern = rf'^\s*"{lipid}":\s*\(moltype,\s*"[^"]*"\),?\s*$'
        
        # 使用多行模式匹配
        matches = re.findall(pattern, lipids_section, re.MULTILINE)
        if matches:
            # 移除匹配的行
            lipids_section = re.sub(pattern, '', lipids_section, flags=re.MULTILINE)
            removed_count += len(matches)
            print(f"从lipids部分移除了 {len(matches)} 个 {lipid} 条目")
    
    # 清理多余的空行（连续的空行）
    lipids_section = re.sub(r'\n\s*\n\s*\n', '\n\n', lipids_section)
    
    # 重新组合文件内容
    modified_content = content[:start_pos] + lipids_section + content[end_pos:]
    
    return modified_content, removed_count

def main():
    # 文件路径
    py_file = "/Users/renxinyu/LNB/LNB-MDT_v1.0/generation/lnb_gener_martini3.py"
    missing_file = "/Users/renxinyu/LNB/LNB-MDT_v1.0/missing_lipids.txt"
    
    # 读取缺失的脂质列表
    with open(missing_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # 提取脂质名称（跳过标题行）
    missing_lipids = []
    for line in lines:
        line = line.strip()
        if line and not line.startswith('在') and not line.startswith('='):
            missing_lipids.append(line)
    
    print(f"需要从lipids部分移除的缺失脂质: {len(missing_lipids)} 个")
    print("缺失脂质列表:", missing_lipids)
    
    # 只移除lipids部分的缺失脂质条目
    print("\n开始从moltype = \"lipids\"部分移除缺失的脂质条目...")
    modified_content, removed_count = remove_missing_lipids_from_lipids_section(py_file, missing_lipids)
    
    # 备份原文件
    backup_file = py_file + ".backup2"
    with open(backup_file, 'w', encoding='utf-8') as f:
        with open(py_file, 'r', encoding='utf-8') as original:
            f.write(original.read())
    print(f"原文件已备份到: {backup_file}")
    
    # 写入修改后的内容
    with open(py_file, 'w', encoding='utf-8') as f:
        f.write(modified_content)
    
    print(f"\n修改完成!")
    print(f"从lipids部分总共移除了 {removed_count} 个脂质条目")
    print(f"修改后的文件已保存到: {py_file}")

if __name__ == "__main__":
    main()

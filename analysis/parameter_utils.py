"""
参数解析工具模块

提供简化的参数输入格式支持，让用户更容易输入residues和gas-group参数。
"""

import ast
import os
import json

def parse_residues_simple(residues_input):
    """
    解析简化的residues输入格式
    
    支持的格式：
    1. 字典字符串: "{'DPPC': ['PO4'], 'CHOL': ['ROH']}"
    2. 简单格式: "DPPC:PO4,CHOL:ROH" 
    3. 多原子格式: "DPPC:PO4+GLY,CHOL:ROH"
    
    Args:
        residues_input: 输入的residues字符串
        
    Returns:
        dict: 解析后的residues字典
    """
    if not residues_input:
        return {}
    
    # 尝试解析为字典字符串
    try:
        result = ast.literal_eval(residues_input)
        if isinstance(result, dict):
            return result
    except (ValueError, SyntaxError):
        pass
    
    # 解析简单格式
    return _parse_simple_format(residues_input)

def parse_gas_group_simple(gas_input):
    """
    解析简化的gas-group输入格式
    
    支持的格式：
    1. 字典字符串: "{'N2': ['N2']}"
    2. 简单格式: "N2:N2" 或 "N2"
    3. 多气体格式: "N2:N2,O2:O2"
    
    Args:
        gas_input: 输入的gas-group字符串
        
    Returns:
        dict: 解析后的gas-group字典
    """
    if not gas_input:
        return {}
    
    # 尝试解析为字典字符串
    try:
        result = ast.literal_eval(gas_input)
        if isinstance(result, dict):
            return result
    except (ValueError, SyntaxError):
        pass
    
    # 解析简单格式
    return _parse_simple_format(gas_input)

def _parse_simple_format(input_str):
    """
    解析简单格式: "RESIDUE:ATOM1+ATOM2,RESIDUE2:ATOM3"
    """
    result = {}
    
    # 按逗号分割
    parts = input_str.split(',')
    
    for part in parts:
        part = part.strip()
        if not part:
            continue
            
        # 检查是否包含冒号
        if ':' in part:
            residue, atoms_str = part.split(':', 1)
            residue = residue.strip()
            atoms_str = atoms_str.strip()
            
            # 处理多个原子（用+分隔）
            if '+' in atoms_str:
                atoms = [atom.strip() for atom in atoms_str.split('+')]
            else:
                atoms = [atoms_str]
            
            result[residue] = atoms
        else:
            # 如果没有冒号，假设residue和atom名称相同
            residue = part.strip()
            result[residue] = [residue]
    
    return result


if __name__ == "__main__":
    # 测试功能
    print("参数解析工具测试")
    print("=" * 40)
    
    # 测试简单格式
    test_cases = [
        "DPPC:PO4,CHOL:ROH",
        "N2:N2,O2:O2",
        "DPPC:PO4+GLY,CHOL:ROH",
        "N2",  # 只有气体名
    ]
    
    for test_input in test_cases:
        print(f"输入: {test_input}")
        result = _parse_simple_format(test_input)
        print(f"输出: {result}")
        print()
    

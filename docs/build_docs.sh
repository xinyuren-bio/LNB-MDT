#!/bin/bash

# LNB-MDT 文档构建脚本

echo "🚀 开始构建LNB-MDT文档..."

# 检查是否在正确的目录
if [ ! -f "source/conf.py" ]; then
    echo "❌ 错误: 请在docs目录中运行此脚本"
    exit 1
fi

# 检查Python环境
if ! command -v python &> /dev/null; then
    echo "❌ 错误: 未找到Python"
    exit 1
fi

# 检查pip
if ! command -v pip &> /dev/null; then
    echo "❌ 错误: 未找到pip"
    exit 1
fi

# 安装依赖
echo "📦 安装文档构建依赖..."
pip install -r requirements.txt

if [ $? -ne 0 ]; then
    echo "❌ 错误: 依赖安装失败"
    exit 1
fi

# 清理之前的构建
echo "🧹 清理之前的构建..."
rm -rf build/

# 构建HTML文档
echo "📚 构建HTML文档..."
make html

if [ $? -eq 0 ]; then
    echo "✅ 文档构建成功!"
    echo "📖 文档位置: build/html/index.html"
    echo "🌐 在浏览器中打开: file://$(pwd)/build/html/index.html"
    
    # 尝试在浏览器中打开
    if command -v open &> /dev/null; then
        echo "🔍 正在浏览器中打开文档..."
        open build/html/index.html
    elif command -v xdg-open &> /dev/null; then
        echo "🔍 正在浏览器中打开文档..."
        xdg-open build/html/index.html
    fi
else
    echo "❌ 错误: 文档构建失败"
    exit 1
fi

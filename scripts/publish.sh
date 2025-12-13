#!/bin/bash

# LNB-MDT PyPI 发布脚本
# 注意：请从项目根目录运行此脚本

set -e  # 遇到错误立即退出

echo "🚀 LNB-MDT PyPI 发布脚本"
echo "=========================="

# 获取脚本所在目录的父目录（项目根目录）
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# 切换到项目根目录
cd "$PROJECT_ROOT"

# 检查是否在正确的目录
if [ ! -f "setup.py" ]; then
    echo "❌ 错误: 未找到 setup.py，请确保在项目根目录运行此脚本"
    exit 1
fi

# 检查必要的工具
echo "📦 检查构建工具..."
if ! command -v python &> /dev/null; then
    echo "❌ 错误: 未找到 Python"
    exit 1
fi

# 安装/升级构建工具
echo "📦 安装/升级构建工具..."
pip install --upgrade build twine -q

# 清理之前的构建
echo "🧹 清理之前的构建文件..."
rm -rf build/ dist/ *.egg-info/ LNB_MDT.egg-info/

# 构建分发包
echo "🔨 构建分发包..."
python -m build

# 检查分发包
echo "✅ 检查分发包..."
twine check dist/*

echo ""
echo "✅ 构建完成！"
echo ""
echo "📦 生成的文件:"
ls -lh dist/
echo ""
echo "📝 下一步操作:"
echo ""
echo "1. 测试上传到 TestPyPI (推荐):"
echo "   twine upload --repository testpypi dist/*"
echo ""
echo "2. 上传到正式 PyPI:"
echo "   twine upload dist/*"
echo ""
echo "3. 测试安装:"
echo "   pip install lnb-mdt"
echo ""


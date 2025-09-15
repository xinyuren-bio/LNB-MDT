#!/bin/bash

# LNB-MDT 文档预览脚本
echo "🚀 构建LNB-MDT文档..."

# 进入docs目录
cd "$(dirname "$0")"

# 快速构建
echo "📚 构建HTML文档..."
sphinx-build -b html source build/html

if [ $? -eq 0 ]; then
    echo "✅ 构建完成！"
    echo ""
    echo "🌐 本地预览地址:"
    echo "   file://$(pwd)/build/html/index.html"
    echo ""
    echo "💡 提示: 修改文件后重新运行此脚本即可更新预览"
    
    # 自动打开浏览器
    if command -v open &> /dev/null; then
        open "file://$(pwd)/build/html/index.html"
    elif command -v xdg-open &> /dev/null; then
        xdg-open "file://$(pwd)/build/html/index.html"
    fi
else
    echo "❌ 构建失败"
    exit 1
fi

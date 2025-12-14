@echo off
REM LNB-MDT PyPI 发布脚本 (Windows)
REM 注意：请从项目根目录运行此脚本

echo 🚀 LNB-MDT PyPI 发布脚本
echo ==========================

REM 获取脚本所在目录的父目录（项目根目录）
set "SCRIPT_DIR=%~dp0"
set "PROJECT_ROOT=%SCRIPT_DIR%.."

REM 切换到项目根目录
cd /d "%PROJECT_ROOT%"

REM 检查是否在正确的目录
if not exist "setup.py" (
    echo ❌ 错误: 未找到 setup.py，请确保在项目根目录运行此脚本
    exit /b 1
)

REM 安装/升级构建工具
echo 📦 安装/升级构建工具...
python -m pip install --upgrade build twine -q

REM 清理之前的构建
echo 🧹 清理之前的构建文件...
if exist "build" rmdir /s /q build
if exist "dist" rmdir /s /q dist
if exist "*.egg-info" rmdir /s /q *.egg-info
if exist "LNB_MDT.egg-info" rmdir /s /q LNB_MDT.egg-info

REM 构建分发包
echo 🔨 构建分发包...
python -m build

REM 检查分发包
echo ✅ 检查分发包...
twine check dist/*

echo.
echo ✅ 构建完成！
echo.
echo 📦 生成的文件:
dir dist
echo.
echo 📝 下一步操作:
echo.
echo 1. 测试上传到 TestPyPI (推荐):
echo    twine upload --repository testpypi dist/*
echo.
echo 2. 上传到正式 PyPI:
echo    twine upload dist/*
echo.
echo 3. 测试安装:
echo    pip install lnb-mdt
echo.

pause


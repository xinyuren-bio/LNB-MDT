# 开发者工具脚本

本目录包含用于发布 LNB-MDT 到 PyPI 的开发者工具脚本。

**注意**：这些脚本仅供项目维护者使用，普通用户不需要这些文件。

## 文件说明

- `publish.sh` - Linux/macOS 发布脚本
- `publish.bat` - Windows 发布脚本
- `PYPI_UPLOAD_GUIDE.md` - PyPI 上传详细指南

## 使用方法

**重要**：这些脚本需要从项目根目录运行，而不是从 `scripts/` 目录运行。

### Linux/macOS

```bash
# 从项目根目录运行
./scripts/publish.sh
```

### Windows

```cmd
REM 从项目根目录运行
scripts\publish.bat
```

详细说明请参考 `PYPI_UPLOAD_GUIDE.md`。


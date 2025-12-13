# PyPI 上传指南

本指南将帮助您将 LNB-MDT 上传到 PyPI。

## 📋 准备工作

### 1. 创建 PyPI 账户

1. 访问 https://pypi.org/account/register/
2. 填写用户名、邮箱和密码
3. 验证邮箱地址
4. （可选）启用双因素认证（2FA）以提高安全性

### 2. 创建 API Token（推荐）

使用 API Token 比使用密码更安全：

1. 登录 PyPI：https://pypi.org/account/login/
2. 进入 **Account settings** → **API tokens**
3. 点击 **Add API token**
4. 设置 Token 名称（如：`lnb-mdt-upload`）
5. 选择作用域：
   - **Entire account**：用于所有项目
   - **Project specific**：仅用于特定项目（选择 `lnb-mdt`）
6. 点击 **Add token**
7. **重要**：复制生成的 token（格式：`pypi-...`），它只会显示一次！

### 3. 安装构建工具

```bash
pip install --upgrade build twine
```

## 🔨 构建分发包

### 方法 1：使用发布脚本（推荐）

```bash
# Linux/macOS
./publish.sh

# Windows
publish.bat
```

### 方法 2：手动构建

```bash
# 1. 清理之前的构建文件
rm -rf build/ dist/ *.egg-info/ LNB_MDT.egg-info/

# 2. 构建分发包
python -m build

# 3. 检查分发包
twine check dist/*
```

构建完成后，`dist/` 目录会包含：
- `lnb-mdt-1.0.0.tar.gz` - 源码分发包
- `lnb_mdt-1.0.0-py3-none-any.whl` - Wheel 分发包

## 🧪 测试上传到 TestPyPI（强烈推荐）

在正式发布前，先上传到 TestPyPI 进行测试：

### 1. 创建 TestPyPI 账户

1. 访问 https://test.pypi.org/account/register/
2. 可以使用与正式 PyPI 相同的用户名和密码

### 2. 上传到 TestPyPI

```bash
twine upload --repository testpypi dist/*
```

系统会提示输入：
- **Username**: `__token__`
- **Password**: 您的 TestPyPI API token（格式：`pypi-...`）

如果没有 API token，可以使用用户名和密码。

### 3. 测试安装

```bash
# 从 TestPyPI 安装
pip install --index-url https://test.pypi.org/simple/ lnb-mdt

# 测试命令
LNB-MDT --help
```

如果测试成功，可以继续上传到正式 PyPI。

## 🚀 上传到正式 PyPI

### 方法 1：使用 API Token（推荐）

```bash
twine upload dist/*
```

输入信息：
- **Username**: `__token__`
- **Password**: 您的 PyPI API token（格式：`pypi-...`）

### 方法 2：使用用户名和密码

```bash
twine upload dist/*
```

输入信息：
- **Username**: 您的 PyPI 用户名
- **Password**: 您的 PyPI 密码

### 方法 3：使用配置文件（可选）

创建 `~/.pypirc` 文件（Linux/macOS）或 `%USERPROFILE%\.pypirc`（Windows）：

```ini
[distutils]
index-servers =
    pypi
    testpypi

[pypi]
username = __token__
password = pypi-您的API_TOKEN

[testpypi]
repository = https://test.pypi.org/legacy/
username = __token__
password = pypi-您的TESTPYPI_API_TOKEN
```

然后直接运行：
```bash
twine upload dist/*
```

## ✅ 验证上传

### 1. 检查 PyPI 页面

等待几分钟后，访问：
https://pypi.org/project/lnb-mdt/

### 2. 测试安装

```bash
# 卸载旧版本（如果已安装）
pip uninstall lnb-mdt -y

# 从 PyPI 安装
pip install lnb-mdt

# 验证安装
LNB-MDT --help
python -c "import LNB_MDT; print('安装成功！')"
```

## 📝 后续版本更新

更新版本时：

1. **更新版本号**：在 `setup.py` 中修改 `version` 字段
   ```python
   version="1.0.1",  # 例如：从 1.0.0 升级到 1.0.1
   ```

2. **重新构建**：
   ```bash
   rm -rf build/ dist/ *.egg-info/
   python -m build
   ```

3. **上传**：
   ```bash
   twine upload dist/*
   ```

## ⚠️ 常见问题

### 1. 包名已存在

如果包名 `lnb-mdt` 已被占用，需要：
- 在 `setup.py` 中修改 `name` 字段
- 选择一个唯一的包名

### 2. 版本号已存在

如果版本号已存在，需要：
- 在 `setup.py` 中更新 `version` 字段
- 使用语义化版本号（如：1.0.0 → 1.0.1）

### 3. 上传失败：认证错误

- 检查 API token 是否正确
- 确保 token 有正确的权限
- 如果使用密码，确保启用了 2FA 时使用应用密码

### 4. 文件太大

PyPI 对文件大小有限制：
- 单个文件最大 60MB
- 总大小建议不超过 100MB

如果文件太大，考虑：
- 移除不必要的文件
- 使用 `.gitignore` 排除大文件
- 检查 `MANIFEST.in` 配置

## 📚 参考资源

- PyPI 官方文档：https://packaging.python.org/en/latest/guides/distributing-packages-using-setuptools/
- TestPyPI：https://test.pypi.org/
- 语义化版本：https://semver.org/

---

**祝您发布顺利！** 🎉


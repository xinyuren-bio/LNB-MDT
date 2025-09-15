from setuptools import setup, find_packages

# 这是一个最小化的 setup.py，仅用于让项目可安装
setup(
    name="LNB-MDT",          # 你的项目名称
    version="0.1.0",        # 任何版本号都可以

    # 自动查找你项目中的所有 Python "包"
    # (它会查找所有包含 __init__.py 文件的文件夹)
    packages=find_packages() 
)
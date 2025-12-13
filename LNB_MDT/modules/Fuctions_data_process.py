from typing import List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

__all__ = ["merge_files_lipid"]

TYPE = {'Area(nm^2)': 0
        , 'Mean Curvature(nm-1)': 0
        , 'Height(nm)': 0
        , 'SZ': 0
        , 'Anisotropy': 1
        , 'Gyration(nm)': 1
        , 'Cluster': 1}


def merge_files_lipid(path: str | List[str], output_path: str, ui_info=None):
    # 按行分割文件路径
    if isinstance(path, str):
        path_list = path.split('\n')
    else:
        path_list = path
    description_set = set()
    merged_df = pd.DataFrame()  # 创建一个空的df来合并多个文件的信息
    column_names = []
    for file in path_list:
        xlsx = pd.ExcelFile(file)
        description = xlsx.parse('Sheet1', nrows=1, header=None).iloc[0, 0]
        description_set.update([description])
        if len(description_set) != 1:
            print("The description of the files is not consistent.")
        df = pd.read_excel(file)
        file_name = file.split('/')[-1].split('.')[0]  # 获取文件名称
        column_names.append(file_name)
        if TYPE[description] == 0:
            results = df.iloc[1:, 5:]
            merged_df[file_name] = results.mean(axis=0)
        elif TYPE[description] == 1:
            results = df.iloc[1:, 1]
            merged_df[file_name] = results
        else:
            print('unsupported type')
    if ui_info and ui_info[1] < ui_info[2]:
        if np.arange(ui_info[1], ui_info[2], ui_info[3]).shape[0] == merged_df.shape[0]:
            merged_df.insert(0, ui_info[0], np.arange(ui_info[1], ui_info[2], ui_info[3]))
        else:
            print(np.arange(ui_info[1], ui_info[2], ui_info[3]).shape[0], '!=', merged_df.shape[0])
            merged_df.insert(0, 'Frames', np.arange(1, merged_df.shape[0] + 1))
        column_names.append(ui_info[0])
    else:
        merged_df.insert(0, 'Frames', np.arange(1, merged_df.shape[0] + 1))
        column_names.append('Frames')

    # 将合并后的 DataFrame 保存为 Excel 文件
    merged_df.to_excel(output_path, index=False, sheet_name='merge')

    print(f"数据已成功合并并保存到 {output_path}")

    for column in column_names[:-1]:
        plt.plot(merged_df[column_names[-1]], merged_df[column])
    plt.show()


if __name__ == '__main__':
    merge_files_lipid(['E:/excel/50.xlsx', 'E:/excel/100.xlsx'], 'E:/excel/merge_50_100.xlsx')
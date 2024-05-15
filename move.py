'''
Author: wds-dxh wdsnpshy@163.com
Date: 2024-05-15 11:09:43
LastEditors: wds-dxh wdsnpshy@163.com
LastEditTime: 2024-05-15 12:16:03
FilePath: \Radar_attitude_detection\move.py
Description: 
微信: 15310638214 
邮箱：wdsnpshy@163.com 
Copyright (c) 2024 by ${wds-dxh}, All Rights Reserved. 
'''
import os
import shutil       # 导入shutil模块, 用于移动文件

# 定义数据集所在的目录和分类后的目标目录
dataset_dir = './images/7 March 2019 West Cumbria Dataset'
target_dir = './dataset/'

# 遍历数据集目录中的所有文件
for filename in os.listdir(dataset_dir):
    if filename.endswith('.jpg') or filename.endswith('.png'):
        # 获取文件名的第一个字母，并构建目标路径
        first_letter = filename[0]
        target_subdir = os.path.join(target_dir, first_letter)
        
        # 如果目标子目录不存在，则创建
        if not os.path.exists(target_subdir):
            os.makedirs(target_subdir)
        
        # 移动文件到目标子目录
        src_path = os.path.join(dataset_dir, filename)
        dest_path = os.path.join(target_subdir, filename)
        shutil.move(src_path, dest_path)

print("分类完成！")

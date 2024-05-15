# Radar_attitude_detection

FMCW雷达人体姿态识别

## 1. 将得到的雷达原始数据转换为多谱勒热图

```matlab
% 获取文件夹中所有.dat文件
% Files = dir(fullfile('D:\nextcloud\Dtat_Set\Radar_attitude_detection\dataset_848\2 March 2017 Dataset\*.dat'));
Files = dir(fullfile('D:\nextcloud\Dtat_Set\Radar_attitude_detection\dataset_848\7 March 2019 West Cumbria Dataset\*.dat'));
lengthFiles = length(Files);
for i=1:lengthFiles
	filename=Files(i).name;
	pathname=Files(i).folder;
	file=[pathname,'/',filename];
	fileID = fopen(file, 'r');
	dataArray = textscan(fileID, '%f');
	fclose(fileID);
	radarData = dataArray{1};
    clearvars fileID dataArray ans;          % 清理不再需要的变量

    % 解析雷达数据参数
    fc = radarData(1);                      % 中心频率
    Tsweep = radarData(2);                  % 扫描时间（毫秒）
    Tsweep = Tsweep / 1000;                 % 转换为秒
    NTS = radarData(3);                     % 每次扫描的采样点数
    Bw = radarData(4);                      % FMCW带宽，对于FSK是频率步进，CW则为0
    Data = radarData(5:end);                % 原始复数数据（实部+虚部）

    
    % 计算辅助参数
    fs = NTS / Tsweep;                      % 采样频率
    record_length = length(Data) / NTS * Tsweep; % 录制时长
    nc = record_length / Tsweep;             % 扫频次数

    %% 数据重塑为扫频单元，并绘制距离-时间图
    Data_time = reshape(Data, [NTS nc]);     % 将数据重塑为扫频单元
    win = ones(NTS, size(Data_time, 2));     % 创建矩形窗函数用于加权

    % FFT处理及IIR滤波
    tmp = fftshift(fft(Data_time .* win), 1); % 对加权后数据进行FFT并中心化
    Data_range = zeros(size(tmp, 1) / 2, size(tmp, 2)); % 初始化距离数据矩阵
    Data_range(1:NTS / 2, :) = tmp(NTS / 2 + 1:NTS, :); % 选取有效范围数据
    ns = oddnumber(size(Data_range, 2)) - 1;            % 确保ns为奇数，用于滤波器设计

    % 动目标显示（MTI）滤波
    [b, a] = butter(4, 0.0075, 'high');           % 设计4阶高通Butterworth滤波器
    [h, f1] = freqz(b, a, ns);                   % 预览滤波器频率响应
    Data_range_MTI = zeros(size(Data_range, 1), ns); % 初始化MTI数据矩阵
    for k = 1:size(Data_range, 1)
        Data_range_MTI(k, 1:ns) = filter(b, a, Data_range(k, 1:ns)); % 应用滤波器
    end

    freq = (0:ns-1) * fs / (2 * ns);            % 频率轴
    range_axis = (freq * 3e8 * Tsweep) / (2 * Bw); % 距离轴转换
    % 裁剪不必要的第一行数据
    Data_range_MTI = Data_range_MTI(2:end, :);
    Data_range = Data_range(2:end, :);

    %% 频谱图处理，以获得多普勒信息
    bin_indl = 10; bin_indu = 30; % 选定分析的距离门范围
    MD.PRF = 1 / Tsweep;        % 脉冲重复频率
    MD.TimeWindowLength = 200;  % 时间窗口长度
    MD.OverlapFactor = 0.95;   % 重叠因子
    MD.OverlapLength = round(MD.TimeWindowLength * MD.OverlapFactor); % 重叠长度
    MD.Pad_Factor = 4;         % 帧填充因子
    MD.FFTPoints = MD.Pad_Factor * MD.TimeWindowLength; % FFT点数
    MD.DopplerBin = MD.PRF / MD.FFTPoints;            % 多普勒分辨率
    MD.DopplerAxis = -MD.PRF / 2 : MD.DopplerBin : MD.PRF / 2 - MD.DopplerBin; % 多普勒轴
    MD.WholeDuration = size(Data_range_MTI, 2) / MD.PRF; % 总持续时间
    MD.NumSegments = floor((size(Data_range_MTI, 2) - MD.TimeWindowLength) / floor(MD.TimeWindowLength * (1 - MD.OverlapFactor))); % 分段数

    % 计算频谱图
    Data_spec_MTI2 = 0; Data_spec2 = 0;
    for RBin = bin_indl:1:bin_indu
        Data_MTI_temp = fftshift(spectrogram(Data_range_MTI(RBin,:), MD.TimeWindowLength, MD.OverlapLength, MD.FFTPoints), 1);
        Data_spec_MTI2 = Data_spec_MTI2 + abs(Data_MTI_temp);
        Data_temp = fftshift(spectrogram(Data_range(RBin,:), MD.TimeWindowLength, MD.OverlapLength, MD.FFTPoints), 1);
        Data_spec2 = Data_spec2 + abs(Data_temp);
    end
    MD.TimeAxis = linspace(0, MD.WholeDuration, size(Data_spec_MTI2, 2)); % 时间轴

    % 反转频谱图以适应MATLAB坐标系
    Data_spec_MTI2 = flipud(Data_spec_MTI2);

    % % 绘制多普勒频谱图
    % figure
    % imagesc(MD.TimeAxis, MD.DopplerAxis .* 3e8 / 2 / 5.8e9, 20*log10(abs(Data_spec_MTI2))); % 绘图，单位转换为速度
    % colormap('jet'); axis xy % 使用彩虹色图并设定坐标轴方向
    % ylim([-6 6]); colorbar % 设置y轴限并添加色标条
    % xlabel('时间[s]', 'FontSize', 16);
    % ylabel('速度 [m/s]','FontSize', 16);
    % set(gca, 'FontSize', 16); % 设置字体大小
    % title(filename) % 图标题为文件名

    out = 20*log10(abs(Data_spec_MTI2));
    imagesc(out); 
    colormap('jet');
    clim = get(gca,'CLim');
    set(gca, 'CLim', clim(2)+[-40,0]);
    [r,c]=size(out);
    set(gca,'xtick',[],'ytick',[]);%去除坐标轴
    set(gca,'LooseInset', get(gca,'TightInset'))    %去除白边
    set(gcf,'innerposition',[0 0 c*16/25 r*16/25])  %调整尺寸
	savepath=['./images/7 March 2019 West Cumbria Dataset/',filename(1:end-4),'.jpg']; %保存地址根据文件名变化
	saveas(gca,savepath,'jpg')%图窗保存
end 


```

## 2. 将得到的图，按照不同的姿态分类，方便后续训练。

```pyhton
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

```

## 3.使用预训练模型resnet训练

```python
import time
import os

from torch.utils.tensorboard import SummaryWriter
from torchvision.models import ResNet18_Weights
from tqdm import tqdm
import pandas as pd
import numpy as np
from torch.utils.data import DataLoader
import torch
import torchvision
import torch.nn as nn
import torch.nn.functional as F
from torchvision import datasets
import matplotlib.pyplot as plt
from torchvision import models
import torch.optim as optim
from torch.optim import lr_scheduler
from torchvision import models
import torch.optim as optim
# import model

# 获取计算硬件
# 有 GPU 就用 GPU，没有就用 CPU
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print('device', device)

from torchvision import transforms
writer = SummaryWriter('logs_train')

# 训练集图像预处理：缩放裁剪、图像增强、转 Tensor、归一化
train_transform = transforms.Compose([transforms.RandomResizedCrop(224),#224*224
                                      transforms.RandomHorizontalFlip(),
                                      transforms.ToTensor(),
                                      transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])
                                     ])

# 测试集图像预处理-RCTN：缩放、裁剪、转 Tensor、归一化
test_transform = transforms.Compose([# 640*480
                                     transforms.RandomResizedCrop(224),
                                     transforms.CenterCrop(224),
                                     transforms.ToTensor(),
                                     transforms.Normalize(
                                         mean=[0.485, 0.456, 0.406],
                                         std=[0.229, 0.224, 0.225])
                                    ])

# 数据集文件夹路径
dataset_dir = './dataset'
train_path = os.path.join(dataset_dir, 'train')
test_path = os.path.join(dataset_dir, 'train')
print('训练集路径', train_path)
print('测试集路径', test_path)
# 载入训练集
train_dataset = datasets.ImageFolder(train_path, train_transform)
# 载入测试集
test_dataset = datasets.ImageFolder(test_path, test_transform)
# 各类别名称
class_names = train_dataset.classes
n_class = len(class_names)
# 映射关系：类别 到 索引号
print(train_dataset.class_to_idx)

# 映射关系：索引号 到 类别
idx_to_labels = {y:x for x,y in train_dataset.class_to_idx.items()}

# 保存为本地的 npy 文件
np.save('idx_to_labels.npy', idx_to_labels)
np.save('labels_to_idx.npy', train_dataset.class_to_idx)

BATCH_SIZE = 1280    

# 训练集的数据加载器
train_loader = DataLoader(train_dataset,
                          batch_size=BATCH_SIZE,
                          shuffle=True,
                          # num_workers=0#线程数
                         )

# 测试集的数据加载器
test_loader = DataLoader(test_dataset,
                         batch_size=BATCH_SIZE,
                         shuffle=False,
                         # num_workers=0
                        )


model = models.resnet18(weights = ResNet18_Weights.IMAGENET1K_V1) # 载入预训练模型
model.fc = nn.Linear(model.fc.in_features, n_class)
# print(model.fc)  # 查看修改后的全连接层
# model = torch.load('checkpoints/mymodel/zzsb.pth', map_location=torch.device('cuda'))
# model = model.CNN_easy(4)
# optimizer = optim.Adam(model.fc.parameters())  # 只优化全连接层的参数
#加载模型训练过的模型，继续训练
# model = torch.load('Radar_attitude_detection.pth')
learning_rate = 5e-3
optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)#优化器
model = model.to(device)
init_image = torch.zeros((1, 3, 224, 224)).to(device)
writer.add_graph(model, init_image)

#启动tensorboard
# tensorboard --logdir=logs_train
#安装tensorboard，pip install tensorboard -i https://pypi.tuna.tsinghua.edu.cn/simple
correct_over = 0
criterion = nn.CrossEntropyLoss() # 交叉熵损失函数
EPOCHS = 500
for epoch in tqdm(range(EPOCHS)):
    model.train()
    for images, labels in train_loader:  # 获取训练集的一个 batch，包含数据和标注
        images = images.to(device)
        labels = labels.to(device)

        outputs = model(images)  # 前向预测，获得当前 batch 的预测结果
        loss = criterion(outputs, labels)  # 比较预测结果和标注，计算当前 batch 的交叉熵损失函数
        writer.add_scalar('Loss/train', loss.item(), epoch) # 记录训练损失
        optimizer.zero_grad()
        loss.backward()  # 损失函数对神经网络权重反向传播求梯度
        optimizer.step()  # 优化更新神经网络权重
    
    # 测试集上的准确率
    model.eval()
    with torch.no_grad():# 不计算梯度
        correct = 0
        total = 0
        for images, labels in test_loader: # 获取测试集  的一个 batch，包含数据和标注
            images = images.to(device)
            labels = labels.to(device)
            outputs = model(images)              # 前向预测，获得当前 batch 的预测置信度
            _, preds = torch.max(outputs, 1)     # 获得最大置信度对应的类别，作为预测结果
            total += labels.size(0)
            correct += (preds == labels).sum()   # 预测正确样本个数
        writer.add_scalar('Accuracy', 100 * correct / total, epoch)  # 记录测试准确率
        print('测试集上的准确率为 {:.3f} %'.format(100 * correct / total))
        if correct > correct_over:
            torch.save(model, './Radar_attitude_detection.pth')
            correct_over = correct
```



## 注：

模型保存为Radar_attitude_detection.pth


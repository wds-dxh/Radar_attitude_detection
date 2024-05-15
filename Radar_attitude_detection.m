%% 数据读取部分
clear         % 清除工作区变量
close all     % 关闭所有图形窗口
[filename, pathname] = uigetfile('*.dat'); % 弹出对话框让用户选择.dat文件
file = [pathname filename];              % 构造完整的文件路径，以便从任意位置读取文件
fileID = fopen(file, 'r');               % 打开文件读取模式
dataArray = textscan(fileID, '%f');      % 读取文件内容为浮点数数组
fclose(fileID);                          % 关闭文件
radarData = dataArray{1};                % 提取雷达数据
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

% 绘制多普勒频谱图
figure
imagesc(MD.TimeAxis, MD.DopplerAxis .* 3e8 / 2 / 5.8e9, 20*log10(abs(Data_spec_MTI2))); % 绘图，单位转换为速度
colormap('jet'); axis xy % 使用彩虹色图并设定坐标轴方向
ylim([-6 6]); colorbar % 设置y轴限并添加色标条
xlabel('时间[s]', 'FontSize', 16);
ylabel('速度 [m/s]','FontSize', 16);
set(gca, 'FontSize', 16); % 设置字体大小
title(filename) % 图标题为文件名



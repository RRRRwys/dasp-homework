# RIR 仿真

## 学习目标

至少掌握一种RIR仿真软件的使用方法

## RIR工具 - RIR-Generator

镜像声源法，是一种由Allen 和 Berkeley[1] 在1979年提出的算法。经常用于各种声学信号处理任务中，来生成房间冲击响应。

RIR-Generator 实现了一个 mex 函数，可以在 MATLAB 中使用，可以使用镜像生源法生成多通道房间冲击响应。

更多的信息可以查看[这里](https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator).

[1] J.B. Allen and D.A. Berkley, "Image method for efficiently simulating small-room acoustics," Journal Acoustic Society of America, 65(4), April 1979, p 943.

## 使用方法

```matlab
c = 340;                     % Sound velocity (m/s)
fs = 16000;                  % Sample frequency (samples/s)
r = [ 2 1.5 2 ];             % Receiver position [ x y z ] (m)
s = [ 2 3.5 2 ];             % Source position [ x y z ] (m)
L = [ 5 4 6 ];               % Room dimensions [ x y z ] (m)
beta = 0 . 4;                % Reverberationtime (s)
nsample = 4096;              % Number of samples
mtype = ’ hypercardioid ’;   % Type of microphone
order = −1;                  % −1 equals maximum reflection order!
dim = 3;                     % Room dimension
orientation = [pi/2 0];      % Microphone orientation [azimuth elevation] in radians
hp_filter = 1;               % Enable high-pass filter

h = rir_generator(c, fs, r, s, L, beta, nsample, mtype, order, dim, orientation, hp_filter);
```

- 输入:
  - c: 声速
  - fs: 采样频率
  - r: 麦克风三维坐标
  - s: 声源坐标
  - L: 房间大小
  - beta: 混响时间 T60 或 1 x 6 向量表示 6 面墙的反射系数
  - nsample: 样本点数(冲击响应长度), 默认是 T60 * fs
  - mtype: 麦克风类型(全向，亚心形，心形，超心形，双向，默认为全向)
  - order: 反射顺序，默认为-1，即最大顺序
  - dim: 房间维数，默认3
  - orientation: 麦克风指向的方向，指定方位角和仰角（以弧度为单位），默认为 [0 0]
  - hp_filter: 即高通滤波器，使用 'false' 禁用高通滤波器，默认启用
- 输出参数:
  - h: M x nsample 矩阵包含计算出的房间脉冲响应
  - beta_hat: 如果混响时间被指定为输入参数，返回对应的反射系数

## 实验设置

- 使用